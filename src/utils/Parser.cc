//
// File:	$URL:
// file:///usr/casc/samrai/repository/AMP/tags/v-2-4-4/source/toolbox/inputdb/Parser.C
// $
// Package:	AMP toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Parser that reads the input database grammar
//

#include "Parser.h"
#include "PIO.h"

#ifdef DEBUG_NO_INLINE
#include "Parser.I"
#endif

#ifndef NULL
#define NULL ( 0 )
#endif

namespace AMP {

extern int yyparse();
extern void yyrestart( FILE * );

extern void parser_static_table_initialize();

Parser *Parser::s_default_parser         = nullptr;
bool Parser::s_static_tables_initialized = 0;

/*
*************************************************************************
*									*
* The constructor creates an unitialized parser object.  All of the	*
* interesting work is done by member function parse().			*
*									*
*************************************************************************
*/

Parser::Parser()
{
    if ( !s_static_tables_initialized ) {
        parser_static_table_initialize();
        s_static_tables_initialized = 1;
    }
    comm = AMP_MPI( AMP_COMM_WORLD );
}

/*
*************************************************************************
*									*
* The destructor automatically deallocates the parser object data.	*
*									*
*************************************************************************
*/

Parser::~Parser() {}

/*
*************************************************************************
*									*
* Begin parsing the input database file.  Return the number of errors	*
* encountered in the parse.						*
*									*
*************************************************************************
*/

int Parser::parse( const std::string &filename, FILE *fstream, AMP::shared_ptr<Database> database )
{
    d_errors   = 0;
    d_warnings = 0;

    // Find the path in the filename, if one exists
    std::string::size_type slash_pos = filename.find_last_of( '/' );
    if ( slash_pos == std::string::npos ) {
        d_pathname = "";
    } else {
        d_pathname = filename.substr( 0, slash_pos + 1 );
    }

    ParseData pd;
    pd.d_filename   = filename;
    pd.d_fstream    = fstream;
    pd.d_linenumber = 1;
    pd.d_cursor     = 1;
    pd.d_nextcursor = 1;
    d_parse_stack.clear();
    d_parse_stack.push_front( pd );

    d_scope_stack.clear();
    d_scope_stack.push_front( database );

    s_default_parser = this;
    yyrestart( nullptr );
    if ( yyparse() && ( d_errors == 0 ) ) {
        error( "Unexpected parse error" );
    }
    s_default_parser = nullptr;

    d_parse_stack.clear();
    d_scope_stack.clear();

    return ( d_errors );
}

/*
*************************************************************************
*									*
* Advance the cursor to the next line in the current input file.	*
*									*
*************************************************************************
*/

void Parser::advanceLine( const int nline )
{
    Parser::ParseData &pd = d_parse_stack.front();
    pd.d_linenumber += nline;
    pd.d_cursor     = 1;
    pd.d_nextcursor = 1;
}

/*
*************************************************************************
*									*
* Advance the cursor position by the token in the specified string.	*
* Tabs are expanded assuming tab stops at eight character markers.	*
*									*
*************************************************************************
*/

void Parser::advanceCursor( const std::string &token )
{
    Parser::ParseData &pd = d_parse_stack.front();
    pd.d_cursor           = pd.d_nextcursor;
    for ( const auto &elem : token ) {
        if ( elem == '\t' ) {
            pd.d_nextcursor = ( ( pd.d_nextcursor + 7 ) & ( ~7 ) ) + 1;
        } else {
            pd.d_nextcursor++;
        }
    }
}

/*
*************************************************************************
*									*
* Print out errors to pout and track the number of errors.		*
*									*
*************************************************************************
*/

void Parser::error( const std::string &message )
{
    Parser::ParseData &pd = d_parse_stack.front();

    pout << "Error in " << pd.d_filename << " at line " << pd.d_linenumber << " column "
         << pd.d_cursor << " : " << message << std::endl
         << std::flush;

    pout << pd.d_linebuffer << std::endl << std::flush;

    for ( int i = 0; i < pd.d_cursor; i++ )
        pout << " ";
    pout << "^\n";

    d_errors++;
}

/*
*************************************************************************
*									*
* Print out warnings to pout and track the number of warnings.		*
*									*
*************************************************************************
*/

void Parser::warning( const std::string &message )
{
    Parser::ParseData &pd = d_parse_stack.front();

    pout << "Warning in " << pd.d_filename << " at line " << pd.d_linenumber << " column "
         << pd.d_cursor << " : " << message << std::endl
         << std::flush;

    pout << pd.d_linebuffer << std::endl << std::flush;

    for ( int i = 0; i < pd.d_cursor; i++ )
        pout << " ";
    pout << "^\n";

    d_warnings++;
}

/*
*************************************************************************
*									*
* Set the input line which is currently being parsed.    		*
*									*
*************************************************************************
*/

void Parser::setLine( const std::string &line )
{
    Parser::ParseData &pd = d_parse_stack.front();
    pd.d_linebuffer       = line;
}

/*
*************************************************************************
*									*
* Iterate through the database scopes, looking for the first match on	*
* the key value.							*
*									*
*************************************************************************
*/

AMP::shared_ptr<Database> Parser::getDatabaseWithKey( const std::string &key )
{
    AMP::shared_ptr<Database> returnPtr;

    std::list<AMP::shared_ptr<Database>>::iterator i = d_scope_stack.begin();
    for ( ; i != d_scope_stack.end(); i++ ) {
        if ( ( *i )->keyExists( key ) )
            return ( ( *i ) );
    }
    return ( returnPtr );
}

/*
*************************************************************************
*									*
* Create a new parse state on the parse stack and open the specified	*
* new file for reading.							*
*									*
*************************************************************************
*/

bool Parser::pushIncludeFile( const std::string &filename )
{
    FILE *fstream = nullptr;

    std::string filename_with_path;

    // If this is not a fully qualified pathname use
    // current search path
    std::string::size_type slash_pos;
    slash_pos = filename.find_first_of( '/' );
    if ( slash_pos == 0 ) {
        filename_with_path = filename;
    } else {
        filename_with_path = d_pathname;
        filename_with_path += filename;
    }

    if ( comm.getRank() == 0 ) {
        fstream = fopen( filename_with_path.c_str(), "r" );
    }

    int worked = ( fstream ? 1 : 0 );
    worked     = comm.bcast( worked, 0 );

    if ( !worked ) {
        error( "Could not open include file ``" + filename_with_path + "''" );
    } else {
        ParseData pd;
        pd.d_filename   = filename_with_path;
        pd.d_fstream    = fstream;
        pd.d_linenumber = 1;
        pd.d_cursor     = 1;
        pd.d_nextcursor = 1;
        d_parse_stack.push_front( pd );
    }

    return ( worked ? true : false );
}

/*
*************************************************************************
*									*
* Close the current input file and pop the parse stack.			*
*									*
*************************************************************************
*/

void Parser::popIncludeFile()
{
    Parser::ParseData &pd = d_parse_stack.front();
    if ( pd.d_fstream )
        fclose( pd.d_fstream );
    d_parse_stack.pop_front();
}

/*
*************************************************************************
*									*
* Manage the input reading for the flex scanner.  If running with MPI,	*
* the node zero reads the data and broadcasts the length and the data	*
* to all processors.							*
*									*
*************************************************************************
*/

int Parser::yyinput( char *buffer, const int max_size )
{
    int byte = 0;
    if ( comm.getRank() == 0 ) {
        byte = fread( buffer, 1, max_size, d_parse_stack.front().d_fstream );
    }
    byte = comm.bcast( byte, 0 );
    if ( byte > 0 ) {
        comm.bcast( buffer, byte, 0 );
    }
    return ( byte );
}


} // namespace AMP
