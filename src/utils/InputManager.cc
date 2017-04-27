//
// File:	$URL:
// file:///usr/casc/samrai/repository/AMP/tags/v-2-4-4/source/toolbox/inputdb/InputManager.C $
// Package:	AMP toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	An input manager singleton class that parses input files
//

#include "InputManager.h"
#include "AMPManager.h"
#include "PIO.h"
#include "Parser.h"
#include "ShutdownRegistry.h"
#include "Utilities.h"
#include "utils/AMP_MPI.h"
#include <stdio.h>
#include <stdlib.h>

#ifndef NULL
#define NULL ( 0 )
#endif

namespace AMP {

InputManager *InputManager::s_manager_instance     = nullptr;
bool InputManager::s_registered_callback           = false;
AMP::shared_ptr<Database> InputManager::s_input_db = AMP::shared_ptr<Database>();

/*
*************************************************************************
 *									*
 * Basic singleton classes to create, set, and destroy the manager	*
 * instance.								*
 *									*
*************************************************************************
*/

InputManager *InputManager::getManager()
{
    if ( !s_manager_instance ) {
        s_manager_instance = new InputManager;
    }
    if ( !s_registered_callback ) {
        ShutdownRegistry::registerShutdownRoutine( freeManager,
                                                   ShutdownRegistry::priorityInputManager );
        s_registered_callback = true;
    }
    return ( s_manager_instance );
}

void InputManager::setManager( InputManager *manager )
{
    if ( s_manager_instance ) {
        delete s_manager_instance;
    }
    if ( !s_registered_callback ) {
        ShutdownRegistry::registerShutdownRoutine( freeManager,
                                                   ShutdownRegistry::priorityInputManager );

        s_registered_callback = true;
    }
    s_manager_instance = manager;
}

void InputManager::freeManager()
{
    if ( s_manager_instance )
        delete s_manager_instance;
    s_manager_instance = ( (InputManager *) nullptr );

    s_input_db.reset();
}

/*
*************************************************************************
 *									*
 * The constructor and destructor are protected and call only be called	*
 * by the singleton class or its subclasses.				*
 *									*
*************************************************************************
*/

InputManager::InputManager() : comm(AMP_COMM_WORLD) { }

InputManager::~InputManager() {}

/*
*************************************************************************
 *									*
 * Return whether or not the manager contains an valid input database.   *
 *									*
*************************************************************************
*/

bool InputManager::inputDatabaseExists() { return ( !( s_input_db.get() == nullptr ) ); }

/*
*************************************************************************
 *									*
 * Parse the specified input file and return the new database.		*
 *									*
*************************************************************************
*/

AMP::shared_ptr<InputDatabase> InputManager::parseInputFile( const std::string &filename )
{
    AMP::shared_ptr<InputDatabase> db( new InputDatabase( "main" ) );
    this->parseInputFile( filename, db );
    return ( db );
}


/*
*************************************************************************
 *									*
 * Accessor method for InputManger's root input database.                *
 *									*
*************************************************************************
*/
AMP::shared_ptr<Database> InputManager::getInputDatabase() { return ( s_input_db ); }

/*
*************************************************************************
 *									*
 * Parse the specified input file into the given database.		*
 *									*
*************************************************************************
*/

void InputManager::parseInputFile( const std::string &filename, AMP::shared_ptr<InputDatabase> db )
{
    FILE *fstream = nullptr;
    if ( comm.getRank() == 0 ) {
        fstream = fopen( filename.c_str(), "r" );
    }
    int worked = ( fstream ? 1 : 0 );
    worked     = comm.bcast( worked, 0 );
    if ( !worked ) {
        AMP_ERROR( "InputManager:: Could not open input file``" << filename.c_str() << "''\n" );
    }

    /*
     * Parse input file.
     */
    Parser *parser     = new Parser();
    const int errors   = parser->parse( filename, fstream, db );
    const int warnings = parser->getNumberWarnings();

    if ( errors > 0 ) {
        AMP_WARNING( "InputManager: Errors = " << errors << ", Warnings = " << warnings
                                               << "\n when parsing input file = "
                                               << filename
                                               << std::endl );
        db->printClassData( plog );
        AMP_ERROR( "InputManager exiting..." << std::endl );
    }
    if ( warnings > 0 ) {
        AMP_WARNING( "InputManager: Warnings  = " << warnings << "\n when parsing input file = "
                                                  << filename
                                                  << std::endl );
    }

    /*
     * Store the root database in the static s_input_db variable.
     */
    s_input_db = db;

    delete parser;
    if ( fstream )
        fclose( fstream );
}
}
