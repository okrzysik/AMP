// NOTE: This is a modification of the PIO class from SAMRAI

#include <string>

#include "AMP/utils/AMP_MPI.h"
#include "PIO.h"
#include "ParallelBuffer.h"
#include "Utilities.h"

#ifndef NULL
#define NULL 0
#endif

namespace AMP {

int PIO::s_rank                  = -1;
std::ofstream *PIO::s_filestream = nullptr;


/************************************************************************
 * Define the parallel buffers and the associated ostream objects.       *
 ************************************************************************/
static ParallelBuffer pout_buffer;
static ParallelBuffer perr_buffer;
static ParallelBuffer plog_buffer;
static ParallelBuffer pnull_buffer;
std::ostream pout( &pout_buffer );
std::ostream perr( &perr_buffer );
std::ostream plog( &plog_buffer );
std::ostream pnull( &pnull_buffer );


/************************************************************************
 * Initialie the parallel I/O streams.  This routine must be called      *
 * before pout, perr, and plog are used for output but after AMP_MPI     *
 * has been initialized.  By default, logging is disabled.               *
 ************************************************************************/
void PIO::initialize()
{
    AMP_MPI comm = AMP_MPI( AMP_COMM_WORLD );
    s_rank       = comm.getRank();
    s_filestream = nullptr;

    // Initialize the standard parallel output stream
    pout_buffer.setActive( s_rank == 0 );
    pout_buffer.setPrefixString( std::string() );
    pout_buffer.setOutputStream1( &std::cout );
    pout_buffer.setOutputStream2( nullptr );

    // Initialize the error parallel output stream

    std::string buffer = "P=" + Utilities::processorToString( s_rank ) + ":";

    perr_buffer.setActive( true );
    perr_buffer.setPrefixString( buffer );
    perr_buffer.setOutputStream1( &std::cerr );
    perr_buffer.setOutputStream2( nullptr );

    // Initialize the parallel log file (disabled by default)
    plog_buffer.setActive( false );
    plog_buffer.setPrefixString( std::string() );
    plog_buffer.setOutputStream1( nullptr );
    plog_buffer.setOutputStream2( nullptr );
}


/************************************************************************
 * Close the output streams.  Flush both cout and cerr.  If logging,     *
 * then flush and close the log stream.                                  *
 ************************************************************************/
void PIO::finalize()
{
    std::cout.flush();
    std::cerr.flush();
    shutdownFilestream();
    pout_buffer.reset();
    perr_buffer.reset();
    plog_buffer.reset();
}


/************************************************************************
 * If the log file stream is open, then shut down the filestream.        *
 * Close and flush the channel and disconnect the output stream buffers. *
 ************************************************************************/
void PIO::shutdownFilestream()
{
    if ( s_filestream ) {
        s_filestream->flush();
        s_filestream->close();

        delete s_filestream;
        s_filestream = nullptr;

        pout_buffer.setOutputStream2( nullptr );
        perr_buffer.setOutputStream2( nullptr );
        plog_buffer.setOutputStream1( nullptr );
        plog_buffer.setActive( false );
    }
}


/************************************************************************
 * Log messages for node zero only.  If a log stream was open, close     *
 * it.  If this is node zero, then open a new log stream and set the     *
 * appropriate buffer streams to point to the log file.                  *
 ************************************************************************/
void PIO::logOnlyNodeZero( const std::string &filename )
{
    // If the filestream was open, then close it and reset streams
    shutdownFilestream();

    // If this is node zero, then open the log stream and redirect output
    if ( s_rank == 0 ) {
        s_filestream = new std::ofstream( filename.c_str() );
        if ( !( *s_filestream ) ) {
            delete s_filestream;
            s_filestream = nullptr;
            perr << "PIO: Could not open log file ``" << filename.c_str() << "''\n";
        } else {
            pout_buffer.setOutputStream2( s_filestream );
            perr_buffer.setOutputStream2( s_filestream );
            plog_buffer.setOutputStream1( s_filestream );
            plog_buffer.setActive( true );
        }
    }
}


/************************************************************************
 * Log messages for all nodes.  If a log stream was open, the close it.  *
 * Open a log stream on every processor.  The filename for the log file  *
 * will be appended with the processor number.                           *
 ************************************************************************/
void PIO::logAllNodes( const std::string &filename )
{
    /*
     * If the filestream was open, then close it and reset streams
     */

    shutdownFilestream();

    /*
     * Open the log stream and redirect output
     */

    std::string full_filename = filename + "." + Utilities::processorToString( s_rank );
    s_filestream              = new std::ofstream( full_filename.c_str() );

    if ( !( *s_filestream ) ) {
        delete s_filestream;
        s_filestream = nullptr;
        perr << "PIO: Could not open log file ``" << full_filename << "''\n";
    } else {
        pout_buffer.setOutputStream2( s_filestream );
        perr_buffer.setOutputStream2( s_filestream );
        plog_buffer.setOutputStream1( s_filestream );
        plog_buffer.setActive( true );
    }
}


/************************************************************************
 * Suspend logging of data to the file stream.  This does not close the  *
 * filestream (assuming it is open) but just disables logging.           *
 ************************************************************************/
void PIO::suspendLogging()
{
    pout_buffer.setOutputStream2( nullptr );
    perr_buffer.setOutputStream2( nullptr );
    plog_buffer.setOutputStream1( nullptr );
    plog_buffer.setActive( false );
}


/************************************************************************
 * Resume logging of the file stream (assuming it was open).  If the     *
 * file stream is NULL, then do nothing.                                 *
 ************************************************************************/
void PIO::resumeLogging()
{
    if ( s_filestream ) {
        pout_buffer.setOutputStream2( s_filestream );
        perr_buffer.setOutputStream2( s_filestream );
        plog_buffer.setOutputStream1( s_filestream );
        plog_buffer.setActive( true );
    }
}


} // namespace AMP
