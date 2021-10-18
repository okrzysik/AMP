#include "AMP/utils/PIO.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Utilities.h"

#include <cstring>
#include <fstream>
#include <string>


namespace AMP {


static ParallelStreamBuffer pout_buffer( ParallelStreamBuffer::StreamOutputType::out );
static ParallelStreamBuffer perr_buffer( ParallelStreamBuffer::StreamOutputType::err );
static ParallelStreamBuffer plog_buffer( ParallelStreamBuffer::StreamOutputType::null );


std::ostream pout( &pout_buffer );
std::ostream perr( &perr_buffer );
std::ostream plog( &plog_buffer );


std::function<void( const char * )> cout_function;
std::function<void( const char * )> cerr_function;


// Get the global rank
inline int getRank()
{
    static int rank = -1;
    if ( rank == -1 ) {
        if ( AMP_MPI::MPI_Active() )
            rank = AMP_MPI( AMP_COMM_WORLD ).getRank();
    }
    return rank;
}


/****************************************************************************
 *  Functions to control logging                                             *
 ****************************************************************************/
std::ofstream *global_filestream = nullptr;
static void shutdownFilestream()
{
    if ( global_filestream != nullptr ) {
        global_filestream->flush();
        global_filestream->close();
        delete global_filestream;
        global_filestream = nullptr;
    }
}
void logOnlyNodeZero( const std::string &filename )
{
    if ( getRank() == 0 )
        logAllNodes( filename, true );
}
void logAllNodes( const std::string &filename, bool singleStream )
{
    // If the filestream was open, then close it and reset streams
    shutdownFilestream();

    // Open the log stream and redirect output
    std::string full_filename = filename;
    if ( !singleStream ) {
        full_filename += AMP::Utilities::stringf( ".%04i", getRank() );
    }
    global_filestream = new std::ofstream( full_filename.c_str() );

    if ( !( *global_filestream ) ) {
        delete global_filestream;
        global_filestream = nullptr;
        perr << "PIO: Could not open log file ``" << full_filename << "''\n";
    } else {
        pout_buffer.setOutputStream( global_filestream );
        perr_buffer.setOutputStream( global_filestream );
        plog_buffer.setOutputStream( global_filestream );
    }
}


/****************************************************************************
 *  ParallelStreamBuffer class                                               *
 ****************************************************************************/
void stopLogging()
{
    std::cout.flush();
    std::cerr.flush();
    pout_buffer.reset();
    perr_buffer.reset();
    plog_buffer.reset();
    shutdownFilestream();
    delete global_filestream;
    global_filestream = nullptr;
}


/****************************************************************************
 *  ParallelStreamBuffer class                                               *
 ****************************************************************************/
ParallelStreamBuffer::ParallelStreamBuffer( StreamOutputType type )
    : d_type( StreamOutputType::null ),
      d_size( 0 ),
      d_buffer_size( 0 ),
      d_buffer( nullptr ),
      d_stream( nullptr )
{
    setOutputType( type );
}
ParallelStreamBuffer::~ParallelStreamBuffer() { delete[] d_buffer; }
void ParallelStreamBuffer::setOutputType( StreamOutputType type ) { d_type = type; }
void ParallelStreamBuffer::setOutputStream( std::ofstream *stream )
{
    setOutputType( static_cast<StreamOutputType>( d_type | StreamOutputType::log ) );
    d_stream = stream;
}
int ParallelStreamBuffer::sync()
{
    if ( d_buffer == nullptr )
        return 0;
    if ( ( d_type & StreamOutputType::out ) && getRank() <= 0 ) {
        // Write to std::cout
        if ( cout_function )
            cout_function( d_buffer );
        else
            std::cout << d_buffer;
    }
    if ( d_type & StreamOutputType::err ) {
        // Write to std::cerr
        if ( cerr_function )
            cerr_function( d_buffer );
        else
            std::cerr << d_buffer;
    }
    if ( d_type & StreamOutputType::log ) {
        // Write to log file
        *d_stream << d_buffer;
        d_stream->flush();
    }
    d_size = 0;
    memset( d_buffer, 0, d_buffer_size );
    return 0;
}
void ParallelStreamBuffer::reserve( size_t size )
{
    if ( size > d_buffer_size ) {
        // Determine the new buffer size
        d_buffer_size = std::max<size_t>( d_buffer_size, 1024 );
        while ( size > d_buffer_size )
            d_buffer_size *= 2;
        // Allocate the new buffer
        char *tmp = d_buffer;
        d_buffer  = new char[d_buffer_size];
        memset( d_buffer, 0, d_buffer_size );
        memcpy( d_buffer, tmp, d_size );
        delete[] tmp;
    }
}
std::streamsize ParallelStreamBuffer::xsputn( const char *text, std::streamsize n )
{
    reserve( d_size + n );
    memcpy( &d_buffer[d_size], text, n );
    d_size += n;
    if ( text[n - 1] == 0 || text[n - 1] == 10 )
        sync();
    return n;
}
int ParallelStreamBuffer::overflow( int ch )
{
    reserve( d_size + 1 );
    d_buffer[d_size] = ch;
    d_size++;
    if ( ch == 0 || ch == 10 )
        sync();
    return std::char_traits<char>::to_int_type( ch );
}
int ParallelStreamBuffer::underflow() { return -1; }
void ParallelStreamBuffer::reset()
{
    sync();
    delete[] d_buffer;
    d_buffer      = nullptr;
    d_buffer_size = 0;
}


/****************************************************************************
 *  Redirect the output/error stream to a user defined function              *
 ****************************************************************************/
void overrideCout( std::function<void( const char * )> fun ) { cout_function = fun; }
void overrideCerr( std::function<void( const char * )> fun ) { cerr_function = fun; }


} // namespace AMP
