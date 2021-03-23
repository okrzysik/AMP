#ifndef included_AMP_PIO
#define included_AMP_PIO

#include <cstdarg>
#include <functional>
#include <iostream>


namespace AMP {


/*!
 * Class ParallelBuffer is a simple I/O stream utility that
 * intercepts output from an ostream and redirects the output as necessary
 * for parallel I/O.  This class defines a stream buffer class for an
 * ostream class.
 */
class ParallelStreamBuffer final : public std::streambuf
{
public:
    enum StreamOutputType : uint8_t {
        null    = 0,
        out     = 1,
        log     = 2,
        err     = 4,
        out_log = 3,
        err_log = 6
    };

    /*!
     * Create a parallel buffer class.  The object will require further
     * initialization to set up the I/O streams and prefix string.
     * @param type      Output type
     */
    explicit ParallelStreamBuffer( StreamOutputType type );

    /*!
     * Set the output stream type
     * @param type      Output type
     */
    void setOutputType( StreamOutputType type );

    /*!
     * Set the output file stream
     * @param stream    Output stream
     */
    void setOutputStream( std::ofstream *stream );

    /*!
     * The destructor simply deallocates any internal data
     * buffers.  It does not modify the output streams.
     */
    virtual ~ParallelStreamBuffer();

    /*!
     * Synchronize the parallel buffer (called from streambuf).
     */
    int sync() override;

    /**
     * Write the specified number of characters into the output stream (called
     * from streambuf).
     */
    std::streamsize xsputn( const char *text, std::streamsize n ) override;

    /*!
     * Write an overflow character into the parallel buffer (called from
     * streambuf).
     */
    int overflow( int ch ) override;

    /*!
     * Read an overflow character from the parallel buffer (called from
     * streambuf).  This is not implemented.  It is needed by the
     * MSVC++ stream implementation.
     */
    int underflow() override;

    /*!
     * Clear the internal buffer's memory
     */
    void reset();

private:
    StreamOutputType d_type;
    size_t d_size;
    size_t d_buffer_size;
    char *d_buffer;
    std::ofstream *d_stream;
    inline void reserve( size_t size );
};


/*!
 * Parallel output stream pout writes to the standard output from node zero
 * only.  Output from other nodes is ignored.  If logging is enabled, then
 * output is mirrored to the log stream, as well.
 */
extern std::ostream pout;

/*!
 * Parallel output stream perr writes to the standard error from all nodes.
 * Output is prepended with the processor number.
 */
extern std::ostream perr;

/*!
 * Parallel output stream plog writes output to the log file.  When logging
 * from multiple processors, the processor number is appended to the filename.
 */
extern std::ostream plog;

/*!
 * Parallel output printp pout writes to the standard output from node zero
 * only.  Output from other nodes is ignored.  If logging is enabled, then
 * output is mirrored to the log stream, as well.
 * The format matches the format for printf
 */
inline int printp( const char *format, ... )
{
    va_list ap;
    va_start( ap, format );
    char tmp[4096];
    int n = vsprintf( tmp, format, ap );
    va_end( ap );
    pout << tmp;
    pout.flush();
    return n;
}


/*!
 * Log messages for node zero only to the specified filename.  All output
 * to pout, perr, and plog on node zero will go to the log file.
 */
void logOnlyNodeZero( const std::string &filename );

/*!
 * Log messages from all nodes.  The diagnostic data for processor XXXXX
 * will be sent to a file with the name filename.XXXXX, where filename is
 * the function argument.
 */
void logAllNodes( const std::string &filename, bool singleStream = false );

/*!
 * Stop logging messages, flush buffers, and reset memory.
 */
void stopLogging();

/*!
 * Redirect the output stream to a user defined function
 */
void overrideCout( std::function<void( const char * )> fun );

/*!
 * Redirect the error stream to a user defined function
 */
void overrideCerr( std::function<void( const char * )> fun );


} // namespace AMP


#endif
