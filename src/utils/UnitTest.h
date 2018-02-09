#ifndef included_AMPUnitTest
#define included_AMPUnitTest

#include "AMP/utils/AMP_MPI.h"

#include <sstream>
#include <string>
#include <vector>


namespace AMP {

/*!
 * @brief Class UnitTest is simple utility for running unit tests.
 * It provides basic routines for tracing success or failure of tests,
 * and reporting the results.
 * \par Code Sample:
 * \code
int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;
    try {
        std::cout << "Testing tstOne" << std::endl;
        tstOne(&ut);
    } catch( ... ) {
        ut.failure("An unknown exception was thrown");
    }
    ut.report();
    return ut.NumFail();
}

void tstOne(AMP::UnitTest *ut)
{
    // Run the test code
    if ( problem1 ) {
        ut.failure("Problem 1 detected");
        return
    }
    // Finished running test
    ut.passes("Test XXX passed");
}
 * \endcode

 */
class UnitTest
{
public:
    //! Constructor
    UnitTest();

    //! Destructor
    ~UnitTest();

    //! Indicate a passed test
    inline void passes( const std::string &in )
    {
        pass_messages.push_back( in );
        printUpdate( "passes", in );
    }

    //! Indicate a failed test
    inline void failure( const std::string &in )
    {
        fail_messages.push_back( in );
        printUpdate( "failure", in );
    }

    //! Indicate an expected failed test
    inline void expected_failure( const std::string &in )
    {
        expected_fail_messages.push_back( in );
        printUpdate( "expected_failure", in );
    }

    //! Return the number of passed tests locally
    size_t NumPassLocal() const { return pass_messages.size(); }

    //! Return the number of failed tests locally
    size_t NumFailLocal() const { return fail_messages.size(); }

    //! Return the number of expected failed tests locally
    size_t NumExpectedFailLocal() const { return expected_fail_messages.size(); }

    //! Return the number of passed tests locally
    size_t NumPassGlobal() const;

    //! Return the number of failed tests locally
    size_t NumFailGlobal() const;

    //! Return the number of expected failed tests locally
    size_t NumExpectedFailGlobal() const;

    //! Return the rank of the current processor
    int rank() const;

    //! Return the number of processors
    int size() const;

    /*!
     * Print a report of the passed and failed tests.
     * Note: This is a blocking call that all processors must execute together.
     * Note: Only rank 0 will print the messages (this is necessary as other ranks may not be able
     * to print correctly).
     * @param level     Optional integer specifying the level of reporting (default: 1)
     *                  0: Report the number of tests passed, failed, and expected failures.
     *                  1: Report the number of passed tests (if <=20) or the number passed
     * otherwise,
     *                     report all failures,
     *                     report the number of expected failed tests (if <=50) or the number passed
     * otherwise.
     *                  2: Report all passed, failed, and expected failed tests.
     */
    void report( const int level = 1 ) const;

    //! Clear the messages
    void reset();

    //! Make the unit test operator verbose?
    void verbose( bool verbose = true ) { d_verbose = verbose; }

private:
    std::vector<std::string> pass_messages;
    std::vector<std::string> fail_messages;
    std::vector<std::string> expected_fail_messages;
    AMP_MPI comm;
    bool d_verbose;

private:
    // Make the copy constructor private
    UnitTest( const UnitTest & ) {}

    // Function to pack the messages into a single data stream and send to the given processor
    // Note: This function does not return until the message stream has been sent
    void pack_message_stream( const std::vector<std::string> &messages,
                              const int rank,
                              const int tag ) const;

    // Function to unpack the messages from a single data stream
    // Note: This function does not return until the message stream has been received
    std::vector<std::string> unpack_message_stream( const int rank, const int tag ) const;

    // Print a status update if we are running in verbose mode
    inline void printUpdate( const char *operation, const std::string &msg ) const
    {
        if ( d_verbose ) {
            printf( "UnitTest: %i %s: %s\n", comm.getRank(), operation, msg.c_str() );
        }
    }
};
} // namespace AMP

#endif
