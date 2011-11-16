#ifndef included_tbox_AMPManager
#define included_tbox_AMPManager

#include "utils/AMP_MPI.h"


namespace AMP {

/*!
 * @brief Class AMPManagerProperties is a class that contains the various startup options for AMP
 */
class AMPManagerProperties
{
public:

    //! Empty constructor.  This creates the default startup parameters.
    AMPManagerProperties();

    /*!
     * Use MPI_Abort.  This determines if we want to use MPI_Abort (default) or
     * throw an exception when we abort.  Note that using MPI_Abort will safely terminate
     * all processes, but is incompatible with try, catch.  A standard exception can be
     * caught, but can easily lead to errors where different processors get stuck when 
     * a single process encounters an error.  It is strongly recommended to use MPI_Abort.
     */
    bool use_MPI_Abort;

    //! Print the time required to initialize or shutdown each package.  Default is false.
    bool print_times;

private:
    friend class AMPManager;
    
};



/*!
 * @brief Class AMPManager is a utility for managing startup and shutdown for AMP 
 * applications and for changing the maximum number of patch data components supported 
 * by AMP patches.  All applications should call AMPManager::startup() at the 
 * beginning of the program.  Startup will initialize all packages including MPI.
 * AMPManager::shutdown() should be called at the end of the program, before calling exit(0). 
 * Note that the shutdown, but it shuts down the packages, MPI, and deallocates memory.
 */
class AMPManager
{
public:
    /*!
     * Initialize the AMP package.  Depending on the architecture and
     * compile flags, this routine sets up MPI, initializes IEEE exception
     * handlers, and other architecture-specific details.
     */
    static void startup(int argc, char *argv[], const AMPManagerProperties &properties=AMPManagerProperties());

    /*!
     * Shutdown the AMP package.  Depending on the compile flags set at
     * compile-time, this routine shuts down MPI and calls registered shutdown
     * handlers.
     */
    static void shutdown();

    /*!
     * Function to check if AMP has been initialized
     */
    static bool isInitialized() { return initialized; }

    /*!
     * Function to return the number command line arguments that were used to initialize AMP.
     */
    static int get_argc( );

    /*!
     * Function to return the command line arguments that were used to initialize AMP.
     * Note: This returns the pointer address for the command line arguments.  The user
     * is responsible to ensure that the arguments are not modified.
     */
    static char** get_argv( );

    /*!
     * Function to return the AMPManagerProperties that was used to initialize AMP
     */
    static AMPManagerProperties getAMPManagerProperties( );

private:
    // Private constructor (we do not actually want to create an object)
    AMPManager() {  }

    // Static variables
    static bool initialized;
    static bool called_MPI_Init;
    static bool called_PetscInitialize;
    static bool use_MPI_Abort;
    static bool print_times;
    static AMP_MPI comm;
    static int argc;
    static char** argv;
    static AMPManagerProperties properties;

    //! abort must be a friend to access use_MPI_Abort to change the abort behavior
    friend void AMP::Utilities::abort(const std::string &, const std::string &, const int);

};


}

#endif
