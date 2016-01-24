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

    //! The default level for the MPI timers.  Default is 2
    int profile_MPI_level;

    /*!
     *  MPI communicator to use for AMP_COMM_WORLD.  By default this should be set to
     *  AMP_COMM_WORLD if MPI is not initialized.  If MPI is initialized, this can be
     *  set to AMP_COMM_WORLD, MPI_COMM_WORLD, or any valid MPI communicator.  If it
     *  is set to AMP_COMM_WORLD then is will default internally to MPI_COMM_WORLD.
     *  Note:  Currently AMP_COMM_WORLD cannot be changed once AMPManager::startup
     *  has been called.
     */
    MPI_Comm COMM_WORLD;

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
    static void startup( int argc,
                         char *argv[],
                         const AMPManagerProperties &properties = AMPManagerProperties() );

    /*!
     * Shutdown the AMP package.  Depending on the compile flags set at
     * compile-time, this routine shuts down MPI and calls registered shutdown
     * handlers.
     */
    static void shutdown();

    /*!
     * Function to check if AMP has been initialized
     */
    static bool isInitialized() { return initialized != 0; }

    /*!
     * Function to return the number command line arguments that were used to initialize AMP.
     */
    static int get_argc();

    /*!
     * Function to return the command line arguments that were used to initialize AMP.
     * Note: This returns the pointer address for the command line arguments.  The user
     * is responsible to ensure that the arguments are not modified.
     */
    static char **get_argv();

    /*!
     * Function to return the AMPManagerProperties that was used to initialize AMP
     */
    static AMPManagerProperties getAMPManagerProperties();

    //! Static function to terminate AMP
    static void terminate_AMP( std::string message );

private:
    // Private constructor (we do not actually want to create an object)
    AMPManager() {}

    // Static variables
    static int initialized;
    static int rank;
    static bool called_MPI_Init;
    static bool called_PetscInitialize;
    static bool use_MPI_Abort;
    static bool print_times;
    static AMP_MPI comm_world;
    static int argc;
    static char **argv;
    static AMPManagerProperties properties;
#ifdef USE_EXT_MPI
    static AMP::shared_ptr<MPI_Errhandler> mpierr;
#endif

    //! abort must be a friend to access use_MPI_Abort to change the abort behavior
    friend void AMP::Utilities::abort( const std::string &, const std::string &, const int );

    //! AMP_MPI must be a friend to access comm_world and the MPI error handler
    friend class AMP::AMP_MPI;

    //! Function to create the arguments to pass to petsc
    static std::vector<char *> getPetscArgs();

    //! Functions to initialize/destroy the mpi error handler
    static void setMPIErrorHandler();
    static void clearMPIErrorHandler();
    static void exitFun();
};
}

#endif
