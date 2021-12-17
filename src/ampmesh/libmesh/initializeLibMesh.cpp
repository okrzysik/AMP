#include "AMP/ampmesh/libmesh/initializeLibMesh.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UtilityMacros.h"

#include <cstring>

// LibMesh include
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"

// Petsc include (needed to fix PETSC_COMM_WORLD problem with libmesh)
#ifdef USE_EXT_PETSC
    #include "petscsys.h"
#endif

namespace AMP {
namespace Mesh {


// Initialize static member variables
volatile int initializeLibMesh::N_copies = 0;
void *initializeLibMesh::lminit          = nullptr;
AMP_MPI initializeLibMesh::d_comm        = AMP_MPI( AMP_COMM_NULL );


/************************************************************
 * Class to wrap libmesh MPI types to properly free datatype *
 ************************************************************/
/*template<class TYPE>
class libMeshWrapperType: public libMesh::Parallel::StandardType<TYPE>
{
public:
    libMeshWrapperType( ): libMesh::Parallel::StandardType<TYPE>() {}
    ~libMeshWrapperType() { free(); }
}
std::shared_ptr<libMeshWrapperType<Hilbert::HilbertIndices> > type_hilbert;*/


/************************************************************
 * Function to alter the command line arguments for libmesh  *
 ************************************************************/
static int add_libmesh_cmdline( const int argc, const char **argv, char ***argv_new )
{
    const int N_add = 0; // Number of additional arguments we want to add
    // Copy the existing command-line arguments (shifting by the number of additional arguments)
    *argv_new = new char *[argc + N_add];
    for ( int i = 0; i < argc; i++ ) {
        ( *argv_new )[i] = new char[strlen( argv[i] ) + 1];
        strcpy( ( *argv_new )[i], argv[i] );
    }
    /*// Add command to keep cout from all processors (not just rank 0)
    (*argv_new)[argc] = new char [12];
    strcpy ( (*argv_new)[argc] , "--keep-cout" );*/
    return argc + N_add;
}


/************************************************************
 * Constructor initilize libmesh on the given comm           *
 ************************************************************/
initializeLibMesh::initializeLibMesh( const AMP_MPI &comm )
{
    if ( N_copies > 0 ) {
        // libmesh is already initialized, check if it is compatible with the current comm
        bool test = canBeInitialized( comm );
        if ( test ) {
            // Add 1 to the count and return
            N_copies++;
            return;
        } else {
            // We can't initialize libmesh
            AMP_ERROR( "libmesh was previously initialized with a different (incompatible) comm" );
        }
    } else {
        // libmesh is not initialized
        if ( lminit != nullptr )
            AMP_ERROR( "Internal error" );
        // Use a barrier to ensure all processors are at the same point
        N_copies = 1;
        d_comm   = comm.dup(); // Create a seperate duplicate comm for libmesh
        d_comm.barrier();
        // Reinitialize LibMesh with the new communicator
        int argc_libmesh    = 0;
        char **argv_libmesh = nullptr;
        const int argc      = AMPManager::get_argc();
        const auto **argv   = (const char **) AMPManager::get_argv();
        argc_libmesh        = add_libmesh_cmdline( argc, argv, &argv_libmesh );
#ifdef USE_EXT_MPI
    #ifdef USE_EXT_PETSC
        MPI_Comm petsc_comm = PETSC_COMM_WORLD;
    #endif
        lminit = new libMesh::LibMeshInit( argc_libmesh, argv_libmesh, d_comm.getCommunicator() );
    #ifdef USE_EXT_PETSC
        PETSC_COMM_WORLD = petsc_comm;
    #endif
#else
        lminit = new libMesh::LibMeshInit( argc_libmesh, argv_libmesh );
#endif
        for ( int i = 0; i < argc_libmesh; i++ )
            delete[] argv_libmesh[i];
        delete[] argv_libmesh;
        // Initialize libmesh MPI types so we can safely free them
        // type_hilbert.reset( new libMeshWrapperType<Hilbert::HilbertIndices>() );
        // Reset the error handlers
        AMP::AMPManager::setHandlers();
    }
}


/************************************************************
 * Deconstructor that will finalize libmesh                  *
 ************************************************************/
initializeLibMesh::~initializeLibMesh()
{
    if ( N_copies <= 0 )
        AMP_ERROR( "Internal error" );
    // Use a barrier to ensure all processors are at the same point
    d_comm.barrier();
    if ( N_copies == 1 ) {
        // Shutdown libmesh
        if ( lminit == nullptr )
            AMP_ERROR( "Internal error" );
        // Free libmesh MPI types
        // type_hilbert.reset();
        // Delete libmesh
        delete (libMesh::LibMeshInit *) lminit;
        lminit   = nullptr;
        d_comm   = AMP_MPI( AMP_COMM_NULL );
        N_copies = 0;
    } else {
        N_copies--;
    }
}


/************************************************************
 * Function check if initiallize can be called successfully  *
 ************************************************************/
bool initializeLibMesh::canBeInitialized( AMP_MPI comm )
{
    if ( N_copies == 0 )
        return true;
    if ( comm == d_comm )
        return true;
    if ( d_comm.compare( comm ) != 0 )
        return true;
    return false;
}


/************************************************************
 * Function to check if libmesh has been initialized         *
 ************************************************************/
bool initializeLibMesh::isInitialized() { return N_copies > 0; }


} // namespace Mesh
} // namespace AMP
