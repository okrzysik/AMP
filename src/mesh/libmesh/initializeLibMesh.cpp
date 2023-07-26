#include "AMP/mesh/libmesh/initializeLibMesh.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UtilityMacros.h"

#include <cstring>

// LibMesh include
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"

// Petsc include (needed to fix PETSC_COMM_WORLD problem with libmesh)
#ifdef AMP_USE_PETSC
    #include "petscsys.h"
#endif

namespace AMP::Mesh {


// Initialize static member variables
volatile int initializeLibMesh::N_copies = 0;
void *initializeLibMesh::lminit          = nullptr;
AMP_MPI initializeLibMesh::d_comm        = AMP_MPI( AMP_COMM_NULL );


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
        auto args              = AMPManager::get_argv();
        char disableRefCount[] = "--disable-refcount-printing";
        args.push_back( disableRefCount );
#ifdef AMP_USE_MPI
    #ifdef AMP_USE_PETSC
        MPI_Comm petsc_comm = PETSC_COMM_WORLD;
    #endif
        lminit = new libMesh::LibMeshInit( args.size(), args.data(), d_comm.getCommunicator() );
    #ifdef AMP_USE_PETSC
        PETSC_COMM_WORLD = petsc_comm;
    #endif
#else
        lminit = new libMesh::LibMeshInit( args.size(), args.data() );
#endif
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
bool initializeLibMesh::canBeInitialized( const AMP_MPI &comm )
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


} // namespace AMP::Mesh
