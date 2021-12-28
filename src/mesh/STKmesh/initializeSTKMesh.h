#ifndef included_AMP_initializeSTKMesh
#define included_AMP_initializeSTKMesh

#include "AMP/utils/AMP_MPI.h"

namespace AMP::Mesh {


/**
 * \class initializeSTKMesh
 * \brief A class to initialize STKMesh
 * \details  This class provides routines for initializing and finalizing STKmesh.
 *    Note: this class is not thread safe yet.
 */
class initializeSTKMesh
{
public:
    /*!
     *  Constructor that inializes STKmesh on the given communicator.
     *  Note: STKmesh can only be initialized on one comm at a given time.
     *  This function can be called more than once.  If the the comms share
     *  the same comm, then this function will succeed.  It the comms are
     *  different, then this function will throw an error.  If a different comm
     *  is needed, then the user needs to call finalize first.
     */
    initializeSTKMesh( AMP_MPI comm );

    /*!
     *  Deconstructor that finalizes STKmesh.  This allows STKmesh to be reinitialized on
     *  a new comm, and frees any memory in use by STKmesh.
     */
    virtual ~initializeSTKMesh();

    /*!
     *  Function to check if STKmesh can be initialized with the given comm.
     *  If this function returns true, then initializeSTKMesh::initialize
     *  can be safely called with the given comm.  Otherwise, initializeSTKMesh::finalize
     * must be called first.
     */
    static bool canBeInitialized( AMP_MPI comm );

    /*!
     *  Function to check if STKmesh has been initialized for any communicator
     *  (are there any copies of the class that have not been destroyed)
     */
    static bool isInitialized();

private:
    // Constructor
    initializeSTKMesh(){};
};


} // namespace AMP::Mesh

#endif
