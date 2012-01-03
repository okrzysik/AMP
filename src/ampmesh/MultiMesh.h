#ifndef included_AMP_MultiMesh
#define included_AMP_MultiMesh

#include "ampmesh/Mesh.h"


namespace AMP {
namespace Mesh {



/**
 * \class MultiMesh
 * \brief A concrete mesh class for a multi-mesh
 *
 * \details  This class provides routines for creating and accessing multimeshes.
 *   The concept of a multimesh is a mesh that is composed of multiple meshes.
 *   This takes care of the need for a mesh manager while allowing all operations
 *   on a given mesh to apply to multiple meshes.  Note: all meshes within a multimesh
 *   are stored in a flat array.  This applies when we have a multimesh that may contain
 *   other multimeshes that may (or may not) overlap.
 */
class MultiMesh: public Mesh
{
public:

    /**
     * \param params Parameters for constructing a mesh from an input database
     * \brief Read in mesh files, partition domain, and prepare environment for simulation
     * \details  For trivial parallelsim, this method reads in the meshes on each processor.  Each
     * processor contains a piece of each mesh.  For massive parallelism, each mesh is on its own
     * communicator.  As such, some math libraries must be initialized accordingly.
     */
    MultiMesh ( const MeshParameters::shared_ptr &params );

    //! Deconstructor
     ~MultiMesh ();

    //! Function to copy the mesh (allows use to proply copy the derived class)
    Mesh copy() const;

    /**
     * \brief   Estimate the number of elements in the mesh 
     * \details  This function will estimate the number of elements in the mesh. 
     *   This is used so that we can properly balance the meshes across multiple processors.
     *   Ideally this should be both an accurate estimate and very fast.  It should not require
     *   any communication and should not have to actually load a mesh.
     * \param params Parameters for constructing a mesh from an input database
     */
    static size_t estimateMeshSize( const MeshParameters::shared_ptr &params );

    /* Return the number of local element of the given type
     * \param type   Geometric type
     */
    virtual size_t  numLocalElements( const GeomType type ) const;


    /* Return the global number of elements of the given type.
     * Note: for a multimesh this will require global communication.
     * To avoid this in the future we would need to cache the value, and register
     *  some type of listener to check if the value changed on any sub meshes.
     * \param type   Geometric type
     */
    virtual size_t  numGlobalElements( const GeomType type ) const;


    /* Return the number of ghost elements of the given type on the current processor
     * \param type   Geometric type
     */
    virtual size_t  numGhostElements( const GeomType type, const int gcw ) const;


    /**
     * \brief    Subset a mesh given a MeshID
     * \details  This function will return the mesh with the given meshID.
     *    Note: for multmeshes, this will return the mesh with the given id.
     *    For a single mesh this will return a pointer to itself if the meshID
     *    matches the meshID of the mesh, and a null pointer otherwise.
     * \param meshID  MeshID of the desired mesh
     */
    virtual boost::shared_ptr<Mesh>  Subset( MeshID meshID );


    /**
     * \brief    Subset a mesh given a mesh name
     * \details  This function will return the mesh with the given name.
     *    Note: for multimeshes, this will return the mesh with the given name.
     *    For a single mesh this will return a pointer to itself if the mesh name
     *    matches the name of the mesh, and a null pointer otherwise.
     *    Note: The mesh name is not gaurenteed to be unique.  If there are multiple
     *    meshes with the same name, the first mesh with the given name will be returned.
     *    It is strongly recommended to use the meshID when possible.
     * \param name  Name of the desired mesh
     */
    virtual boost::shared_ptr<Mesh>  Subset ( std::string name );


    /**
     * \brief    Return an MeshIterator over the given geometric objects
     * \details  Return an MeshIterator over the given geometric objects
     * \param type   Geometric type to iterate over
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator getIterator ( const GeomType type, const int gcw=0 );


    /**
     * \brief    Return an MeshIterator over the given geometric objects on the surface
     * \details  Return an MeshIterator over the given geometric objects on the surface
     * \param type   Geometric type to iterate over
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator  getSurfaceIterator ( const GeomType type, const int gcw=0 );


    /**
     * \brief    Return the list of all ID sets in the mesh
     * \details  Return the list of all ID sets in the mesh
     *   Note: for a multimesh this will require global communication.
     *   To avoid this in the future we would need to cache the value, and register
     *   some type of listener to check if the value changed on any sub meshes.
     */
    virtual std::vector<int> getIDSets ( );


    /**
     * \brief    Return an MeshIterator over the given geometric objects on the given ID set
     * \details  Return an MeshIterator over the given geometric objects on the given ID set
     * \param type   Geometric type to iterate over
     * \param id     id for the elements (example: nodeset id)
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator getIDsetIterator ( const GeomType type, const int id, const int gcw=0 );


    /**
     * \brief    Displace the entire mesh
     * \details  This function will displace the entire mesh by a scalar value.
     *   This function is a blocking call for the mesh communicator, and requires
     *   the same value on all processors.  The displacement vector should be the 
     *   size of the physical dimension.
     * \param x  Displacement vector
     */
    virtual void displaceMesh( std::vector<double> x );


    /**
     *  Get the meshIDs of all meshes that compose the current mesh (including its self)
     *  Note: This function will require global communication
     */
    virtual std::vector<MeshID> getAllMeshIDs() const;


    /**
     *  Get the meshIDs of all the basic meshes that compose the current mesh (excluding multimeshes and subset meshes)
     *  Note: This function will require global communication
     */
    virtual std::vector<MeshID> getBaseMeshIDs() const;


private:

    //! Empty constructor for a mesh
    MultiMesh ( ) {};

    //! Function to create the databases for the meshes within the multimesh
    static std::vector<boost::shared_ptr<AMP::Database> >  createDatabases(boost::shared_ptr<AMP::Database> database);

    //! A list of all meshes in the multimesh
    std::vector<AMP::Mesh::Mesh::shared_ptr> d_meshes;

    /**
     * \brief    A function to compute the comms given a weight
     * \details  This function computes the sub communicators given the weights 
     *   for each submesh.  Note that this function requires global communication
     *   on the current comm.
     * \param weights   Standard vector with the relative weights for each group
     *                  Note: the weights do not need to be normalized, 
     *                  but need to be the same on all processors.
     * \param method    Method to use for the load balance calculation
     *                  0: All meshes will be on the same communication
     *                  1: Minimize comm size and split the meshes.
     *                     If there are more processor than meshes, then every
     *                     mesh will be on a different set of processors regardless
     *                     of the effect on the load balancer.
     *                  2: Non-overlapping communicators that tries to achieve a good
     *                     load balancing.  This will attempt to split the meshes onto different
     *                     communicators, but may combine them to achieve a better load
     *                     balance.  This is not implimented yet.
     */
    std::vector<AMP_MPI> loadBalancer( std::vector<double> &weights, int method=1 );


    // Structure used to create communication groups
    struct comm_groups{
        int N_procs;                // Number of processor in the group
        std::vector<int> ids;       // ids in the groups
    };
    
    // Function to distribute N groups with weights onto M processors (M>N) with the greatest number of groups possible (smallest comms)
    std::vector<comm_groups>  independentGroups1( int N_procs, std::vector<std::pair<double,int> >  &ids);

    // Function to distribute N groups with weights onto M processors (N>M) with the greatest number of groups possible (comm size = 1)
    std::vector<comm_groups>  independentGroups2( int N_procs, std::vector<std::pair<double,int> >  &ids );

};

} // Mesh namespace
} // AMP namespace

#endif

