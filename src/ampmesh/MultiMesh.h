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
     * \brief Default constructor
     * \details  This constructor works with the input parameters to create the mesh.
     * For trivial parallelsim, this method reads in the meshes on each processor.  Each
     * processor contains a piece of each mesh.  For massive parallelism, each mesh is on its own
     * communicator.  As such, some math libraries must be initialized accordingly.
     * \param params Parameters for constructing a mesh from an input database
     */
    MultiMesh ( const MeshParameters::shared_ptr &params );


    /**
     * \brief Contructor to create a MultiMesh from existing meshes
     * \details  This constructor takes a list of meshes and a communicator
     *    and generates the appropriate multimesh
     * \param comm      Desired communicator for the multimesh
     * \param meshes    Meshes to be used as part of the multimesh
     */
    MultiMesh ( const AMP_MPI &comm, const std::vector<Mesh::shared_ptr> &meshes );


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
    virtual boost::shared_ptr<Mesh>  Subset( MeshID meshID ) const;


    /**
     * \brief    Subset a mesh given a MeshIterator
     * \details  This function will subset a mesh over a given iterator.
     *   This will return a new mesh object.
     * \param iterator  MeshIterator used to subset
     */
    virtual boost::shared_ptr<Mesh>  Subset ( const MeshIterator &iterator ) const;


    /**
     * \brief    Subset a mesh given a mesh name
     * \details  This function will return the mesh with the given name.
     *    For a single mesh this will return a pointer to itself if the mesh name
     *    matches the name of the mesh, and a null pointer otherwise.
     *    Note: The mesh name is not gaurenteed to be unique.  If there are multiple
     *    meshes with the same name, all meshed with the given name will be returned 
     *    within a new multimesh.
     *    It is strongly recommended to use the meshID when possible.
     * \param name  Name of the desired mesh
     */
    virtual boost::shared_ptr<Mesh>  Subset ( std::string name ) const;


    /**
     * \brief    Return an MeshIterator over the given geometric objects
     * \details  Return an MeshIterator over the given geometric objects
     * \param type   Geometric type to iterate over
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator getIterator ( const GeomType type, const int gcw=0 ) const;


    /**
     * \brief    Return an MeshIterator over the given geometric objects on the surface
     * \details  Return an MeshIterator over the given geometric objects on the surface
     * \param type   Geometric type to iterate over
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator  getSurfaceIterator ( const GeomType type, const int gcw=0 ) const;


    /**
     * \brief    Return the list of all boundary ID sets in the mesh
     * \details  Return the list of all boundary ID sets in the mesh
     * Note: depending on the mesh this routine may require global communication across the mesh.
     */
    virtual std::vector<int> getBoundaryIDs ( ) const;


    /**
     * \brief    Return an MeshIterator over the given geometric objects on the given boundary ID set
     * \details  Return an MeshIterator over the given geometric objects on the given boundary ID set
     * \param type   Geometric type to iterate over
     * \param id     Boundary id for the elements (example: sideset id)
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator getBoundaryIDIterator ( const GeomType type, const int id, const int gcw=0 ) const;

    /**
     * \brief    Return the list of all boundary ID sets in the mesh
     * \details  Return the list of all boundary ID sets in the mesh
     * Note: depending on the mesh this routine may require global communication across the mesh.
     */
    virtual std::vector<int> getBlockIDs ( ) const;


    /**
     * \brief    Return an MeshIterator over the given geometric objects on the given block ID set
     * \details  Return an MeshIterator over the given geometric objects on the given block ID set
     * \param type   Geometric type to iterate over
     * \param id     Block id for the elements (example: block id in cubit, subdomain in libmesh)
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator getBlockIDIterator ( const GeomType type, const int id, const int gcw=0 ) const;


    /**
     * \brief    Displace the entire mesh
     * \details  This function will displace the entire mesh by a scalar value.
     *   This function is a blocking call for the mesh communicator, and requires
     *   the same value on all processors.  The displacement vector should be the 
     *   size of the physical dimension.
     * \param x  Displacement vector
     */
    virtual void displaceMesh( std::vector<double> x );


    //! Is the current mesh a base mesh
    virtual inline bool isBaseMesh() const { return false; }


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


    /**
     *  Get the meshIDs of all meshes that compose the current mesh (including its self)
     *  on the current processor.
     */
    virtual std::vector<MeshID> getLocalMeshIDs() const;


    /**
     *  Get the meshIDs of all the basic meshes that compose the current mesh 
     *  (excluding multimeshes and subset meshes) on the current processor.
     */
    virtual std::vector<MeshID> getLocalBaseMeshIDs() const;


    /**
     *  Get the meshes composing the multimesh
     */
    virtual std::vector<Mesh::shared_ptr> getMeshes() const;


#ifdef USE_AMP_VECTORS
    /**
     * \brief    Displace the entire mesh
     * \details  This function will displace the entire mesh by displacing
     *   each node by the values provided in the vector.  This function is 
     *   a blocking call for the mesh communicator
     * \param x  Displacement vector.  Must have N DOFs per node where N 
     *           is the physical dimension of the mesh.
     */
    virtual void displaceMesh ( boost::shared_ptr<const AMP::LinearAlgebra::Vector> x );
#endif

    // Needed to prevent problems with virtual functions
    using Mesh::Subset;

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

