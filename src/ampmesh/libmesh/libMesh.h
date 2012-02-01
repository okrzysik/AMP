#ifndef included_AMP_LibMesh
#define included_AMP_LibMesh

#include "ampmesh/Mesh.h"
#include "ampmesh/libmesh/initializeLibMesh.h"

// LibMesh include
#include "mesh.h"
#include "mesh_data.h"


namespace AMP {
namespace Mesh {


class libMeshElement;


/**
 * \class libMesh
 * \brief A concrete mesh class for libMesh
 *
 * \details  This class provides routines for reading, accessing and writing libMesh meshes.
 * The generation of the mesh is controlled by the database passed in through the params object.
 * The database fields control the mesh and provide several options: 
 * @code
 *    dim - required integer specifying the physical dimension
 *    FileName - If specified this will load the mesh from the given file
 *    Generator - If specified this will generate a new mesh using the optional parameters in the database
 *       This field must be a string specifying the generator to use.  Valid gerators are:
 *          "cube" - Will generate a cube mesh
 *       Additional areguments:
 *          size - Integer array specifying the number of elements in each direction
 * @endcode
 * The parallel decomposition of the mesh is controlled by libmesh and occurs on the communicator
 * specified through the params object.  Note that libmesh does not support meshes on overlapping 
 * communicators.  If multiple meshes are used, they must either share communicators or have unique 
 * communicators.
 */
class libMesh: public Mesh
{
public:

    /**
     * \brief Read in mesh files, partition domain, and prepare environment for simulation
     * \details  For trivial parallelsim, this method reads in the meshes on each processor.  Each
     * processor contains a piece of each mesh.  For massive parallelism, each mesh is on its own
     * communicator.  As such, some math libraries must be initialized accordingly.
     * \param params Parameters for constructing a mesh from an input database
     */
    libMesh ( const MeshParameters::shared_ptr &params );

    /**
     * \brief Contructor to create a libMesh object from a libMesh mesh.
     * \details This constructor allows the user to construct a mesh directly in libmesh
     * and use it to create the new mesh object.  This function is intended for advanced
     * users only.  Note: the user is responsible for ensuring that libMesh is properly 
     * initialized (using initializeLibMesh), and that the mesh is created properly.
     * \param mesh The mesh in libmesh we want to use to construct the new mesh object
     * \param name The name of the new mesh object
     */
    libMesh ( boost::shared_ptr< ::Mesh> mesh, std::string name );

    //! Deconstructor
     ~libMesh ();

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


    /* Return the global number of elements of the given type
     * \param type   Geometric type
     */
    virtual size_t  numGlobalElements( const GeomType type ) const;


    /* Return the number of ghost elements of the given type on the current processor
     * \param type   Geometric type
     * \param gcw    Desired ghost cell width
     */
    virtual size_t  numGhostElements( const GeomType type, const int gcw ) const;


    /**
     * \brief    Return an MeshIterator over the given geometric objects
     * \details  Return an MeshIterator over the given geometric objects
     * \param type   Geometric type to iterate over
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator  getIterator ( const GeomType type, const int gcw=0 ) const;


    /**
     * \brief    Return an MeshIterator over the given geometric objects on the surface
     * \details  Return an MeshIterator over the given geometric objects on the surface
     * \param type   Geometric type to iterate over
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator  getSurfaceIterator ( const GeomType type, const int gcw=0 ) const;


    /**
     * \brief    Return the list of all ID sets in the mesh
     * \details  Return the list of all ID sets in the mesh
     */
    virtual std::vector<int>  getIDSets ( ) const;


    /**
     * \brief    Return an MeshIterator over the given geometric objects on the given ID set
     * \details  Return an MeshIterator over the given geometric objects on the given ID set
     * \param type   Geometric type to iterate over
     * \param id     id for the elements (example: nodeset id)
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator  getIDsetIterator ( const GeomType type, const int id, const int gcw=0 ) const;


    /**
     * \brief    Displace the entire mesh
     * \details  This function will displace the entire mesh by a scalar value.
     *   This function is a blocking call for the mesh communicator, and requires
     *   the same value on all processors.  The displacement vector should be the 
     *   size of the physical dimension.
     * \param x  Displacement vector
     */
    virtual void displaceMesh( std::vector<double> x );


    //! Return the underlying libMesh object
    inline boost::shared_ptr< ::Mesh> getlibMesh( ) const { return d_libMesh; }


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


protected:

    /**  
     *  Function to return the neighbors for a node.  The neighbors are defined
     *  as those nodes that share an element with the given node.
     *  Note: the nodes returns are returned in unsorted order.
     */
    std::vector< ::Node* > getNeighborNodes( MeshElementID ) const;

    // Friend functions to access protected functions    
    friend class libMeshElement;

private:

    //!  Empty constructor for a mesh
    libMesh( ) {};

    //!  Function to properly initialize the internal data once a libmesh mesh is loaded
    void initialize( );

    // libMesh objects
    boost::shared_ptr< ::Mesh>          d_libMesh;
    boost::shared_ptr< ::MeshData>      d_libMeshData;

    // Some basic internal data
    std::vector<size_t> n_local, n_global, n_ghost;
    boost::shared_ptr<initializeLibMesh> libmeshInit;

    // Data used to store the node neighbor lists
    std::vector<unsigned int> neighborNodeIDs;
    std::vector< std::vector< ::Node* > > neighborNodes;

    // Data used to elements that libmesh doesn't create
    std::map< GeomType, boost::shared_ptr<std::vector<MeshElement> > >  d_localElements;
    std::map< GeomType, boost::shared_ptr<std::vector<MeshElement> > >  d_ghostElements;

    // Data used to store the boundary elements
    std::map< std::pair<int,GeomType>, boost::shared_ptr<std::vector<MeshElement> > >  d_boundarySets;

    // Data used to store the surface elements
    std::vector< boost::shared_ptr<std::vector<MeshElement> > >  d_localSurfaceElements;
    std::vector< boost::shared_ptr<std::vector<MeshElement> > >  d_ghostSurfaceElements;

};

} // Mesh namespace
} // AMP namespace

#endif

