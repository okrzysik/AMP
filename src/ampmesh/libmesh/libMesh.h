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
 */
class libMesh: public Mesh
{
public:

    /**
     * \param params Parameters for constructing a mesh from an input database
     * \brief Read in mesh files, partition domain, and prepare environment for simulation
     * \details  For trivial parallelsim, this method reads in the meshes on each processor.  Each
     * processor contains a piece of each mesh.  For massive parallelism, each mesh is on its own
     * communicator.  As such, some math libraries must be initialized accordingly.
     */
    libMesh ( const MeshParameters::shared_ptr &params );

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
     */
    virtual size_t  numGhostElements( const GeomType type, const int gcw ) const;


    /**
     * \brief    Return an MeshIterator over the given geometric objects
     * \details  Return an MeshIterator over the given geometric objects
     * \param type   Geometric type to iterate over
     * \param gcw    Desired ghost 
     */
    virtual MeshIterator getIterator ( const GeomType type, const int gcw=0 );


    /**
     * \brief    Return the list of all ID sets in the mesh
     * \details  Return the list of all ID sets in the mesh
     */
    virtual std::vector<int> getIDSets ( );


    /**
     * \brief    Return an MeshIterator over the given geometric objects on the given ID set
     * \details  Return an MeshIterator over the given geometric objects on the given ID set
     * \param type   Geometric type to iterate over
     * \param type   id for the elements (example: nodeset id)
     */
    virtual MeshIterator getIDsetIterator ( const GeomType type, const int id, const int gcw=0 );



    //! Return the underlying libMesh object
    inline boost::shared_ptr< ::Mesh> getlibMesh( ) const { return d_libMesh; }

protected:

    /**  
     *  Function to return the neighbors for a node.  The neighbors are defined
     *  as those nodes that share an element with the given node.
     *  Note: the nodes returns are returned in unsorted order.
     */
    std::vector< ::Node* > getNeighborNodes( MeshElementID );

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

    // Data use to store the node neighbor lists
    std::vector<unsigned int> neighborNodeIDs;
    std::vector< std::vector< ::Node* > > neighborNodes;

};

} // Mesh namespace
} // AMP namespace

#endif

