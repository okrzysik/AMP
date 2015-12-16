#ifndef included_AMP_STKMesh
#define included_AMP_STKMesh

#include "ampmesh/Mesh.h"
#include "ampmesh/MeshIterator.h"
#include "ampmesh/STKmesh/initializeSTKMesh.h"

#include "Teuchos_RCP.hpp"
#include "stk_io/util/Gmesh_STKmesh_Fixture.hpp"
#include "stk_io/util/IO_Fixture.hpp"
#include "stk_mesh/fem/FEMMetaData.hpp"


namespace AMP {
namespace Mesh {


class STKMeshElement;


/**
 * \class STKMesh
 * \brief A concrete mesh class for STKMesh
 *
 * \details  This class provides routines for reading, accessing and writing STKMesh meshes.
 * The generation of the mesh is controlled by the database passed in through the params object.
 * The database fields control the mesh and provide several options:
 * @code
 *    dim - required integer specifying the physical dimension
 *    FileName - If specified this will load the mesh from the given file
 *    Generator - If specified this will generate a new mesh using the optional parameters in the
 * database
 *       This field must be a string specifying the generator to use.  Valid gerators are:
 *          "cube" - Will generate a cube mesh
 *       Additional areguments:
 *          size - Integer array specifying the number of elements in each direction
 * @endcode
 * The parallel decomposition of the mesh is controlled by STKmesh and occurs on the communicator
 * specified through the params object.  Note that STKmesh does not support meshes on overlapping
 * communicators.  If multiple meshes are used, they must either share communicators or have unique
 * communicators.
 */
class STKMesh : public Mesh
{
public:
    /**
     * \brief Read in mesh files, partition domain, and prepare environment for simulation
     * \details  For trivial parallelsim, this method reads in the meshes on each processor.  Each
     * processor contains a piece of each mesh.  For massive parallelism, each mesh is on its own
     * communicator.  As such, some math STKraries must be initialized accordingly.
     * \param params Parameters for constructing a mesh from an input database
     */
    STKMesh( const MeshParameters::shared_ptr &params );

    /**
     * \brief Contructor to create a STKMesh object from a STKMesh mesh.
     * \details This constructor allows the user to construct a mesh directly in STKmesh
     * and use it to create the new mesh object.  This function is intended for advanced
     * users only.  Note: the user is responsible for ensuring that STKMesh is properly
     * initialized (using initializeSTKMesh), and that the mesh is created properly.
     * \param mesh The mesh in STKmesh we want to use to construct the new mesh object
     * \param name The name of the new mesh object
     */
    STKMesh( AMP::shared_ptr<stk::mesh::BulkData> mesh, std::string name );

    //! Deconstructor
    virtual ~STKMesh();

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


    /**
     * \brief   Return the maximum number of processors that can be used with the mesh
     * \details  This function will return the maximum number of processors that can
     *   be used with the mesh.
     * \param params Parameters for constructing a mesh from an input database
     */
    static size_t maxProcs( const MeshParameters::shared_ptr &params );


    /* Return the number of local element of the given type
     * \param type   Geometric type
     */
    virtual size_t numLocalElements( const GeomType type ) const;


    /* Return the global number of elements of the given type
     * \param type   Geometric type
     */
    virtual size_t numGlobalElements( const GeomType type ) const;


    /* Return the number of ghost elements of the given type on the current processor
     * \param type   Geometric type
     * \param gcw    Desired ghost cell width
     */
    virtual size_t numGhostElements( const GeomType type, const int gcw ) const;


    /**
     * \brief    Return an MeshIterator over the given geometric objects
     * \details  Return an MeshIterator over the given geometric objects
     * \param type   Geometric type to iterate over
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator getIterator( const GeomType type, const int gcw = 0 ) const;


    /**
     * \brief    Return an MeshIterator over the given geometric objects on the surface
     * \details  Return an MeshIterator over the given geometric objects on the surface
     * \param type   Geometric type to iterate over
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator getSurfaceIterator( const GeomType type, const int gcw = 0 ) const;


    /**
     * \brief    Return the list of all boundary ID sets in the mesh
     * \details  Return the list of all boundary ID sets in the mesh
     * Note: depending on the mesh this routine may require global communication across the mesh.
     */
    virtual std::vector<int> getBoundaryIDs() const;


    /**
     * \brief    Return an MeshIterator over the given geometric objects on the given boundary ID
     * set
     * \details  Return an MeshIterator over the given geometric objects on the given boundary ID
     * set
     * \param type   Geometric type to iterate over
     * \param id     Boundary id for the elements (example: sideset id)
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator
    getBoundaryIDIterator( const GeomType type, const int id, const int gcw = 0 ) const;

    /**
     * \brief    Return the list of all boundary ID sets in the mesh
     * \details  Return the list of all boundary ID sets in the mesh
     * Note: depending on the mesh this routine may require global communication across the mesh.
     */
    virtual std::vector<int> getBlockIDs() const;


    /**
     * \brief    Return an MeshIterator over the given geometric objects on the given block ID set
     * \details  Return an MeshIterator over the given geometric objects on the given block ID set
     * \param type   Geometric type to iterate over
     * \param id     Block id for the elements (example: block id in cubit, subdomain in STKmesh)
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator
    getBlockIDIterator( const GeomType type, const int id, const int gcw = 0 ) const;


    /**
     * \brief    Return a mesh element given it's id.
     * \details  This function queries the mesh to get an element given the mesh id.
     *    This function is only required to return an element if the id is local.
     *    This function uses STKmesh for elements it owns, and uses an O(ln(N)) for
     *    elements that were constructed internally.
     * \param id    Mesh element id we are requesting.
     */
    virtual MeshElement getElement( const MeshElementID &id ) const;


    /**
     * \brief    Displace the entire mesh
     * \details  This function will displace the entire mesh by a scalar value.
     *   This function is a blocking call for the mesh communicator, and requires
     *   the same value on all processors.  The displacement vector should be the
     *   size of the physical dimension.
     * \param x  Displacement vector
     */
    virtual void displaceMesh( std::vector<double> x );


    //! Return the underlying STKMesh object
    inline AMP::shared_ptr<stk::mesh::BulkData> getSTKMeshBulk() const { return d_STKMeshBulk; }
    inline AMP::shared_ptr<stk::mesh::fem::FEMMetaData> getSTKMeshMeta() const
    {
        return d_STKMeshMeta;
    }


#ifdef USE_AMP_VECTORS
    /**
     * \brief    Displace the entire mesh
     * \details  This function will displace the entire mesh by displacing
     *   each node by the values provided in the vector.  This function is
     *   a blocking call for the mesh communicator
     * \param x  Displacement vector.  Must have N DOFs per node where N
     *           is the physical dimension of the mesh.
     */
    virtual void displaceMesh( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> x );
#endif


protected:
    /**
     *  Function to return the neighbors for a node.  The neighbors are defined
     *  as those nodes that share an element with the given node.
     *  Note: the nodes returns are returned in unsorted order.
     */
    std::vector<stk::mesh::Entity *> getNeighborNodes( MeshElementID ) const;

    // Friend functions to access protected functions
    friend class STKMeshElement;

private:
    //!  Empty constructor for a mesh
    STKMesh();

    //!  Function to properly initialize the internal data once a STKmesh mesh is loaded
    void initialize();

    // STKMesh objects
    AMP::shared_ptr<stk::mesh::BulkData> d_STKMeshBulk;
    AMP::shared_ptr<stk::mesh::fem::FEMMetaData> d_STKMeshMeta;
    Teuchos::RCP<Ioss::Region> d_STKIORegion;
    AMP::shared_ptr<stk::io::util::IO_Fixture> d_STKIOFixture;
    AMP::shared_ptr<stk::io::util::Gmesh_STKmesh_Fixture> d_STKGMeshFixture;

    // Some basic internal data
    std::vector<unsigned> n_local;
    std::vector<unsigned> n_global;
    std::vector<unsigned> n_ghost;
    AMP::shared_ptr<initializeSTKMesh> STKmeshInit;

    // Data used to store the node neighbor lists
    std::vector<unsigned int> neighborNodeIDs;
    std::vector<std::vector<stk::mesh::Entity *>> neighborNodes;

    // Data used to elements that STKmesh doesn't create
    std::map<GeomType, AMP::shared_ptr<std::vector<MeshElement>>> d_localElements;
    std::map<GeomType, AMP::shared_ptr<std::vector<MeshElement>>> d_ghostElements;

    // Data used to store the boundary elements
    std::map<std::pair<int, GeomType>, AMP::shared_ptr<std::vector<MeshElement>>> d_boundarySets;

    // Data used to store the surface elements
    std::vector<AMP::shared_ptr<std::vector<MeshElement>>> d_localSurfaceElements;
    std::vector<AMP::shared_ptr<std::vector<MeshElement>>> d_ghostSurfaceElements;

    // Data used to store block info
    std::vector<int> d_block_ids;
};

} // Mesh namespace
} // AMP namespace

#endif
