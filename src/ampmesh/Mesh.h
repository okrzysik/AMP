#ifndef included_AMP_Mesh
#define included_AMP_Mesh

#include "AMP/ampmesh/MeshID.h"
#include "AMP/ampmesh/MeshIterator.h"
#include "AMP/ampmesh/MeshParameters.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/enable_shared_from_this.h"
#include <memory>


// Forward declerations
namespace AMP {
namespace Geometry {
class Geometry;
}
#ifdef USE_AMP_VECTORS
namespace LinearAlgebra {
class Vector;
}
#endif
} // namespace AMP


namespace AMP {
namespace Mesh {


//! Enumeration for basic mesh-based quantities
enum class SetOP { Union, Intersection, Complement };


/**
 * \class Mesh
 * \brief A class used to abstract away mesh from an application
 *
 * \details  This class provides routines for reading, accessing and writing meshes.
 * The database fields control the mesh and will differ for each mesh type.  However,
 * there are some common fields to all meshes:
 *     MeshName - The name to associate with the mesh
 *     MeshType - The mesh type (libMesh, Multimesh, AMP)
 *     x_offset - Optional argument specifying the offset in the x-direction
 *     y_offset - Optional argument specifying the offset in the y-direction
 *     z_offset - Optional argument specifying the offset in the z-direction
 *     NumberOfElements - Optional argument indicating the number of elements in the mesh (will
 * override all other
 * calulations)
 *     Weight - Optional argument indicating the relative weight of the mesh for the domain
 * decomposition (relative to
 * 1.0)
 */
class Mesh : public AMP::enable_shared_from_this<AMP::Mesh::Mesh>
{
public:
    /**
     *\typedef shared_ptr
     *\brief  Name for the shared pointer.
     *\details  Use this typedef for a reference counted pointer to a mesh manager object.
     */
    typedef std::shared_ptr<AMP::Mesh::Mesh> shared_ptr;

    /**
     *\typedef const_shared_ptr
     *\brief  Name for the const shared pointer.
     *\details  Use this typedef for a reference counted pointer to a mesh manager object.
     */
    typedef std::shared_ptr<const AMP::Mesh::Mesh> const_shared_ptr;

    /**
     *\typedef generator
     *\brief  Generator for meshes
     *\details  This is a user-supplied function to generate a mesh.  Users may register their
     *     own mesh generators using registerGenerator and the mesh builder will call them.
     */
    typedef std::function<Mesh::shared_ptr( MeshParameters::shared_ptr )> generatorType;

    //! Enumeration for basic mesh-based quantities
    enum class Movable : uint8_t { Fixed = 0, Displace = 1, Deform = 2 };


    /**
     * \brief Read in mesh files, partition domain, and prepare environment for simulation
     * \details  For trivial parallelsim, this method reads in the meshes on each processor.  Each
     * processor contains a piece of each mesh.  For massive parallelism, each mesh is on its own
     * communicator.  As such, some math libraries must be initialized accordingly.
     * \param params  Parameters for constructing a mesh from an input database
     */
    explicit Mesh( const MeshParameters::shared_ptr &params );


    /**
     * \brief Construct a new mesh from an existing mesh.
     * \details  This constructor will construct a new mesh from an existing mesh
     * using an iterator over the existing mesh.
     * This is designed as a path to create a new mesh object of one type from
     * an existing mesh of a different type.  It also allows creating a new single mesh
     * from a subset or superset of other meshes.  Note that instantion of this routine
     * may not be able to create it's mesh from any arbitrary mesh, and may throw an
     * error.
     * \param old_mesh  Existing mesh that we will use to construct the new mesh
     * \param iterator  Iterator over the existing mesh
     */
    Mesh( const Mesh::shared_ptr &old_mesh, MeshIterator::shared_ptr &iterator );


    /**
     * \brief   Create a mesh
     * \details  This function will create a mesh (or series of meshes) based on
     *   the input database.
     * \param params Parameters for constructing a mesh from an input database
     */
    static std::shared_ptr<AMP::Mesh::Mesh> buildMesh( const MeshParameters::shared_ptr &params );


    /**
     * \brief   Create a mesh
     * \details  This function will create a mesh (or series of meshes) based on
     *   the input database.
     * \param name      Name of mesh generator
     * \param gen       Mesh generator to use
     */
    static inline void registerGenerator( const std::string &name, const generatorType &gen )
    {
        d_generators[name] = gen;
    }


    /**
     * \brief   Return the geometry of the mesh
     * \details  This function will return the geometry for the mesh if it exists.
     *    Not all meshes will have a geometry associated with them.
     */
    inline auto getGeometry() { return d_geometry; }


    /**
     * \brief   Return the geometry of the mesh
     * \details  This function will return the geometry for the mesh if it exists.
     *    Not all meshes will have a geometry associated with them.
     */
    inline auto getGeometry() const { return d_geometry; }


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


    //! Deconstructor
    virtual ~Mesh();


    //! Virtual function to clone the mesh (allows us to properly copy the derived class)
    virtual std::unique_ptr<Mesh> clone() const = 0;


    /**
     * \brief    Subset a mesh given a MeshID
     * \details  This function will return the mesh with the given meshID.
     *    Note: for multimeshes, this will return the mesh with the given id.
     *    For a single mesh this will return a pointer to itself if the meshID
     *    matches the meshID of the mesh, and a null pointer otherwise.
     * \param meshID  MeshID of the desired mesh
     */
    virtual std::shared_ptr<Mesh> Subset( MeshID meshID ) const;


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
    virtual std::shared_ptr<Mesh> Subset( std::string name ) const;


    /**
     * \brief    Subset a mesh given a MeshIterator
     * \details  This function will subset a mesh over a given iterator.
     *   This will return a new mesh object.
     * \param iterator  MeshIterator used to subset
     * \param isGlobal  Is the new subset mesh global over the entire mesh (true,default),
     *                  or do we only want to keep the local mesh (false)
     */
    virtual std::shared_ptr<Mesh> Subset( const MeshIterator &iterator,
                                          bool isGlobal = true ) const;


    /**
     * \brief        Subset a mesh given another mesh
     * \details      This function will subset a mesh given another mesh
     * \param mesh   Mesh used to subset
     */
    virtual std::shared_ptr<Mesh> Subset( Mesh &mesh ) const;


    /* Return the number of local element of the given type
     * \param type   Geometric type
     */
    virtual size_t numLocalElements( const GeomType type ) const;


    /* Return the global number of elements of the given type
     * Note: depending on the mesh this routine may require global communication across the mesh.
     * \param type   Geometric type
     */
    virtual size_t numGlobalElements( const GeomType type ) const;


    /* Return the number of ghost elements of the given type on the current processor
     * \param type   Geometric type
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
     * \param id     Block id for the elements (example: block id in cubit, subdomain in libmesh)
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator
    getBlockIDIterator( const GeomType type, const int id, const int gcw = 0 ) const;


    /**
     * \brief    Return an MeshIterator constructed through a set operation of two other
     * MeshIterators.
     * \details  Return an MeshIterator constructed through a set operation of two other
     * MeshIterators.
     * \param OP Set operation to perform.
     *           SetOP::Union - Perform a union of the iterators ( A U B )
     *           SetOP::Intersection - Perform an intersection of the iterators ( A n B )
     *           SetOP::Complement - Perform a compliment of the iterators ( A - B )
     * \param A  Pointer to MeshIterator A
     * \param B  Pointer to MeshIterator B
     */
    static MeshIterator getIterator( SetOP OP, const MeshIterator &A, const MeshIterator &B );


    /**
     * \brief    Check if an element is in the mesh
     * \details  This function queries the mesh to determine if the given element is a member of the
     * mesh
     * \param id    Mesh element id we are querying.
     */
    virtual bool isMember( const MeshElementID &id ) const;


    /**
     * \brief    Return a mesh element given it's id.
     * \details  This function queries the mesh to get an element given the mesh id.
     *    This function is only required to return an element if the id is local.
     *    Ideally, this should be done in O(1) time, but the implimentation is up to
     *    the underlying mesh.  The base class provides a basic implimentation, but
     *    uses mesh iterators and requires O(N) time on the number of elements in the mesh.
     * \param id    Mesh element id we are requesting.
     */
    virtual MeshElement getElement( const MeshElementID &id ) const;


    /**
     * \brief    Return the parent elements of the given mesh element
     * \details  This function queries the mesh to get an element given the mesh id,
     *    then returns the parent elements that have the element as a child
     * \param elem  Mesh element of interest
     * \param type  Element type of the parents requested
     */
    virtual std::vector<MeshElement> getElementParents( const MeshElement &elem,
                                                        const GeomType type ) const;


    //! Get the largest geometric type in the mesh
    inline GeomType getGeomType() const { return GeomDim; }


    //! Get the physical dimension of the mesh
    inline uint8_t getDim() const { return PhysicalDim; }


    //! Get the communicator for the mesh
    inline const AMP_MPI &getComm() const { return d_comm; }


    //! Get the maximum ghost width
    inline uint8_t getMaxGhostWidth() const { return d_max_gcw; }


    //! Get the mesh ID
    inline MeshID meshID() const { return d_meshID; }


    //! Is the current mesh a base mesh
    virtual inline bool isBaseMesh() const { return true; }


    /**
     *  Get the meshIDs of all meshes that compose the current mesh (including its self)
     *  Note: This function may require global communication depending on the implimentation
     */
    virtual std::vector<MeshID> getAllMeshIDs() const;


    /**
     *  Get the meshIDs of all the basic meshes that compose the current mesh
     *     (excluding multimeshes and subset meshes)
     *  Note: This function may require global communication depending on the implimentation
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


    //! Get the mesh name
    virtual inline std::string getName() const { return d_name; }


    //! Set the mesh name
    virtual inline void setName( std::string name ) { d_name = name; }


    /**
     * \brief    Get the bounding box for the mesh
     * \details  This function will return the bounding box for the entire mesh.
     *   The vector returned contains the box that contains the mesh in the form
     *   [ x_min  x_max  y_min  y_max  z_min  z_max ].
     */
    virtual std::vector<double> getBoundingBox() const { return d_box; }


    /**
     * \brief    Get the bounding box for the local part of the mesh
     * \details  This function will return the bounding box for the local part of the mesh.
     *   The vector returned contains the box that contains the mesh in the form
     *   [ x_min  x_max  y_min  y_max  z_min  z_max ].
     */
    virtual std::vector<double> getLocalBoundingBox() const { return d_box_local; }


    /**
     * \brief    Is the mesh movable
     * \details  This function will check if the mesh can be displaced.
     * @return   enum indicating the extent the mesh can be moved
     */
    virtual Movable isMeshMovable() const = 0;


    /**
     * \brief    Identify if the position has moved
     * \details  This function will return a hash that can be used to
     *    identify if the mesh has been moved.  Any time that displaceMesh
     *    is called, the hash value should change.  There is no requirement
     *    that dispacing a mesh and returning it back to the original position
     *    will return the original hash.
     * @return   hash value with current position id
     */
    virtual uint64_t positionHash() const = 0;


    /**
     * \brief    Displace the entire mesh
     * \details  This function will displace the entire mesh by a scalar value.
     *   This function is a blocking call for the mesh communicator, and requires
     *   the same value on all processors.  The displacement vector should be the
     *   size of the physical dimension.
     * \param x  Displacement vector
     */
    virtual void displaceMesh( const std::vector<double> &x ) = 0;


#ifdef USE_AMP_VECTORS
    /**
     * \brief    Displace the entire mesh
     * \details  This function will displace the entire mesh by displacing
     *   each node by the values provided in the vector.  This function is
     *   a blocking call for the mesh communicator
     * \param x  Displacement vector.  Must have N DOFs per node where N
     *           is the physical dimension of the mesh.
     */
    virtual void displaceMesh( std::shared_ptr<const AMP::LinearAlgebra::Vector> x ) = 0;


    /**
     * \brief    Get a vector of the coordinates of the nodes
     * \details  This function will return a const vector containing the coordinates of
     *           all the nodes.
     * \param name   Name of the vector
     * \param gcw    Desired ghost cell width
     */
    virtual std::shared_ptr<AMP::LinearAlgebra::Vector>
    getPositionVector( std::string name, const int gcw = 0 ) const;
#endif

    std::shared_ptr<AMP::Database> DB() const { return d_db; }

protected:
    //!  Empty constructor for a mesh
    Mesh() {}

    //! The mesh parameters
    MeshParameters::shared_ptr d_params;

    //! The geometry parameters
    std::shared_ptr<Geometry::Geometry> d_geometry;

    //! The geometric dimension (equivalent to the highest geometric object that could be
    //! represented)
    GeomType GeomDim;

    //! The physical dimension
    uint8_t PhysicalDim;

    //! The physical dimension
    uint8_t d_max_gcw;

    //! The communicator over which the mesh is stored
    AMP_MPI d_comm;

    //! A pointer to an AMP database containing the mesh info
    std::shared_ptr<AMP::Database> d_db;

    //! A unique id for each mesh
    MeshID d_meshID;

    //! A name for the mesh
    std::string d_name;

    //! The bounding box for the mesh
    std::vector<double> d_box, d_box_local;

    //! A list of mesh generators to use
    static std::map<std::string, generatorType> d_generators;

    /**
     *  A function to create a unique id for the mesh (requires the comm to be set)
     *  Note: this requires a global communication across the mesh communicator.
     *  Note: this function is NOT thread safe, and will need to be modified before threads are
     * used.
     */
    void setMeshID();

    // Private copy constructor
    explicit Mesh( const Mesh &old_mesh );

    // Private assigment operator
    Mesh &operator=( const Mesh &old_mesh ) = delete;
};

} // namespace Mesh
} // namespace AMP

#endif
