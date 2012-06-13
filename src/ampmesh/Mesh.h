#ifndef included_AMP_Mesh
#define included_AMP_Mesh

#include "ampmesh/MeshParameters.h"
#include "ampmesh/MeshID.h"
#include "ampmesh/MeshIterator.h"
#include "utils/AMP_MPI.h"

#ifdef USE_AMP_VECTORS
    namespace AMP {
    namespace LinearAlgebra {
        class Vector;
    }
    }
#endif

#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

namespace AMP {
namespace Mesh {


//! Enumeration for basic mesh-based quantities
enum SetOP { Union, Intersection, Complement };


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
 *     NumberOfElements - Optional argument indicating the number of elements in the mesh (will override all other calulations)
 *     Weight - Optional argument indicating the relative weight of the mesh for the domain decomposition (relative to 1.0)
 */
class Mesh: public boost::enable_shared_from_this<AMP::Mesh::Mesh>
{
public:

    /**
     *\typedef shared_ptr
     *\brief  Name for the shared pointer.
     *\details  Use this typedef for a reference counted pointer to a mesh manager object.
     */
    typedef boost::shared_ptr<AMP::Mesh::Mesh>  shared_ptr;


    /**
     * \brief Read in mesh files, partition domain, and prepare environment for simulation
     * \details  For trivial parallelsim, this method reads in the meshes on each processor.  Each
     * processor contains a piece of each mesh.  For massive parallelism, each mesh is on its own
     * communicator.  As such, some math libraries must be initialized accordingly.
     * \param params  Parameters for constructing a mesh from an input database
     */
    Mesh ( const MeshParameters::shared_ptr &params );


    /**
     * \brief Construct a new mesh from an existing mesh
     * \details  This constructor will construct a new mesh from an existing mesh.
     * This is designed as a path to create a new mesh object of one type from
     * an existing mesh of a different type.  It also allows creating a new single mesh
     * from a subset or superset of other meshes.  Note that instantion of this routine 
     * may not be able to create it's mesh from any arbitrary mesh, and may throw an 
     * error.
     * \param old_mesh  Existing mesh that we will use to construct the new mesh
     */
    Mesh ( const Mesh::shared_ptr &old_mesh );


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
    Mesh ( const Mesh::shared_ptr &old_mesh, MeshIterator::shared_ptr &iterator);


    /**
     * \brief   Create a mesh 
     * \details  This function will create a mesh (or series of meshes) based on
     *   the input database.  
     * \param params Parameters for constructing a mesh from an input database
     */
    static boost::shared_ptr<AMP::Mesh::Mesh> buildMesh( const MeshParameters::shared_ptr &params );


    /**
     * \brief   Estimate the number of elements in the mesh 
     * \details  This function will estimate the number of elements in the mesh. 
     *   This is used so that we can properly balance the meshes across multiple processors.
     *   Ideally this should be both an accurate estimate and very fast.  It should not require
     *   any communication and should not have to actually load a mesh.
     * \param params Parameters for constructing a mesh from an input database
     */
    static size_t estimateMeshSize( const MeshParameters::shared_ptr &params );


    //! Deconstructor
     ~Mesh ();


    //! Assignment operator
    virtual Mesh operator=(const Mesh&);


    //! Virtual function to copy the mesh (allows use to proply copy the derived class)
    virtual Mesh copy() const;


    /**
     * \brief    Subset a mesh given a MeshID
     * \details  This function will return the mesh with the given meshID.
     *    Note: for multimeshes, this will return the mesh with the given id.
     *    For a single mesh this will return a pointer to itself if the meshID
     *    matches the meshID of the mesh, and a null pointer otherwise.
     * \param meshID  MeshID of the desired mesh
     */
    virtual boost::shared_ptr<Mesh>  Subset ( MeshID meshID ) const;


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
     * \brief    Subset a mesh given a MeshIterator
     * \details  This function will subset a mesh over a given iterator.
     *   This will return a new mesh object.
     * \param iterator  MeshIterator used to subset
     */
    virtual boost::shared_ptr<Mesh>  Subset ( const MeshIterator &iterator ) const;


    /**
     * \brief        Subset a mesh given another mesh
     * \details      This function will subset a mesh given another mesh
     * \param mesh   Mesh used to subset
     */
    virtual boost::shared_ptr<Mesh>  Subset ( Mesh &mesh ) const;


    /* Return the number of local element of the given type
     * \param type   Geometric type
     */
    virtual size_t  numLocalElements( const GeomType type ) const;


    /* Return the global number of elements of the given type
     * Note: depending on the mesh this routine may require global communication across the mesh.
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
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator getIterator ( const GeomType type, const int gcw=0 ) const;


    /**
     * \brief    Return an MeshIterator over the given geometric objects on the surface
     * \details  Return an MeshIterator over the given geometric objects on the surface
     * \param type   Geometric type to iterate over
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator getSurfaceIterator ( const GeomType type, const int gcw=0 ) const;


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
     * \brief    Return an MeshIterator constructed through a set operation of two other MeshIterators.
     * \details  Return an MeshIterator constructed through a set operation of two other MeshIterators.
     * \param OP Set operation to perform.
     *           Union - Perform a union of the iterators ( A U B )
     *           Intersection - Perform an intersection of the iterators ( A n B )
     *           Complement - Perform a compliment of the iterators ( A - B )
     * \param A  Pointer to MeshIterator A
     * \param B  Pointer to MeshIterator B
     */
    static MeshIterator getIterator ( SetOP OP, const MeshIterator &A, const MeshIterator &B);
 

    /**
     * \brief    Return a mesh element given it's id.
     * \details  This function queries the mesh to get an element given the mesh id.
     *    This function is only required to return an element if the id is local.
     *    Ideally, this should be done in O(1) time, but the implimentation is up to 
     *    the underlying mesh.  The base class provides a basic implimentation, but 
     *    uses mesh iterators and requires O(N) time on the number of elements in the mesh.
     * \param id    Mesh element id we are requesting.
     */
    virtual MeshElement getElement ( const MeshElementID &id ) const;
 

    //! Get the largest geometric type in the mesh
    virtual GeomType getGeomType() const { return GeomDim; } 


    //! Get the physical dimension of the mesh
    virtual unsigned char getDim() const { return PhysicalDim; } 


    //! Get the largest geometric type in the mesh
    virtual AMP_MPI getComm() const { return d_comm; }


    //! Get the maximum ghost width
    virtual unsigned char getMaxGhostWidth() const { return d_max_gcw; }


    //! Get the mesh ID
    virtual inline MeshID meshID() const { return d_meshID; }


    //! Is the current mesh a base mesh
    virtual inline bool isBaseMesh() const { return true; }


    /**
     *  Get the meshIDs of all meshes that compose the current mesh (including its self)
     *  Note: This function may require global communication depending on the implimentation
     */
    virtual std::vector<MeshID> getAllMeshIDs() const;


    /**
     *  Get the meshIDs of all the basic meshes that compose the current mesh (excluding multimeshes and subset meshes)
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
    virtual inline void setName(std::string name) { d_name = name; }


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
     * \brief    Displace the entire mesh
     * \details  This function will displace the entire mesh by a scalar value.
     *   This function is a blocking call for the mesh communicator, and requires
     *   the same value on all processors.  The displacement vector should be the 
     *   size of the physical dimension.
     * \param x  Displacement vector
     */
    virtual void displaceMesh( std::vector<double> x );


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


    /**
     * \brief    Get a vector of the coordinates of the nodes
     * \details  This function will return a const vector containing the coordinates of 
     *           all the nodes.  
     * \param name   Name of the vector
     * \param gcw    Desired ghost cell width
     */
    virtual boost::shared_ptr<AMP::LinearAlgebra::Vector>  getPositionVector( std::string name, const int gcw=0 ) const;
#endif


    /**
     * \struct MultiMesh
     * \brief Structure used to contain simulated mesh load
     * \details  This structure provides info that can be used to simulate loading 
     *   a mesh, and checking the resulting load balance
     */
    struct simulated_mesh_struct {
        std::string name;                               //!< The mesh name
        std::string type;                               //!< The mesh type
        size_t N_elements;                              //!< The number of elements in the mesh
        boost::shared_ptr<const Database> db;           //!< The database used for the mesh
        std::vector<int> ranks;                         //!< The ranks of the processors that own a piece of the mesh
        std::vector<Mesh::simulated_mesh_struct> submeshes;   //!< Sub-meshes to the current mesh
        simulated_mesh_struct() {};                     //!< Empty constructor
        simulated_mesh_struct(const simulated_mesh_struct&);  //!< Copy constructor
        void print();                                   //!< Function to print the mesh hierarchy
    };

    
    /**
     * \brief    Simulate the mesh build process
     * \details  This function will simulate the loading and load balancing of the mesh hierarchy
     * \param params        Parameters to use for the mesh construction
     * \param comm_ranks    Simulated ranks that are used to create the mesh
     */
    static simulated_mesh_struct  simulateBuildMesh( const MeshParameters::shared_ptr &params, std::vector<int> &comm_ranks );


protected:

    //!  Empty constructor for a mesh
    Mesh() {}

    //! The mesh parameters
    MeshParameters::shared_ptr params;

    //! The geometric dimension (equivalent to the highest geometric object that could be represented)
    GeomType GeomDim;

    //! The physical dimension
    unsigned char PhysicalDim;

    //! The physical dimension
    unsigned char d_max_gcw;

    //! The communicator over which the mesh is stored
    AMP_MPI d_comm;

    //! A pointer to an AMP database containing the mesh info
    boost::shared_ptr<AMP::Database>  d_db;

    //! A unique id for each mesh
    MeshID d_meshID;

    //! A name for the mesh
    std::string d_name;

    //! The bounding box for the mesh
    std::vector<double> d_box, d_box_local;

    /**
     *  A function to create a unique id for the mesh (requires the comm to be set)
     *  Note: this requires a global communication across the mesh communicator.
     *  Note: this function is NOT thread safe, and will need to be modified before threads are used.
     */
    void setMeshID();


};

} // Mesh namespace
} // AMP namespace

#endif

