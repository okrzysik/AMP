#ifndef included_AMP_Mesh
#define included_AMP_Mesh

#include "MeshIterator.h"
#include "MeshElement.h"


namespace AMP {
namespace Mesh {


//! Enumeration for basic mesh-based quantities
enum SetOP { Union, Intersection, Complement };


/**
 * \class Mesh
 * \brief A class used to abstract away mesh from an application
 *
 * \details  This class provides routines for reading, accessing and writing meshes.
 */
class Mesh
{
public:
    /**
     *\typedef shared_ptr
     *\brief  Name for the shared pointer.
     *\details  Use this typedef for a reference counted pointer to a mesh manager object.
     */
    typedef boost::shared_ptr<Mesh>  shared_ptr;

    /**
     * \param params Parameters for constructing a mesh from an input database
     * \brief Read in mesh files, partition domain, and prepare environment for simulation
     * \details  For trivial parallelsim, this method reads in the meshes on each processor.  Each
     * processor contains a piece of each mesh.  For massive parallelism, each mesh is on its own
     * communicator.  As such, some math libraries must be initialized accordingly.
     */
    virtual Mesh ( const MeshParameters::shared_ptr &params );

    /**
     * \brief Construct a new mesh from an existing mesh
     * \details  This constructor will construct a new mesh from an existing mesh.
     * This is designed as a path to create a new mesh object of one type from
     * an existing mesh of a different type.  It also allows creating a new single mesh
     * from a subset or superset of other meshes.  Note that instantion of this routine 
     * may not be able to create it's mesh from any arbitrary mesh, and may throw an 
     * error.
     * \param old_mesh Existing mesh that we will use to construct the new mesh
     */
    virtual Mesh ( const Mesh::shared_ptr &old_mesh );

    /**
     * \brief Construct a new mesh from an existing mesh.
     * \details  This constructor will construct a new mesh from an existing mesh
     * using an iterator over the existing mesh.
     * This is designed as a path to create a new mesh object of one type from
     * an existing mesh of a different type.  It also allows creating a new single mesh
     * from a subset or superset of other meshes.  Note that instantion of this routine 
     * may not be able to create it's mesh from any arbitrary mesh, and may throw an 
     * error.
     * \param old_mesh Existing mesh that we will use to construct the new mesh
     */
    virtual Mesh ( const Mesh::shared_ptr &old_mesh, MeshIterator::shared_ptr &iterator);

    /**
     * \brief    Subset a mesh given a MeshIterator
     * \details  This function will subset a mesh over a given iterator.
     *   This will return a new mesh object.
     * \param iterator  MeshIterator used to subset
     */
    virtual boost::shared_ptr<Mesh>  Subset ( MeshIterator::shared_ptr &iterator );

    /**
     * \brief        Subset a mesh given another mesh
     * \details      This function will subset a mesh given another mesh
     * \param mesh   Mesh used to subset
     */
    virtual boost::shared_ptr<Mesh>  Subset ( Mesh &mesh );

    /**
     * \brief    Return an MeshIterator over the given geometric objects
     * \details  Return an MeshIterator over the given geometric objects
     * \param type   Geometric type to iterate over
     * \param gcw    Desired ghost 
     */
    virtual MeshIterator getIterator ( GeomType &type, int gcw=0 );

    /**
     * \brief    Return an MeshIterator over the given geometric objects on the surface
     * \details  Return an MeshIterator over the given geometric objects on the surface
     * \param type   Geometric type to iterate over
     */
    virtual MeshIterator getSurfaceIterator ( GeomType &type );

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
    virtual MeshIterator getIterator ( SetOP &OP, MeshIterator::shared_ptr &A, MeshIterator::shared_ptr &B);
 

protected:

    //!  Empty constructor for a mesh
    Mesh ( );

    //! The mesh parameters
    const MeshParameters::shared_ptr &params;

    //! The geometric dimension (equivalent to the highest geometric object that could be represented)
    const GeomType &GeomDim;

    //! The communicator over which the mesh is stored
    AMP_MPI &comm;

};

} // Mesh namespace
} // AMP namespace

#endif

