#ifndef included_AMP_moabMesh
#define included_AMP_moabMesh

#include "AMP/ampmesh/Mesh.h"

// MOAB include
#include "moab/Core.hpp"


namespace AMP {
namespace Mesh {


class moabMeshElement;


/**
 * \class moabMesh
 * \brief A concrete mesh class for moabMesh
 *
 * \details  This class provides routines for reading, accessing and writing moabMesh meshes.
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
 * The parallel decomposition of the mesh is controlled by moabMesh and occurs on the communicator
 * specified through the params object.  Note that moabMesh does not support meshes on overlapping
 * communicators.  If multiple meshes are used, they must either share communicators or have unique
 * communicators.
 */
class moabMesh : public Mesh
{
public:
    /**
     * \brief Read in mesh files, partition domain, and prepare environment for simulation
     * \details  For trivial parallelsim, this method reads in the meshes on each processor.  Each
     * processor contains a piece of each mesh.  For massive parallelism, each mesh is on its own
     * communicator.  As such, some math libraries must be initialized accordingly.
     * \param params Parameters for constructing a mesh from an input database
     */
    moabMesh( const std::shared_ptr<MeshParameters> &params );


    //! Deconstructor
    virtual ~moabMesh();


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
    static size_t estimateMeshSize( const std::shared_ptr<MeshParameters> &params );


    /* Return the number of local element of the given type
     * \param type   Geometric type
     */
    size_t numLocalElements( const GeomType type ) const override;


    /* Return the global number of elements of the given type
     * \param type   Geometric type
     */
    size_t numGlobalElements( const GeomType type ) const override;


    /* Return the number of ghost elements of the given type on the current processor
     * \param type   Geometric type
     * \param gcw    Desired ghost cell width
     */
    size_t numGhostElements( const GeomType type, const int gcw ) const override;


    /**
     * \brief    Return an MeshIterator over the given geometric objects
     * \details  Return an MeshIterator over the given geometric objects
     * \param type   Geometric type to iterate over
     * \param gcw    Desired ghost cell width
     */
    MeshIterator getIterator( const GeomType type, const int gcw = 0 ) const override;


    /**
     * \brief    Return an MeshIterator over the given geometric objects on the surface
     * \details  Return an MeshIterator over the given geometric objects on the surface
     * \param type   Geometric type to iterate over
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator getSurfaceIterator( const GeomType type,
                                             const int gcw = 0 ) const override;


    //! Check if two meshes are equal
    bool operator==( const Mesh &mesh ) const override;


    /**
     * \brief    Return the list of all ID sets in the mesh
     * \details  Return the list of all ID sets in the mesh
     */
    std::vector<int> getBoundaryIDs() const override;


    /**
     * \brief    Return an MeshIterator over the given geometric objects on the given ID set
     * \details  Return an MeshIterator over the given geometric objects on the given ID set
     * \param type   Geometric type to iterate over
     * \param id     id for the elements (example: nodeset id)
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator
    getBoundaryIDIterator( const GeomType type, const int id, const int gcw = 0 ) const override;


    /**
     * \brief    Displace the entire mesh
     * \details  This function will displace the entire mesh by a scalar value.
     *   This function is a blocking call for the mesh communicator, and requires
     *   the same value on all processors.  The displacement vector should be the
     *   size of the physical dimension.
     * \param x  Displacement vector
     */
    void displaceMesh( std::vector<double> x ) override;


#ifdef USE_AMP_VECTORS
    /**
     * \brief    Displace the entire mesh
     * \details  This function will displace the entire mesh by displacing
     *   each node by the values provided in the vector.  This function is
     *   a blocking call for the mesh communicator
     * \param x  Displacement vector.  Must have N DOFs per node where N
     *           is the physical dimension of the mesh.
     */
    void displaceMesh( std::shared_ptr<const AMP::LinearAlgebra::Vector> x ) override;
#endif


protected:
private:
    //!  Empty constructor for a mesh
    moabMesh(){};

    //!  Function to properly initialize the internal data once a moabMesh mesh is loaded
    void initialize();

    //  Internal variables
    std::shared_ptr<moab::Core> d_core;
};

} // namespace Mesh
} // namespace AMP

#endif
