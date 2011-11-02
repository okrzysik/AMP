#ifndef included_AMP_LibMesh
#define included_AMP_LibMesh

#include "ampmesh/Mesh.h"

// LibMesh include
#include "mesh.h"
#include "mesh_data.h"


namespace AMP {
namespace Mesh {



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

    /* Return the number of local element of the given type
     * \param type   Geometric type
     */
    virtual size_t  numLocalElements( GeomType &type );


    /* Return the global number of elements of the given type
     * \param type   Geometric type
     */
    virtual size_t  numTotalElements( GeomType &type );


    /**
     * \brief    Return an MeshIterator over the given geometric objects
     * \details  Return an MeshIterator over the given geometric objects
     * \param type   Geometric type to iterate over
     * \param gcw    Desired ghost 
     */
    virtual MeshIterator getIterator ( GeomType &type, int gcw=0 );

private:

    //!  Empty constructor for a mesh
    libMesh ( ) {};

    // libMesh objects
    boost::shared_ptr< ::Mesh>          d_libMesh;
    boost::shared_ptr< ::MeshData>      d_libMeshData;

};

} // Mesh namespace
} // AMP namespace

#endif

