#ifndef included_AMP_StructuredMeshHelper
#define included_AMP_StructuredMeshHelper

#include "AMP/mesh/Mesh.h"


namespace AMP::Mesh {


/**
 * \class StructuredMeshHelper
 * \brief A class used to provide utility functions for simple structured meshes
 * \details  This class provides utility routines when working with simple structured meshes
 *   in the x-y-z planes.  Most of the functions require the underlying mesh to have it's
 *   nodes aligned in the x, y, and z directions.
 */
class StructuredMeshHelper
{
public:
    /**
     * \brief Get the x, y, z, coordinates of a cube mesh
     * \details  For a simple cube mesh with nodes aligned in the x, y, z directions,
     *   this function returns the list of the x, y, and z coordinates across all processors.
     *   This function will require communication on the mesh and will throw an error if
     *   the mesh nodes are not aligned
     * \param mesh  Mesh that we want to use for the iterator
     * \param x         x-coordinates
     * \param y         y-coordinates
     * \param z         z-coordinates
     * \param check     Check that the total number of nodes is unchanged
     */
    static void getXYZCoordinates( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                   std::vector<double> &x,
                                   std::vector<double> &y,
                                   std::vector<double> &z,
                                   bool check = true );

    /**
     * \brief Get an iterator over the faces in the x-y planes
     * \details  For a simple mesh with nodes aligned in the x, y, z directions,
     *   this function returns an iterators that lie on the x-y planes
     *   (all nodes for the given element share the same z-coordinate).
     * \param mesh  Mesh that we want to use for the iterator
     * \param gcw   Desired ghost width
     */
    static AMP::Mesh::MeshIterator getXYFaceIterator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                                      int gcw = 0 );

    /**
     * \brief Get an iterator over the faces in the x-z planes
     * \details  For a simple mesh with nodes aligned in the x, y, z directions,
     *   this function returns an iterators that lie on the x-z planes
     *   (all nodes for the given element share the same y-coordinate).
     * \param mesh  Mesh that we want to use for the iterator
     * \param gcw   Desired ghost width
     */
    static AMP::Mesh::MeshIterator getXZFaceIterator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                                      int gcw = 0 );

    /**
     * \brief Get an iterator over the faces in the y-z planes
     * \details  For a simple mesh with nodes aligned in the x, y, z directions,
     *   this function returns an iterators that lie on the y-z planes
     *   (all nodes for the given element share the same x-coordinate).
     * \param mesh  Mesh that we want to use for the iterator
     * \param gcw   Desired ghost width
     */
    static AMP::Mesh::MeshIterator getYZFaceIterator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                                      int gcw );

    /**
     * \brief Get an iterator over the faces
     * \details  For a simple mesh with nodes aligned in the x, y, z directions,
     *   this function returns an iterator of the faces that lie on the faces
     *   perpendicular to the given direction.
     * \param mesh      Mesh that we want to use for the iterator
     * \param gcw       Desired ghost width
     * \param direction Direction of faces to iterate
     */
    static AMP::Mesh::MeshIterator
    getFaceIterator( std::shared_ptr<AMP::Mesh::Mesh> mesh, int gcw, int direction );


    static AMP::Mesh::MeshIterator getGapFaceIterator( std::shared_ptr<AMP::Mesh::Mesh> subChannel,
                                                       int ghostWidth );

protected:
};
} // namespace AMP::Mesh
#endif
