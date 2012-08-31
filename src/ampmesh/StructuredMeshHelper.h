
#ifndef included_AMP_StructuredMeshHelper
#define included_AMP_StructuredMeshHelper

#include "utils/Utilities.h"
#include "ampmesh/Mesh.h"
#include "ampmesh/MeshElementVectorIterator.h"

namespace AMP{
namespace Mesh{


/**
 * \class StructuredMeshHelper
 * \brief A class used to provide utility functions for simple structured meshes
 * \details  This class provides utility routines when working with simple structured meshes
 *   in the x-y-z planes.  Most of the functions require the underlying mesh to have it's 
 *   nodes aligned in the x, y, and z directions.
 */
class StructuredMeshHelper 
{
public : 

    /**
     * \brief Get an iterator over the faces in the x-y planes
     * \details  For a simple mesh with nodes aligned in the x, y, z directions,
     *   this function returns an iterators that lie on the x-y planes 
     *   (all nodes for the given element share the same z-coordinate).
     * \param mesh  Mesh that we want to use for the iterator
     * \param gcw   Desired ghost width
     */
    static AMP::Mesh::MeshIterator getXYFaceIterator(AMP::Mesh::Mesh::shared_ptr mesh, int gcw=0);

    /**
     * \brief Get an iterator over the faces in the x-z planes
     * \details  For a simple mesh with nodes aligned in the x, y, z directions,
     *   this function returns an iterators that lie on the x-z planes 
     *   (all nodes for the given element share the same y-coordinate).
     * \param mesh  Mesh that we want to use for the iterator
     * \param gcw   Desired ghost width
     */
    static AMP::Mesh::MeshIterator getXZFaceIterator(AMP::Mesh::Mesh::shared_ptr mesh, int gcw=0);

    /**
     * \brief Get an iterator over the faces in the y-z planes
     * \details  For a simple mesh with nodes aligned in the x, y, z directions,
     *   this function returns an iterators that lie on the y-z planes 
     *   (all nodes for the given element share the same x-coordinate).
     * \param mesh  Mesh that we want to use for the iterator
     * \param gcw   Desired ghost width
     */
    static AMP::Mesh::MeshIterator getYZFaceIterator(AMP::Mesh::Mesh::shared_ptr mesh, int gcw);


    static AMP::Mesh::MeshIterator getGapFaceIterator(AMP::Mesh::Mesh::shared_ptr subChannel, int ghostWidth);

protected:

    static AMP::Mesh::MeshIterator getFaceIterator(AMP::Mesh::Mesh::shared_ptr mesh, int gcw, int direction);

};


}
}
#endif
