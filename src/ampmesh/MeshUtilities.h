#ifndef included_AMP_MeshUtilities
#define included_AMP_MeshUtilities


#include "AMP/ampmesh/Mesh.h"


namespace AMP::Mesh {


/**
 * \brief    Get points in the mesh
 * \details  This function returns points in the mesh at the given resolution.
 *     Note: The nodes are always used.
 * \param dx    Resolution to use.
 *              Note: if dx=0, then an automatic resolution is used.
 */
std::tuple<std::vector<Point>, std::vector<MeshElementID>> sample( const Mesh &mesh, double dx );


/**
 * \brief    Get points in the element
 * \details  This function returns points in the mesh element at the given resolution.
 * \param dx    Resolution to use.
 */
std::vector<Point> sample( const MeshElement &elem, double dx );

} // namespace AMP::Mesh

#endif
