#ifndef included_AMP_MeshUtilities
#define included_AMP_MeshUtilities


#include "AMP/ampmesh/Geometry.h"
#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MeshPoint.h"
#include "AMP/utils/Array.h"
#include "AMP/utils/kdtree.h"


namespace AMP::Mesh {


/**
 * \brief    Get points in the mesh
 * \details  This function returns points in the mesh at the given resolution.
 *     Note: The nodes are always used.
 * \param mesh  Mesh to sample
 * \param dx    Resolution to use
 *              Note: if dx=0, then an automatic resolution is used.
 */
std::tuple<std::vector<Point>, std::vector<MeshElementID>> sample( const Mesh &mesh, double dx );


/**
 * \brief    Get points in the element
 * \details  This function returns points in the mesh element at the given resolution.
 * \param elem  Element to sample
 * \param dx    Resolution to use.
 */
std::vector<Point> sample( const MeshElement &elem, double dx );


/**
 * \brief    Get a grid with the overlap with the given geometry
 * \details  This function returns an Array with the given resolution in which
 *     each cell contains the volume of the underlying geometry that intersects with the cell.
 * \param geom  Geometry to use
 * \param N     Number of cells to use
 */
Array<double> volumeOverlap( const AMP::Geometry::Geometry &geom, const std::vector<int> &N );


/**
 * \class ElementFinder
 * \brief A class used to help find the nearest element to a point
 * \details  This class provides functions to help find the nearest element to a point
 */
class ElementFinder final
{
public:
    //! Empty constructor
    ElementFinder() : d_pos_hash( 0 ), d_type( AMP::Mesh::GeomType::Vertex ) {}

    //! Empty constructor
    ElementFinder(
        std::shared_ptr<AMP::Mesh::Mesh> mesh,
        std::vector<AMP::Mesh::MeshElementID> elements = std::vector<AMP::Mesh::MeshElementID>() );

    //! Copy constructor
    ElementFinder( const ElementFinder & ) = delete;

    //! Move constructor
    ElementFinder( ElementFinder && ) = default;

    //! Move operator
    ElementFinder &operator=( ElementFinder && ) = default;

    //! Get the nearest elements to a point
    std::vector<AMP::Mesh::MeshElement> getNearestElements( const Point &x ) const;

    //! Get the nearest element and point
    std::pair<AMP::Mesh::MeshElement, Point> getNearestPoint( const Point &x ) const;

private:
    void initialize() const;

private:
    std::shared_ptr<AMP::Mesh::Mesh> d_mesh;
    mutable uint64_t d_pos_hash;
    AMP::Mesh::GeomType d_type;
    mutable kdtree d_tree;
    mutable std::vector<AMP::Mesh::MeshElementID> d_ids;
    std::vector<AMP::Mesh::MeshElementID> d_elements;
};


} // namespace AMP::Mesh


#endif
