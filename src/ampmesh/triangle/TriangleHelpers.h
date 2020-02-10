#ifndef included_AMP_TriangleMeshHelpers
#define included_AMP_TriangleMeshHelpers

#include "AMP/ampmesh/Geometry.h"
#include "AMP/ampmesh/Mesh.h"
#include "AMP/utils/AMP_MPI.h"

#include <array>
#include <memory>
#include <vector>


namespace AMP::Mesh::TriangleHelpers {


//! Count the number of unique triangles
template<size_t NG>
size_t count( const std::vector<std::array<int64_t, NG + 1>> &tri );


//! Read an STL file
std::vector<std::array<std::array<double, 3>, 3>> readSTL( const std::string &filename,
                                                           double scale );

//! Read the header for an STL file
size_t readSTLHeader( const std::string &filename );


//! Create triangles/verticies from a set of triangles specified by their coordinates
template<size_t NG, size_t NP>
void createTriangles( const std::vector<std::array<std::array<double, NP>, NG + 1>> &tri_list,
                      std::vector<std::array<double, NP>> &verticies,
                      std::vector<std::array<int64_t, NG + 1>> &triangles,
                      double tol );

//! Create triangles neighbors from the triangles
template<size_t NG>
std::vector<std::array<int64_t, NG + 1>>
create_tri_neighbors( const std::vector<std::array<int64_t, NG + 1>> &tri );

//! Create triangles neighbors from the triangles
template<size_t NG>
std::vector<std::vector<std::array<int64_t, NG + 1>>>
    splitDomains( std::vector<std::array<int64_t, NG + 1>> tri );

//! Read an STL file and generate a mesh (triangle mesh or multi-mesh)
std::shared_ptr<AMP::Mesh::Mesh> generateSTL( std::shared_ptr<AMP::Mesh::MeshParameters> params );

//! Generate a triangle mesh (or multi-mesh) from a geometry
std::shared_ptr<AMP::Mesh::Mesh>
generate( std::shared_ptr<AMP::Geometry::Geometry> geom, const AMP_MPI &comm, double resolution );

} // namespace AMP::Mesh::TriangleHelpers

#endif
