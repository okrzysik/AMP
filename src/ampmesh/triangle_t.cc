#include <ampmesh/triangle_t.h>
#include <ampmesh/euclidean_geometry_tools.h>

#include <cassert>

triangle_t::triangle_t(double const * A, double const * B, double const * C) {
  set_support_points(A, B, C);
}

void triangle_t::set_support_points(double const * A, double const * B, double const * C) {
  if (support_points_ptr.size() == 0) { support_points_ptr.resize(3); }
  support_points_ptr[0] = A; 
  support_points_ptr[1] = B; 
  support_points_ptr[2] = C; 
  normal_updated = false;
  centroid_updated = false;
  edges_updated = false;
}

void triangle_t::set_support_points(double const * * ptr) {
  set_support_points(ptr[0], ptr[1], ptr[2]);
}

double const * triangle_t::get_support_point_ptr(unsigned int i) const {
  assert(i < 3);
  return support_points_ptr[i];
}

double const * triangle_t::get_normal() {
  if (!normal_updated) { compute_normal(); }
  return &(normal[0]);
}

double const * triangle_t::get_centroid() {
  if (!centroid_updated) { compute_centroid(); }
  return &(centroid[0]);
}

void triangle_t::compute_centroid() {
  if (centroid.size() == 0) { centroid.resize(3); }
  assert(!centroid_updated);
  for (unsigned int i = 0; i < 3; ++i) {
    centroid[i] = 0.0;
    for (unsigned int j = 0; j < 3; ++j) {
      centroid[i] += support_points_ptr[j][i];
    } // end for j
    centroid[i] /= 3.0;
  } // end for i
  centroid_updated = true;
}

void triangle_t::compute_normal() {
  if (normal.size() == 0) { normal.resize(3); tmp.resize(6); }
  assert(!normal_updated);
  make_vector_from_two_points(support_points_ptr[0], support_points_ptr[1], &(tmp[0])+0);
  make_vector_from_two_points(support_points_ptr[0], support_points_ptr[2], &(tmp[0])+3);
  compute_cross_product(&(tmp[0])+0, &(tmp[0])+3, &(normal[0]));
  normalize_vector(&(normal[0]));
  normal_updated= true;
}

// check whether plane supported by the triangle ABC is above the point with respect to the normal
bool triangle_t::above_point(double const * p, double tolerance) {
  if (!normal_updated) { 
    compute_normal();
  } // end if
  if (!centroid_updated) { 
    compute_centroid();
  } // end if
  make_vector_from_two_points(&(centroid[0]), p, &(tmp[0]));
  return (compute_scalar_product(&(tmp[0]), &(normal[0])) < tolerance);
}

edge_t * triangle_t::get_edge(unsigned int i) {
  assert(i < 3);
  if (!edges_updated) { build_edges(); }
  return &(edges[i]);
}
bool triangle_t::contains_point(double const * p, double tolerance) {
  if (!edges_updated) { build_edges(); }
  for (unsigned int i = 0; i < 3; ++i) {
    if (!edges[i].above_point(p, tolerance)) { 
      return false;
    } // end if
  } // end for i
  return true; 
}

void triangle_t::build_edges() {
  if (edges.size() == 0) { edges.reserve(3); }
  assert(!edges_updated);
  edges.clear();
  if (!normal_updated) { compute_normal(); }
  edges.push_back(edge_t(support_points_ptr[0], support_points_ptr[1], &(normal[0])));
  edges.push_back(edge_t(support_points_ptr[1], support_points_ptr[2], &(normal[0])));
  edges.push_back(edge_t(support_points_ptr[2], support_points_ptr[0], &(normal[0])));
  edges_updated = true;
}
