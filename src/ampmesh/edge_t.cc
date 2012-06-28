#include <vector>
#include <cassert>

#include "ampmesh/euclidean_geometry_tools.h"
#include "ampmesh/edge_t.h"

edge_t::edge_t(double const * A, double const * B, double const * ABC) {
  set_support_points(A, B);
  set_containing_plane(ABC);
}
void edge_t::set_support_points(double const * A, double const * B) {
  if (support_points_ptr.size() == 0) { support_points_ptr.resize(2); }
  support_points_ptr[0] = A;
  support_points_ptr[1] = B;
  center_updated = false;
  normal_updated = false;
}
void edge_t::set_containing_plane(double const * ABC) {
  containing_plane_ptr = ABC;
  normal_updated = false;
}
double const * edge_t::get_support_point_ptr(unsigned int i) const {
  assert(i < 2);
  return support_points_ptr[i];
}
double const * edge_t::get_normal() {
  if (!normal_updated) {
    compute_normal();
  } // end if
  return &(normal[0]);
}
double const * edge_t::get_center() {
  if (!center_updated) {
    compute_center();
  } // end if
  return &(center[0]);
}
bool edge_t::above_point(double const * p, double tolerance) {
  if (!normal_updated) {
    compute_normal();
  } // end if
  if (!center_updated) {
    compute_center();
  } // end if
  make_vector_from_two_points(&(center[0]), p, &(tmp[0]));
  return (compute_scalar_product(&(tmp[0]), &(normal[0])) < tolerance);
}

void edge_t::compute_normal() {
  if (normal.size() == 0) { normal.resize(3); tmp.resize(3); }
  assert(!normal_updated);
  make_vector_from_two_points(support_points_ptr[0], support_points_ptr[1], &(tmp[0])); 
  compute_cross_product(&(tmp[0]), containing_plane_ptr, &(normal[0]));
  normalize_vector(&(normal[0]));
  normal_updated = true;
}
void edge_t::compute_center() {
  if (center.size() == 0) { center.resize(3); }
  assert(!center_updated);
  for (unsigned int i = 0; i < 3; ++i) {
    center[i] = 0.0;
    for (unsigned int j = 0; j < 2; ++j) {
      center[i] += support_points_ptr[j][i];
    } // end for j
    center[i] /= 2.0;
  } // end for i
  center_updated = true;
}

