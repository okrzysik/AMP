
#include <ampmesh/euclidean_geometry_tools.h>
#include <ampmesh/edge_t.h>
#include <utils/Utilities.h>

#include <vector>
#include <cassert>
#include <cmath>


edge_t::edge_t(double const * A, double const * B, double const * ABC) {
  set_support_points(A, B);
  set_containing_plane(ABC);
}

void edge_t::set_support_points(double const * A, double const * B) {
  support_points_ptr[0] = A;
  support_points_ptr[1] = B;
  center_updated = false;
  direction_updated = false;
  normal_updated = false;
}

void edge_t::set_containing_plane(double const * ABC) {
  containing_plane_ptr = ABC;
  normal_updated = false;
}

double const * edge_t::get_support_point_ptr(unsigned int i) const {
  AMP_CHECK_ASSERT(i < 2);
  return support_points_ptr[i];
}

double const * edge_t::get_normal() {
  if (!normal_updated) { compute_normal(); }
  return &(normal[0]);
}

double const * edge_t::get_direction() {
  if (!direction_updated) { compute_direction(); }
  return &(direction[0]);
}

double const * edge_t::get_center() {
  if (!center_updated) { compute_center(); }
  return &(center[0]);
}

void edge_t::compute_normal() {
  AMP_CHECK_ASSERT(!normal_updated);
  if (normal.size() == 0) { normal.resize(3); }
  if (!direction_updated) { compute_direction(); }
  compute_cross_product(&(direction[0]), containing_plane_ptr, &(normal[0]));
  assert(fabs(compute_vector_norm(&(normal[0])) - 1.0) < 1.0e-14);
  //  normalize_vector(&(normal[0]));
  normal_updated = true;
}

void edge_t::compute_direction() {
  AMP_CHECK_ASSERT(!direction_updated);
  if (direction.size() == 0) { direction.resize(3); }
  make_vector_from_two_points(support_points_ptr[0], support_points_ptr[1], &(direction[0])); 
  normalize_vector(&(direction[0]));
  direction_updated = true;
}

void edge_t::compute_center() {
  AMP_CHECK_ASSERT(!center_updated);
  if (center.size() == 0) { center.resize(3); }
  for (unsigned int i = 0; i < 3; ++i) {
    center[i] = 0.0;
    for (unsigned int j = 0; j < 2; ++j) {
      center[i] += support_points_ptr[j][i];
    } // end for j
    center[i] /= 2.0;
  } // end for i
  center_updated = true;
}

double edge_t::compute_distance_to_containing_line(double const * point_in_containing_plane) {
  if (tmp.size() == 0) { tmp.resize(3); }
  if (!normal_updated) { compute_normal(); }
  if (!center_updated) { compute_center(); }
  // vector from center of the edge to point
  make_vector_from_two_points(&(center[0]), point_in_containing_plane, &(tmp[0]));
  // scalar product with normal to the edge in containing plane
  return compute_scalar_product(&(tmp[0]), &(normal[0]));
}

bool edge_t::above_point(double const * point, double tolerance) {
  double distance_to_containing_line = compute_distance_to_containing_line(point);
  return (distance_to_containing_line < tolerance);
}

int edge_t::project_point(double const * point_in_containing_plane, double * projection, double tolerance) {
  if (tmp.size() == 0) { tmp.resize(3); }
  double distance_to_containing_line = compute_distance_to_containing_line(point_in_containing_plane);

  // ensure that point was already projected onto containing plane
  make_vector_from_two_points(&(center[0]), point_in_containing_plane, &(tmp[0]));
  assert(compute_scalar_product(&(tmp[0]), containing_plane_ptr) < tolerance);

  if (distance_to_containing_line + tolerance > 0.0) {
    for (unsigned int i = 0; i < 3; ++i) {
      projection[i] = point_in_containing_plane[i] - distance_to_containing_line * normal[i];
    } // end for i

    make_vector_from_two_points(support_points_ptr[0], projection, &(tmp[0]));
    double position_on_containing_line = compute_scalar_product(&(tmp[0]), &(direction[0]));
    make_vector_from_two_points(support_points_ptr[0], support_points_ptr[1], &(tmp[0]));
    position_on_containing_line /= compute_vector_norm(&(tmp[0]));
    if (position_on_containing_line < 0.0 + tolerance) {
      for (unsigned int i = 0; i < 3; ++i) { projection[i] = support_points_ptr[0][i]; }
      return 0;
    } else if (position_on_containing_line > 1.0 - tolerance) {
      for (unsigned int i = 0; i < 3; ++i) { projection[i] = support_points_ptr[1][i]; }
      return 1;
    } // end if
    return 2;
  } // end if

  for (unsigned int i = 0; i < 3; ++i) { projection[i] = point_in_containing_plane[i]; }
  return -1;
}
