#include <ampmesh/triangle_t.h>
#include <ampmesh/euclidean_geometry_tools.h>

#include <cassert>
#include <iostream>

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
  if (normal.size() == 0) { normal.resize(3); }
  if (tmp.size() == 0) { tmp.resize(6); }
  assert(!normal_updated);
  make_vector_from_two_points(support_points_ptr[0], support_points_ptr[1], &(tmp[0])+0);
  make_vector_from_two_points(support_points_ptr[0], support_points_ptr[2], &(tmp[0])+3);
  compute_cross_product(&(tmp[0])+0, &(tmp[0])+3, &(normal[0]));
  normalize_vector(&(normal[0]));
  normal_updated= true;
}

double triangle_t::compute_distance_to_containing_plane(double const * point) {
  if (tmp.size()==0) { tmp.resize(6); }
  if (!normal_updated) { compute_normal(); }
  if (!centroid_updated) { compute_centroid(); }
  // vector from triangle centroid to point
  make_vector_from_two_points(&(centroid[0]), point, &(tmp[0]));
  // scalar product with normal to the triangle
  return compute_scalar_product(&(tmp[0]), &(normal[0]));
}

// check whether plane supported by the triangle ABC is above the point with respect to the normal
bool triangle_t::above_point(double const * point, double tolerance) {
  double distance_to_containing_plane = compute_distance_to_containing_plane(point);
  return (distance_to_containing_plane < tolerance);
}

edge_t * triangle_t::get_edge(unsigned int i) {
  assert(i < 3);
  if (!edges_updated) { build_edges(); }
  return &(edges[i]);
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

bool triangle_t::contains_point(double const * point, double tolerance) {
  if (!edges_updated) { build_edges(); }
  for (unsigned int i = 0; i < 3; ++i) {
    if (!edges[i].above_point(point, tolerance)) { 
      return false;
    } // end if
  } // end for i
  return true; 
}

bool triangle_t::project_point(double const * point, double * projection, double tolerance) {
  if (tmp.size()==0) { tmp.resize(6); }
  double distance_to_containing_plane = compute_distance_to_containing_plane(point);
  for (unsigned int i = 0; i < 3; ++i) { tmp[i] = point[i] - distance_to_containing_plane * normal[i]; }
  if (!edges_updated) { build_edges(); }
  for (unsigned int i = 0; i < 3; ) {
    int status = edges[i].project_point(&(tmp[0]), projection, tolerance);
    // -1 -> edge is above the point, we cannot conclude and go to the next edge
    if (status == -1) {
      ++i;
    // 2 -> projection onto the edge is normal, we have found the closest point on the triangle  
    } else if (status == 2) {
      return false;
    // 0 -> point was projected onto the first support point
    //      if we are on edge 0 we need to check edge 2
    //      otherwise we are done
    } else if (status == 0) {
      if (i == 0) { i = 2; continue; }
      return false;
    // 1 -> point was projected onto the second support point
    //      if we are on edge 2 we are done
    //      otherwise we check the next edge
    } else if (status == 1) {
      if (i == 2) { return false; }
      ++i;
    // just making sure nothing unexpected happened
    } else {
      std::cerr<<"how did you end up here in the first place?"<<std::endl;
      assert(false);
    } // end if
  } // end for i

  return true;
}
