#ifndef TRIANGLE_T
#define TRIANGLE_T

#include <cassert>
#include <cmath>
#include <numeric>
#include <vector>

void scale_points(unsigned int direction, double scaling_factor, unsigned int n_points, double* points);
void scale_points(const std::vector<double> &scaling_factors, unsigned int n_points, double* points);
void translate_points(unsigned int direction, double distance, unsigned int n_points, double* points);
void translate_points(std::vector<double> translation_vector, unsigned int n_points, double* points);
void rotate_points(unsigned int rotation_axis, double rotation_angle, unsigned int n_points, double* points);

std::vector<double> compute_cross_product(const std::vector<double> &u, const std::vector<double> &v);
double compute_scalar_product(const std::vector<double> &u, const std::vector<double> &v);
std::vector<double> make_vector_from_two_points(const std::vector<double> &start_point, const std::vector<double> &end_point);
void compute_cross_product(double const * u, double const * v, double * w);
double compute_scalar_product(double const * u, double const * const v);
void make_vector_from_two_points(double const * start_point, double const * end_point, double * vector);
void normalize_vector(double * vector);

class edge_t {
public:
  edge_t(double const * A, double const * B, double const * ABC) {
    set_support_points(A, B);
    set_containing_plane(ABC);
  }
  void set_support_points(double const * A, double const * B) {
    if (support_points_ptr.size() == 0) { support_points_ptr.resize(2); }
    support_points_ptr[0] = A;
    support_points_ptr[1] = B;
    center_updated = false;
    normal_updated = false;
  }
  void set_containing_plane(double const * ABC) {
    containing_plane_ptr = ABC;
    normal_updated = false;
  }
  double const * get_support_point_ptr(unsigned int i) const {
    assert(i < 2);
    return support_points_ptr[i];
  }
  double const * get_normal() {
    if (!normal_updated) {
      compute_normal();
    } // end if
    return &(normal[0]);
  }
  double const * get_center() {
    if (!center_updated) {
      compute_center();
    } // end if
    return &(center[0]);
  }
  bool above_point(double const * p, double tolerance = 1.0e-12) {
    if (!normal_updated) {
      compute_normal();
    } // end if
    if (!center_updated) {
      compute_center();
    } // end if
    make_vector_from_two_points(&(center[0]), p, &(tmp[0]));
    return (compute_scalar_product(&(tmp[0]), &(normal[0])) < tolerance);
  }

private:
  std::vector<double const *> support_points_ptr;
  double const * containing_plane_ptr;
  std::vector<double> normal, center, tmp;
  bool normal_updated, center_updated;
  void compute_normal() {
    if (normal.size() == 0) { normal.resize(3); tmp.resize(3); }
    assert(!normal_updated);
    make_vector_from_two_points(support_points_ptr[0], support_points_ptr[1], &(tmp[0])); 
    compute_cross_product(&(tmp[0]), containing_plane_ptr, &(normal[0]));
    normalize_vector(&(normal[0]));
    normal_updated = true;
  }
  void compute_center() {
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
};

class triangle_t {
public:
  triangle_t(double const * A, double const * B, double const * C);
  void set_support_points(double const * A, double const * B, double const * C);
  void set_support_points(double const * * ptr);
  bool above_point(double const * p, double tolerance = 1.0e-12);
  double const * get_support_point_ptr(unsigned int i) const;
  double const * get_normal();
  double const * get_centroid();
  edge_t * get_edge(unsigned int i) {
    assert(i < 3);
    if (!edges_updated) { build_edges(); }
    return &(edges[i]);
  }
  bool contains_point(double const * p, double tolerance = 1.0e-12) {
    if (!edges_updated) { build_edges(); }
    for (unsigned int i = 0; i < 3; ++i) {
      if (!edges[i].above_point(p, tolerance)) { 
        return false;
      } // end if
    } // end for i
    return true; 
  }

private:
  std::vector<double const *> support_points_ptr;
  std::vector<double> normal, centroid, tmp;
  bool normal_updated, centroid_updated;
  void compute_normal();
  void compute_centroid();

  std::vector<edge_t> edges;
  bool edges_updated;
  void build_edges() {
    if (edges.size() == 0) { edges.reserve(3); }
    assert(!edges_updated);
    edges.clear();
    if (!normal_updated) { compute_normal(); }
    edges.push_back(edge_t(support_points_ptr[0], support_points_ptr[1], &(normal[0])));
    edges.push_back(edge_t(support_points_ptr[1], support_points_ptr[2], &(normal[0])));
    edges.push_back(edge_t(support_points_ptr[2], support_points_ptr[0], &(normal[0])));
    edges_updated = true;
  }
};

#endif // TRIANGLE_T
