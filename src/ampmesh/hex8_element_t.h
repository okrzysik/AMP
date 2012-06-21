#ifndef HEX8_ELEMENT_T
#define HEX8_ELEMENT_T

#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <cmath>

#include <ampmesh/triangle_t.h>

// deprecated
std::vector<double> get_basis_functions_values(const std::vector<double> &x);
std::vector<double> get_basis_functions_derivatives(const std::vector<double> &x);
//
void get_basis_functions_values(double const * const x, double * const basis_functions_values);
void get_basis_functions_derivatives( double const * const x, double * const basis_functions_derivatives);

class hex8_element_t {
public:
  hex8_element_t(const std::vector<double> &p);
  void set_support_points(const std::vector<double> &p);
  std::vector<double> get_support_points() const;
  std::vector<double> get_support_point(unsigned int i) const;
  std::vector<double> get_bounding_box();
  bool within_bounding_box(const std::vector<double> &p);
  bool within_bounding_polyhedron(const std::vector<double> &p);
  // this the user responsability to call first within_bounding_box(...) and within_bounding_polyhedron(...)
  // before trying to map a point so that newton won't fail
  std::vector<double> map_global_to_local(const std::vector<double> &global_coordinates);
  std::vector<double> map_local_to_global(const std::vector<double> &local_coordinates);
  void map_global_to_local(double const * const global_coordinates, double * const local_coordinates);
  void map_local_to_global(double const * const local_coordinates, double * const global_coordinates);
  bool contains_point(const std::vector<double> &coordinates, bool coordinates_are_local = false);
  std::pair<unsigned int, std::vector<double> > project_on_face(double a, double b, double c);
  std::pair<unsigned int, std::vector<double> > project_on_face(const std::vector<double> &p);
private:
  // numbering of the 8 support points (or nodes) follows libmesh hex8 convention which is as follows
  //
  //       7        6
  //        o--------o
  //       /:       /|
  //      / :      / |
  //   4 /  :   5 /  |
  //    o--------o   |
  //    |   o....|...o 2
  //    |  .3    |  /
  //    | .      | /
  //    |.       |/
  //    o--------o
  //    0        1
  //
  //
  // reference frame xyz 
  //    z   y
  //    |  .
  //    | .
  //    |.
  //    o------ x
  //
  //
  // node 0 -> x y z = -1 -1 -1
  //      1            +1 -1 -1
  //      2            +1 +1 -1
  //      3            -1 +1 -1
  //      4            -1 -1 +1
  //      5            +1 -1 +1
  //      6            +1 +1 +1
  //      7            -1 +1 +1
  //
  // numbering of the faces is
  // face 0 -> at z=-1 supported by nodes 0321 
  //      1       y=-1                    0154 
  //      2       x=+1                    1265 
  //      3       y=+1                    2376 
  //      4       x=-1                    3047 
  //      5       z=+1                    4567 
  //
  std::vector<double> support_points, point_candidate;
  bool bounding_box_updated, bounding_polyhedron_updated;
  std::vector<double> bounding_box;
  std::vector<triangle_t> bounding_polyhedron;
  bool center_of_element_data_updated;
  std::vector<double> center_of_element, inverse_jacobian_matrix_at_center_of_element;
  void compute_center_of_element_data();
  void build_bounding_box();
  void build_bounding_polyhedron();
  // residual vector x = (sum_i x_i b_i, sum_i y_i b_i, sum_i z_i b_i)^t
  // where the x_i y_i z_i, i=0...7, are the coordinates of the support points and the b_i are basis functions
  std::vector<double> compute_residual_vector(const std::vector<double> &x) const;
  // jacobian matrix J_ij = df_i / dj
  std::vector<double> compute_jacobian_matrix(const std::vector<double> &x) const;
  // A is a 3x3-matrix and the numbering of its elements is as follows
  // A[0] A[1] A[2] ^-1
  // A[3] A[4] A[5]
  // A[6] A[7] A[8]
  std::vector<double> compute_inverse_matrix(const std::vector<double> &A) const;
  // A is a 3x3-matrix and x is a column vector (3x1)
  // A[0] A[1] A[2]   x[0]   b[0]
  // A[3] A[4] A[5] * x[1] = b[1]
  // A[6] A[7] A[8]   x[2]   b[2]
  std::vector<double> compute_matrix_times_vector(const std::vector<double> &A, const std::vector<double> &x) const;
  std::vector<double> compute_inverse_jacobian_times_residual(const std::vector<double> &J, const std::vector<double> &f) const;
  std::vector<double> compute_initial_guess();
  // map the coordinates of the point candidate onto the reference frame of the volume element defined by the support points
  std::vector<double> solve_newton(double abs_tol = 1.0e-14, double rel_tol = 1.0e-14, unsigned int max_iter = 100, bool verbose = false);
};

#endif // HEX8_ELEMENT_T
