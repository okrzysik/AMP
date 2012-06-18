#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <cmath>

std::vector<double> get_basis_functions_values(const std::vector<double> &x) {
  assert(x.size() == 3);
  std::vector<double> basis_functions_values(8);
  basis_functions_values[0] = 0.125*(1.0-x[0])*(1.0-x[1])*(1.0-x[2]);
  basis_functions_values[1] = 0.125*(1.0+x[0])*(1.0-x[1])*(1.0-x[2]);
  basis_functions_values[2] = 0.125*(1.0+x[0])*(1.0+x[1])*(1.0-x[2]);
  basis_functions_values[3] = 0.125*(1.0-x[0])*(1.0+x[1])*(1.0-x[2]);
  basis_functions_values[4] = 0.125*(1.0-x[0])*(1.0-x[1])*(1.0+x[2]);
  basis_functions_values[5] = 0.125*(1.0+x[0])*(1.0-x[1])*(1.0+x[2]);
  basis_functions_values[6] = 0.125*(1.0+x[0])*(1.0+x[1])*(1.0+x[2]);
  basis_functions_values[7] = 0.125*(1.0-x[0])*(1.0+x[1])*(1.0+x[2]);

  return basis_functions_values;
}

std::vector<double> get_basis_functions_derivatives(const std::vector<double> &x) {
  assert(x.size() == 3);
  std::vector<double> basis_functions_derivatives(24);
  basis_functions_derivatives[ 0] = 0.125*(-1.0)*(1.0-x[1])*(1.0-x[2]);
  basis_functions_derivatives[ 1] = 0.125*(+1.0)*(1.0-x[1])*(1.0-x[2]);
  basis_functions_derivatives[ 2] = 0.125*(+1.0)*(1.0+x[1])*(1.0-x[2]);
  basis_functions_derivatives[ 3] = 0.125*(-1.0)*(1.0+x[1])*(1.0-x[2]);
  basis_functions_derivatives[ 4] = 0.125*(-1.0)*(1.0-x[1])*(1.0+x[2]);
  basis_functions_derivatives[ 5] = 0.125*(+1.0)*(1.0-x[1])*(1.0+x[2]);
  basis_functions_derivatives[ 6] = 0.125*(+1.0)*(1.0+x[1])*(1.0+x[2]);
  basis_functions_derivatives[ 7] = 0.125*(-1.0)*(1.0+x[1])*(1.0+x[2]);

  basis_functions_derivatives[ 8] = 0.125*(1.0-x[0])*(-1.0)*(1.0-x[2]);
  basis_functions_derivatives[ 9] = 0.125*(1.0+x[0])*(-1.0)*(1.0-x[2]);
  basis_functions_derivatives[10] = 0.125*(1.0+x[0])*(+1.0)*(1.0-x[2]);
  basis_functions_derivatives[11] = 0.125*(1.0-x[0])*(+1.0)*(1.0-x[2]);
  basis_functions_derivatives[12] = 0.125*(1.0-x[0])*(-1.0)*(1.0+x[2]);
  basis_functions_derivatives[13] = 0.125*(1.0+x[0])*(-1.0)*(1.0+x[2]);
  basis_functions_derivatives[14] = 0.125*(1.0+x[0])*(+1.0)*(1.0+x[2]);
  basis_functions_derivatives[15] = 0.125*(1.0-x[0])*(+1.0)*(1.0+x[2]);

  basis_functions_derivatives[16] = 0.125*(1.0-x[0])*(1.0-x[1])*(-1.0);
  basis_functions_derivatives[17] = 0.125*(1.0+x[0])*(1.0-x[1])*(-1.0);
  basis_functions_derivatives[18] = 0.125*(1.0+x[0])*(1.0+x[1])*(-1.0);
  basis_functions_derivatives[19] = 0.125*(1.0-x[0])*(1.0+x[1])*(-1.0);
  basis_functions_derivatives[20] = 0.125*(1.0-x[0])*(1.0-x[1])*(+1.0);
  basis_functions_derivatives[21] = 0.125*(1.0+x[0])*(1.0-x[1])*(+1.0);
  basis_functions_derivatives[22] = 0.125*(1.0+x[0])*(1.0+x[1])*(+1.0);
  basis_functions_derivatives[23] = 0.125*(1.0-x[0])*(1.0+x[1])*(+1.0);

  return basis_functions_derivatives;
}

class hex8_element_t {
public:
  hex8_element_t(const std::vector<double> &p) { set_support_points(p); }

  void set_support_points(const std::vector<double> &p) { 
    assert(p.size() == 24); 
    support_points = p; 
    build_bounding_box();
  }

  std::vector<double> get_support_points() const { return support_points; } 

  std::vector<double> get_bounding_box() const { return bounding_box; }

  bool within_bounding_box(const std::vector<double> &p) {
    assert(p.size() == 3);
    for (unsigned int j = 0; j < 3; ++j) {
      if ((bounding_box[j+0] > p[j]) || (bounding_box[j+3] < p[j])) { return false; }
    } // end for j
    return true;
  }

  std::vector<double> map_global_to_local(const std::vector<double> &global_coordinates) {
    assert(global_coordinates.size() == 3); 
    point_candidate = global_coordinates;
    std::vector<double> local_coordinates = solve_newton();
    return local_coordinates;
  }

  std::vector<double> map_local_to_global(const std::vector<double> &local_coordinates) {
    assert(local_coordinates.size() == 3); 
    point_candidate = std::vector<double>(3, 0.0);
    std::vector<double> global_coordinates = compute_residual_vector(local_coordinates);
    for (unsigned int i = 0; i < 3; ++i) { global_coordinates[i] *= -1.0; }
    return global_coordinates;
  }

  bool contains_point(const std::vector<double> &coordinates, bool coordinates_are_local = false) {
    assert(coordinates.size() == 3); 
    std::vector<double> local_coordinates;
    if (!coordinates_are_local) {
      point_candidate = coordinates;
      if (!within_bounding_box(point_candidate)) { return false; }
      local_coordinates = solve_newton();
    } else {
      local_coordinates = coordinates;
    } // end if
    for (unsigned int i = 0; i < 3; ++i) {
      if (fabs(local_coordinates[i]) > 1.0) {
        return false;
      } // end if
    } // end for i
    return true;
  }

  std::pair<unsigned int, std::vector<double> > project_on_face(double a, double b, double c) {
    double p[3] = { a, b, c };
    return project_on_face(std::vector<double>(p, p+3));
  }

  std::pair<unsigned int, std::vector<double> > project_on_face(const std::vector<double> &p) {
    assert(p.size() == 3); 
    if (!within_bounding_box(p)) { return std::pair<unsigned int, std::vector<double> >(99, std::vector<double>(2, 0.0)); }
    point_candidate = p;

    // compute the coordinates of the point candidate in the frame of the volume element
    std::vector<double> x = solve_newton();

    double distance_to_face[6];
    distance_to_face[0] = 1.0 + x[2]; 
    distance_to_face[1] = 1.0 + x[1]; 
    distance_to_face[2] = 1.0 - x[0]; 
    distance_to_face[3] = 1.0 - x[1]; 
    distance_to_face[4] = 1.0 + x[0]; 
    distance_to_face[5] = 1.0 - x[2]; 

    double min_distance = 1.0;
    unsigned int closest_face = 33;
    for (unsigned int i = 0; i < 6; ++i) {
      if (distance_to_face[i] < 0.0) {
        closest_face = 99;
        break;
      } else if (distance_to_face[i] < min_distance) {
        min_distance = distance_to_face[i];
        closest_face = i;
      } // end if
    } // end for i

    // faces are oriented and defined by their 4 support nodes
    //   3          2    
    //    o--------o    
    //    |        |
    //    |        |   
    //    |        |  
    //    |        | 
    //    o--------o
    //   0          1
    //
    // coordinates on the face are given in the following xy reference frame
    // 
    //    y
    //    |
    //    |
    //    |
    //    o------ x
    //
    std::vector<double> coordinates_on_face(2, 0.0);
    if (closest_face == 0) {
      coordinates_on_face[0] = x[1];
      coordinates_on_face[1] = x[0];
    } else if (closest_face == 1) {
      coordinates_on_face[0] = x[0];
      coordinates_on_face[1] = x[2];
    } else if (closest_face == 2) {
      coordinates_on_face[0] = x[1];
      coordinates_on_face[1] = x[2];
    } else if (closest_face == 3) {
      coordinates_on_face[0] = -x[0];
      coordinates_on_face[1] = x[2];
    } else if (closest_face == 4) {
      coordinates_on_face[0] = -x[1];
      coordinates_on_face[1] = x[2];
    } else if (closest_face == 5) {
      coordinates_on_face[0] = x[0];
      coordinates_on_face[1] = x[1];
    } // end if

    assert(closest_face != 33);
    return std::pair<unsigned int, std::vector<double> >(closest_face, coordinates_on_face);
  }
  
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
  std::vector<double> support_points, point_candidate, bounding_box;

  void build_bounding_box() {
    bounding_box = std::vector<double>(6, 0.0);
    for (unsigned int j = 0; j < 3; ++j) {
      bounding_box[j+0] = support_points[3*0+j];
      bounding_box[j+3] = support_points[3*0+j];
    } // end for j
    for (unsigned int i = 1; i < 8; ++i) {
      for (unsigned int j = 0; j < 3; ++j) {
        if (bounding_box[j+0] > support_points[3*i+j]) {
          bounding_box[j+0] = support_points[3*i+j];
        } // end if
        if (bounding_box[j+3] < support_points[3*i+j]) {
          bounding_box[j+3] = support_points[3*i+j];
        } // end if
      } // end for j
    } // end for i
  }

  // residual vector x = (sum_i x_i b_i, sum_i y_i b_i, sum_i z_i b_i)^t
  // where the x_i y_i z_i, i=0...7, are the coordinates of the support points and the b_i are basis functions
  std::vector<double> compute_residual_vector(const std::vector<double> &x) {
    assert(x.size() == 3);
    std::vector<double> f(3, 0.0);

    // basis functions i is one at node i and zero at node j not i
    std::vector<double> basis_functions_values = get_basis_functions_values(x);

    for (unsigned int i = 0; i < 3; ++i) {
      f[i] = point_candidate[i]; 
      for (unsigned int j = 0; j < 8; ++j) {
        f[i] -= support_points[3*j+i] * basis_functions_values[j];
      } // end for j
    } // end for i

    return f;
  }

  // jacobian matrix J_ij = df_i / dj
  std::vector<double> compute_jacobian_matrix(const std::vector<double> &x) {
    assert(x.size() == 3);
    std::vector<double> J(9, 0.0);

    // derivatives of the basis function i with respect to j stored in i+8*j
    std::vector<double> basis_functions_derivatives = get_basis_functions_derivatives(x);

    for (unsigned int i = 0; i < 3; ++i) {
      for (unsigned int j = 0; j < 3; ++j) { 
        for (unsigned int k = 0; k < 8; ++k) {
          J[3*i+j] -= support_points[3*k+i] * basis_functions_derivatives[8*j+k];
        } // end for k
      } // end for j
    } // end for i
  
    return J;
  }

  // A is a 3x3-matrix and the numbering of its elements is as follows
  // A[0] A[1] A[2] ^-1
  // A[3] A[4] A[5]
  // A[6] A[7] A[8]
  std::vector<double> compute_inverse_matrix(const std::vector<double> &A) {
    assert(A.size() == 9);
    double determinant = A[0]*(A[4]*A[8]-A[5]*A[7])-A[1]*(A[3]*A[8]-A[5]*A[6])+A[2]*(A[3]*A[7]-A[4]*A[6]);
    assert(fabs(determinant) > 1.0e-15);
    double one_over_determinant = 1.0 / determinant;
    std::vector<double> inverse_of_A(9, 0.0);
    inverse_of_A[0] = one_over_determinant*(A[4]*A[8]-A[5]*A[7]);
    inverse_of_A[1] = one_over_determinant*(A[2]*A[7]-A[1]*A[8]);
    inverse_of_A[2] = one_over_determinant*(A[1]*A[5]-A[2]*A[4]);
    inverse_of_A[3] = one_over_determinant*(A[5]*A[6]-A[3]*A[8]);
    inverse_of_A[4] = one_over_determinant*(A[0]*A[8]-A[2]*A[6]);
    inverse_of_A[5] = one_over_determinant*(A[2]*A[3]-A[0]*A[5]);
    inverse_of_A[6] = one_over_determinant*(A[3]*A[7]-A[4]*A[6]);
    inverse_of_A[7] = one_over_determinant*(A[1]*A[6]-A[0]*A[7]);
    inverse_of_A[8] = one_over_determinant*(A[0]*A[4]-A[1]*A[3]);
    return inverse_of_A;
  }

  // A is a 3x3-matrix and x is a column vector (3x1)
  // A[0] A[1] A[2]   x[0]   b[0]
  // A[3] A[4] A[5] * x[1] = b[1]
  // A[6] A[7] A[8]   x[2]   b[2]
  std::vector<double> compute_matrix_times_vector(const std::vector<double> &A, const std::vector<double> &x) {
    assert(A.size() == 9);
    assert(x.size() == 3);
    std::vector<double> b(3, 0.0);
    b[0] = A[0]*x[0]+A[1]*x[1]+A[2]*x[2];
    b[1] = A[3]*x[0]+A[4]*x[1]+A[5]*x[2];
    b[2] = A[6]*x[0]+A[7]*x[1]+A[8]*x[2];
    return b;
  }

  std::vector<double> compute_inverse_jacobian_times_residual(const std::vector<double> &J, const std::vector<double> &f) {
    assert(J.size() == 9);
    assert(f.size() == 3);
    return compute_matrix_times_vector(compute_inverse_matrix(J), f);
  }

  // map the coordinates of the point candidate onto the reference frame of the volume element defined by the support points
  std::vector<double> solve_newton(double abs_tol = 1.0e-14, double rel_tol = 1.0e-14, unsigned int max_iter = 100, bool verbose = false) {
    if (verbose) { std::cout<<"solve newton with line search\n"; }
    std::vector<double> x(3, 0.0);
    std::vector<double> residual_vector = compute_residual_vector(x);
    double residual_norm = sqrt(inner_product(residual_vector.begin(), residual_vector.end(), residual_vector.begin(), 0.0));
    double tol = abs_tol + rel_tol * residual_norm; 
    for (unsigned int iter = 0; iter < max_iter; ++iter) {
      std::vector<double> inverse_jacobian_matrix_times_residual_vector = compute_inverse_jacobian_times_residual(compute_jacobian_matrix(x), residual_vector);
      for (double alpha = 1.0; alpha > 1.0e-4; alpha /= 1.1) {
        std::vector<double> tmp(3);
        for (unsigned int i = 0; i < 3; ++i) { tmp[i] = x[i] - alpha * inverse_jacobian_matrix_times_residual_vector[i]; };
        residual_vector = compute_residual_vector(tmp);
        double tmp_residual_norm = sqrt(inner_product(residual_vector.begin(), residual_vector.end(), residual_vector.begin(), 0.0));
        if (tmp_residual_norm < residual_norm) {
          x = tmp;
          residual_norm = tmp_residual_norm;
          break;
        } // end if
      } // end for
      if (verbose) { std::cout<<iter<<"  "<<residual_norm<<std::endl; }
      if (residual_norm < tol) { 
        if (verbose) { std::cout<<"converged at iteration "<<iter<<" with residual norm "<<residual_norm<<"\n"; }
        return x; 
      } // end if
    } // end for
    std::cerr<<"failed to converge with tolerance "<<tol<<" after "<<max_iter-1<<" iterations (residual norm was "<<residual_norm<<")"<<std::endl;
    abort(); 
  }
};

void scale_points(unsigned int direction, double scaling_factor, unsigned int n_points, double* points) {
  assert(direction < 3);
  for (unsigned int i = 0; i < n_points; ++i) { 
    points[3*i+direction] *= scaling_factor; 
  } // end for i
}

void scale_points(const std::vector<double> &scaling_factors, unsigned int n_points, double* points) {
  assert(scaling_factors.size() == 3);
  for (unsigned int i = 0; i < 3; ++i) { 
    scale_points(i, scaling_factors[i], n_points, points);
  } // end for i
}

void translate_points(unsigned int direction, double distance, unsigned int n_points, double* points) {
  assert(direction < 3);
  for (unsigned int i = 0; i < n_points; ++i) { 
    points[3*i+direction] += distance; 
  } // end for i
}

void translate_points(std::vector<double> translation_vector, unsigned int n_points, double* points) {
  assert(translation_vector.size() == 3);
  for (unsigned int i = 0; i < 3; ++i) { 
   translate_points(i, translation_vector[i], n_points, points); 
  } // end for i
}

void rotate_points(unsigned int rotation_axis, double rotation_angle, unsigned int n_points, double* points) {
  assert(rotation_axis < 3);
  unsigned int non_fixed_directions[2];
  unsigned int i = 0;
  for (unsigned int j = 0; j < 3; ++j) { 
    if (j != rotation_axis) {
      non_fixed_directions[i++] = j;
    } // end if
  } // end for j
  double tmp[3];
  for (unsigned int j = 0; j < n_points; ++j) {
    tmp[non_fixed_directions[0]] = cos(rotation_angle)*points[3*j+non_fixed_directions[0]]-sin(rotation_angle)*points[3*j+non_fixed_directions[1]];
    tmp[non_fixed_directions[1]] = sin(rotation_angle)*points[3*j+non_fixed_directions[0]]+cos(rotation_angle)*points[3*j+non_fixed_directions[1]];
    tmp[rotation_axis] = points[3*j+rotation_axis];
    std::copy(tmp, tmp+3, points+3*j);
  } // end for j
}

