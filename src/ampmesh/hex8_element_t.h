#ifndef HEX8_ELEMENT_T
#define HEX8_ELEMENT_T

#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <cmath>

#include "triangle_t.h"

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

  std::vector<double> get_support_point(unsigned int i) const { 
    assert(i < 8);
    return std::vector<double>(&support_points[3*i], &support_points[3*i]+3);
  } 

  std::vector<double> get_bounding_box() const { return bounding_box; }

  bool within_bounding_box(const std::vector<double> &p) {
    assert(p.size() == 3);
    for (unsigned int j = 0; j < 3; ++j) {
      if ((bounding_box[j+0] > p[j]) || (bounding_box[j+3] < p[j])) { return false; }
    } // end for j
    return true;
  }

  bool contained_by_triangles_on_faces(const std::vector<double> &p) {
    unsigned int faces[24] = {
      0, 3, 2, 1, 
      0, 1, 5, 4, 
      1, 2, 6, 5, 
      2, 3, 7, 6, 
      3, 0, 4, 7, 
      4, 5, 6, 7
    };
    std::vector<triangle_t> triangles(4);
    for (unsigned int i = 0; i < 6; ++i) {
      // first configuration when splitting face into two triangles
      //   3          2    
      //    o--------o    
      //    | .      |
      //    |   .    |   
      //    |     .  |  
      //    |       .| 
      //    o--------o
      //   0          1
      //
      //   3        
      //    o
      //    | .
      //    |   .
      //    |     .
      //    |       .
      //    o--------o
      //   0          1
      triangles[0].set_support_points(get_support_point(faces[4*i+0]), get_support_point(faces[4*i+1]), get_support_point(faces[4*i+3])); 
      //   3          2    
      //    o--------o    
      //      .      |
      //        .    |   
      //          .  |  
      //            .| 
      //             o
      //              1
      triangles[1].set_support_points(get_support_point(faces[4*i+2]), get_support_point(faces[4*i+3]), get_support_point(faces[4*i+1])); 

      // second configuration
      //   3          2    
      //    o--------o    
      //    |      . |
      //    |    .   |   
      //    |  .     |  
      //    |.       | 
      //    o--------o
      //   0          1
      //
      //              2    
      //             o    
      //           . |
      //         .   |   
      //       .     |  
      //     .       | 
      //    o--------o
      //   0          1
      triangles[2].set_support_points(get_support_point(faces[4*i+1]), get_support_point(faces[4*i+2]), get_support_point(faces[4*i+0])); 
      //   3          2    
      //    o--------o    
      //    |      .  
      //    |    .       
      //    |  .        
      //    |.         
      //    o         
      //   0           
      triangles[3].set_support_points(get_support_point(faces[4*i+3]), get_support_point(faces[4*i+0]), get_support_point(faces[4*i+2])); 

      if (((!triangles[0].above_point(p)) || (!triangles[1].above_point(p))) && ((!triangles[2].above_point(p)) || (!triangles[3].above_point(p)))) { return false; }
    } // end for i
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
    std::vector<double> tmp = point_candidate;
    point_candidate = std::vector<double>(3, 0.0);
    std::vector<double> global_coordinates = compute_residual_vector(local_coordinates);
    for (unsigned int i = 0; i < 3; ++i) { global_coordinates[i] *= -1.0; }
    point_candidate = tmp; 
    return global_coordinates;
  }

  bool contains_point(const std::vector<double> &coordinates, bool coordinates_are_local = false) {
    assert(coordinates.size() == 3); 
    std::vector<double> local_coordinates;
    if (!coordinates_are_local) {
      point_candidate = coordinates;
      if (!within_bounding_box(point_candidate)) { return false; }
      if (!contained_by_triangles_on_faces(point_candidate)) { return false; }
      local_coordinates = solve_newton();
    } else {
      local_coordinates = coordinates;
    } // end if
    for (unsigned int i = 0; i < 3; ++i) {
      if (fabs(local_coordinates[i]) > 1.0 + epsilon) {
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
    if (!contained_by_triangles_on_faces(p)) { return std::pair<unsigned int, std::vector<double> >(99, std::vector<double>(2, 0.0)); }
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
    if (fabs(determinant) < 1.0e-16) {
      std::cerr<<"determinant="<<determinant<<"\n";
      std::cerr<<"matrix=\n";
      for (unsigned int i = 0; i < 9; ++i) {
        std::cerr<<A[i]<<((i%3 == 2) ? "\n" : "  ");  
      } // end for i
      for (unsigned int i = 0; i < 1000; ++i) { std::cerr<<std::flush; }
    } // end if
    assert(fabs(determinant) > 1.0e-16);
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

  std::vector<double> compute_initial_guess() {
    // matlab code
    // M=[-1 1 1 -1 -1 1 1 -1; -1 -1 1 1 -1 -1 1 1; -1 -1 -1 -1 1 1 1 1; 1 1 1 1 1 1 1 1]
    // K=transpose(M)*inv(M*transpose(M))
    // Y=[...]
    // KA=K(:,1:3)
    // Kb=K(:,4)
    // A=Y*KA
    // b=Y*Kb
    // p=[...]
    // inv(A)*(p-b)
    double KA[24] = {
      -0.125, -0.125, -0.125, 
       0.125, -0.125, -0.125, 
       0.125,  0.125, -0.125, 
      -0.125,  0.125, -0.125, 
      -0.125, -0.125,  0.125, 
       0.125, -0.125,  0.125, 
       0.125,  0.125,  0.125, 
      -0.125,  0.125,  0.125, 
    };

    double Kb[8] = {
      0.125,
      0.125,
      0.125,
      0.125,
      0.125,
      0.125,
      0.125,
      0.125,
    };

    double Y[24];
    for (unsigned int i = 0; i < 8; ++i) {
      for (unsigned int j = 0; j < 3; ++j) {
        Y[j*8+i] = support_points[i*3+j];
      } // end for j
    } // end for i


    std::vector<double> A(9, 0.0);
    for (unsigned int i = 0; i < 3; ++i) {
      for (unsigned int j = 0; j < 3; ++j) {
        A[3*i+j] = 0.0;
        for (unsigned int k = 0; k < 8; ++k) {
          A[3*i+j] += Y[8*i+k] * KA[3*k+j];
        } // end for k
      } // end for j
    } // end for i

    std::vector<double> b(3, 0.0);
    for (unsigned int i = 0; i < 3; ++i) {
      b[i] = 0.0;
      for (unsigned int j = 0; j < 8; ++j) {
        b[i] += Y[8*i+j] * Kb[j];
      } // end for j
      b[i] = point_candidate[i] - b[i];
    } // end for i

    return compute_matrix_times_vector(compute_inverse_matrix(A), b);
  }

  // map the coordinates of the point candidate onto the reference frame of the volume element defined by the support points
  std::vector<double> solve_newton(double abs_tol = 1.0e-14, double rel_tol = 1.0e-14, unsigned int max_iter = 100, bool verbose = false) {
    if (verbose) { std::cout<<"solve newton with line search\n"; }
    std::vector<double> x(3, 0.0);
    x = compute_initial_guess();
//    std::cout<<"initial guess=\n";
//    for (unsigned int i = 0; i < 3; ++i) { std::cout<<x[i]<<"\n"; }
    std::vector<double> residual_vector = compute_residual_vector(x);
    double residual_norm = sqrt(std::inner_product(residual_vector.begin(), residual_vector.end(), residual_vector.begin(), 0.0));
    double tol = abs_tol + rel_tol * residual_norm; 
    for (unsigned int iter = 0; iter < max_iter; ++iter) {
      if (verbose) { std::cout<<iter<<"  "<<residual_norm<<std::endl; }
      if (residual_norm < tol) { 
        if (verbose) { std::cout<<"converged at iteration "<<iter<<" with residual norm "<<residual_norm<<"\n"; }
        return x; 
      } // end if
      std::vector<double> inverse_jacobian_matrix_times_residual_vector = compute_inverse_jacobian_times_residual(compute_jacobian_matrix(x), residual_vector);
      bool line_search_passed = false;
      for (double alpha = 1.0; alpha > 1.0e-14; alpha /= 2.0) {
        std::vector<double> tmp(3);
        for (unsigned int i = 0; i < 3; ++i) { tmp[i] = x[i] - alpha * inverse_jacobian_matrix_times_residual_vector[i]; };
        residual_vector = compute_residual_vector(tmp);
        double tmp_residual_norm = sqrt(std::inner_product(residual_vector.begin(), residual_vector.end(), residual_vector.begin(), 0.0));
        if (tmp_residual_norm < residual_norm) {
          x = tmp;
          residual_norm = tmp_residual_norm;
          line_search_passed = true;
          break;
        } // end if
      } // end for
      assert(line_search_passed);
    } // end for
    std::cerr<<"failed to converge with tolerance "<<tol<<" after "<<max_iter-1<<" iterations (residual norm was "<<residual_norm<<")"<<std::endl;
    std::cerr<<"support_points=\n";
    for (unsigned int i = 0; i < 8; ++i) {
      std::cerr<<i<<"  ["<<support_points[3*i]<<", "<<support_points[3*i+1]<<", "<<support_points[3*i+2]<<"]\n";  
    } // end for i
    std::cerr<<"bounding_box=\n";
    for (unsigned int i = 0; i < 3; ++i) { std::cerr<<bounding_box[i]<<"  "<<bounding_box[i+3]<<"\n"; }
    std::cerr<<"point_candidate=["<<point_candidate[0]<<", "<<point_candidate[1]<<", "<<point_candidate[2]<<"]\n";
    for (unsigned int i = 0; i < 1000; ++i) { std::cerr<<std::flush; }
    abort(); 
  }
};

void test_inverted_element(const std::vector<double> &support_points) {
  assert(support_points.size() == 24);
  unsigned int faces[24] = {
    0, 3, 2, 1, 
    0, 1, 5, 4, 
    1, 2, 6, 5, 
    2, 3, 7, 6, 
    3, 0, 4, 7, 
    4, 5, 6, 7
  };
  double scalar_product = 0.0;
  for (unsigned int i = 0; i < 6; ++i) {
    scalar_product = 0.0;
    for (unsigned int j = 0; j < 3; ++j) {
      scalar_product += (support_points[3*faces[4*i+0]+j] - support_points[3*faces[4*i+1]+j]) *
        (support_points[3*faces[4*i+3]+j] - support_points[3*faces[4*i+2]+j]);
    } // end for j
    assert(scalar_product > 0.0);
    scalar_product = 0.0;
    for (unsigned int j = 0; j < 3; ++j) {
      scalar_product += (support_points[3*faces[4*i+1]+j] - support_points[3*faces[4*i+2]+j]) *
        (support_points[3*faces[4*i+0]+j] - support_points[3*faces[4*i+3]+j]);
    } // end for j
    assert(scalar_product > 0.0);
  } // end for i
}

bool quick_contains_point(const hex8_element_t &ve, const std::vector<double> &point_candidate) {
  unsigned int faces[24] = {
    0, 3, 2, 1, 
    0, 1, 5, 4, 
    1, 2, 6, 5, 
    2, 3, 7, 6, 
    3, 0, 4, 7, 
    4, 5, 6, 7
  };
  triangle_t triangle;
  for (unsigned int i = 0; i < 6; ++i) {
    bool first_configuration_passed = true;
    // first configuration when splitting face into two triangles
    //   3          2    
    //    o--------o    
    //    | .      |
    //    |   .    |   
    //    |     .  |  
    //    |       .| 
    //    o--------o
    //   0          1
    //
    //   3        
    //    o
    //    | .
    //    |   .
    //    |     .
    //    |       .
    //    o--------o
    //   0          1
    triangle.set_support_points(ve.get_support_point(faces[4*i+0]), ve.get_support_point(faces[4*i+1]), ve.get_support_point(faces[4*i+3])); 
    if (!triangle.above_point(point_candidate)) { first_configuration_passed = false; }
    //   3          2    
    //    o--------o    
    //      .      |
    //        .    |   
    //          .  |  
    //            .| 
    //             o
    //              1
    triangle.set_support_points(ve.get_support_point(faces[4*i+2]), ve.get_support_point(faces[4*i+3]), ve.get_support_point(faces[4*i+1])); 
    if (!triangle.above_point(point_candidate)) { first_configuration_passed = false; }


    bool second_configuration_passed = true;
    // second configuration
    //   3          2    
    //    o--------o    
    //    |      . |
    //    |    .   |   
    //    |  .     |  
    //    |.       | 
    //    o--------o
    //   0          1
    //
    //              2    
    //             o    
    //           . |
    //         .   |   
    //       .     |  
    //     .       | 
    //    o--------o
    //   0          1
    triangle.set_support_points(ve.get_support_point(faces[4*i+1]), ve.get_support_point(faces[4*i+2]), ve.get_support_point(faces[4*i+0])); 
    if (!triangle.above_point(point_candidate)) { second_configuration_passed = false; }
    //   3          2    
    //    o--------o    
    //    |      .  
    //    |    .       
    //    |  .        
    //    |.         
    //    o         
    //   0           
    triangle.set_support_points(ve.get_support_point(faces[4*i+3]), ve.get_support_point(faces[4*i+0]), ve.get_support_point(faces[4*i+2])); 
    if (!triangle.above_point(point_candidate)) { second_configuration_passed = false; }


    if ((!first_configuration_passed) && (!second_configuration_passed)) { return false; }
  } // end for i
  return true;
}

#endif // HEX8_ELEMENT_T
