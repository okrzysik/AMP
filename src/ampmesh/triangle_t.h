#ifndef TRIANGLE_T
#define TRIANGLE_T

#include <cassert>
#include <cmath>
#include <numeric>
#include <vector>

const double epsilon = 1.0e-12;

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

std::vector<double> compute_cross_product(const std::vector<double> &u, const std::vector<double> &v) {
  assert(u.size() == 3);
  assert(v.size() == 3);
  std::vector<double> w(3);
  w[0] = u[1]*v[2]-u[2]*v[1];
  w[1] = u[2]*v[0]-u[0]*v[2];
  w[2] = u[0]*v[1]-u[1]*v[0];
  return w;
}

double compute_scalar_product(const std::vector<double> &u, const std::vector<double> &v) {
  assert(u.size() == 3);
  assert(v.size() == 3);
  return std::inner_product(u.begin(), u.end(), v.begin(), 0.0);
}

std::vector<double> make_vector_from_two_points(const std::vector<double> &start_point, const std::vector<double> &end_point) {
  assert(start_point.size() == 3);
  assert(end_point.size() == 3);
  std::vector<double> vector(3);
  for (unsigned int i = 0; i < 3; ++i) { vector[i] = end_point[i] - start_point[i]; }
  return vector;
}

class triangle_t {
public:
  triangle_t() {
    double A[3] = { 0.0, 0.0, 0.0 };
    double B[3] = { 1.0, 0.0, 0.0 };
    double C[3] = { 0.0, 1.0, 0.0 };
    set_support_points(std::vector<double>(A, A+3), std::vector<double>(B, B+3), std::vector<double>(C, C+3)); 
  }
  triangle_t(const std::vector<double> &A, const std::vector<double> &B, const std::vector<double> &C) {
    set_support_points(A, B, C); 
  }
  void set_support_points(const std::vector<double> &A, const std::vector<double> &B, const std::vector<double> &C) { 
    assert(A.size() == 3);
    assert(B.size() == 3);
    assert(C.size() == 3);
    _A = A;
    _B = B;
    _C = C;
    compute_n_ABC();
  }

  // plane supported by the triangle ABC is above the point with respect to the normal
  bool above_point(const std::vector<double> &p) {
    assert(p.size() == 3);
    std::vector<double> A_p = make_vector_from_two_points(_A, p);
    std::vector<double> B_p = make_vector_from_two_points(_B, p);
    std::vector<double> C_p = make_vector_from_two_points(_C, p);
    std::vector<double> ABC_p(3);
    for (unsigned int i = 0; i < 3; ++i) { ABC_p[i] = A_p[i] + B_p[i] + C_p[i]; }
    return (compute_scalar_product(ABC_p, _n_ABC) < epsilon);
  }

private:
  std::vector<double> _A, _B, _C, _n_ABC;

  void compute_n_ABC() {
    _n_ABC = compute_cross_product(make_vector_from_two_points(_A, _B), make_vector_from_two_points(_A, _C));
    double normalizing_factor = sqrt(std::inner_product(_n_ABC.begin(), _n_ABC.end(), _n_ABC.begin(), 0.0));
    assert(normalizing_factor > 1.0e-12);
    scale_points(std::vector<double>(3, 1.0/normalizing_factor), 1, &(_n_ABC[0]));
  }
};

#endif // TRIANGLE_T
