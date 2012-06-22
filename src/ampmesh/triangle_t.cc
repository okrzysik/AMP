#include <ampmesh/triangle_t.h> 

triangle_t::triangle_t(double const * A, double const * B, double const * C) {
  memory_allocated = false;
  set_support_points(A, B, C);
}

void triangle_t::set_support_points(double const * A, double const * B, double const * C) {
  _A = A;
  _B = B;
  _C = C;
  _n_ABC_updated = false;
}

// check whether plane supported by the triangle ABC is above the point with respect to the normal
bool triangle_t::above_point(double const * p, double tolerance) {
  if (!memory_allocated) {
    _AB.resize(3);
    _AC.resize(3);
    _n_ABC.resize(3);
    _Ap.resize(3);
    _Bp.resize(3);
    _Cp.resize(3);
    _ABCp.resize(3);
  } // end if
  if (!_n_ABC_updated) { 
    compute_n_ABC();
    _n_ABC_updated = true;
  } // end if
  make_vector_from_two_points(_A, p, &(_Ap[0]));
  make_vector_from_two_points(_B, p, &(_Bp[0]));
  make_vector_from_two_points(_C, p, &(_Cp[0]));
  for (unsigned int i = 0; i < 3; ++i) { _ABCp[i] = (_Ap[i] + _Bp[i] + _Cp[i]) / 3.0; }
  return (compute_scalar_product(_ABCp, _n_ABC) < tolerance);
}

void triangle_t::compute_n_ABC() {
  make_vector_from_two_points(_A, _B, &(_AB[0]));
  make_vector_from_two_points(_A, _C, &(_AC[0]));
  compute_cross_product(&(_AB[0]), &(_AC[0]), &(_n_ABC[0]));
  double normalizing_factor = 1.0 / sqrt(std::inner_product(&(_n_ABC[0]), &(_n_ABC[0])+3, &(_n_ABC[0]), 0.0));
  assert(normalizing_factor < 1.0e12);
  for (unsigned int i = 0; i < 3; ++i) { scale_points(i, normalizing_factor, 1, &(_n_ABC[0])); }
}


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
  compute_cross_product(&(u[0]), &(v[0]), &(w[0]));
  return w;
}

void compute_cross_product(double const * u, double const * v, double * w) {
  w[0] = u[1]*v[2]-u[2]*v[1];
  w[1] = u[2]*v[0]-u[0]*v[2];
  w[2] = u[0]*v[1]-u[1]*v[0];
}

double compute_scalar_product(double const * u, double const * v) {
  return std::inner_product(&(u[0]), &(u[0])+3, &(v[0]), 0.0);
}

double compute_scalar_product(const std::vector<double> &u, const std::vector<double> &v) {
  assert(u.size() == 3);
  assert(v.size() == 3);
  return compute_scalar_product(&(u[0]), &(v[0]));
}

std::vector<double> make_vector_from_two_points(const std::vector<double> &start_point, const std::vector<double> &end_point) {
  assert(start_point.size() == 3);
  assert(end_point.size() == 3);
  std::vector<double> vector(3);
  make_vector_from_two_points(&(start_point[0]), &(end_point[0]), &(vector[0]));
  return vector;
}

void make_vector_from_two_points(double const * start_point, double const * end_point, double * vector) {
  for (unsigned int i = 0; i < 3; ++i) { vector[i] = end_point[i] - start_point[i]; }
}
