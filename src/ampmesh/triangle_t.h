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
void compute_cross_product(double const * const u, double const * const v, double * const w);
double compute_scalar_product(double const * const u, double const * const v);
void make_vector_from_two_points(double const * const start_point, double const * const end_point, double * const vector);

class triangle_t {
public:
  triangle_t(double const * A, double const * B, double const * C);
  void set_support_points(double const * A, double const * B, double const * C);
  bool above_point(double const * p, double tolerance = 1.0e-12);

private:
  double const * _A;
  double const * _B;
  double const * _C;
  bool memory_allocated;
  std::vector<double> _AB, _AC, _n_ABC;
  bool _n_ABC_updated;
  std::vector<double> _Ap, _Bp, _Cp, _ABCp;
  void compute_n_ABC();
};

#endif // TRIANGLE_T
