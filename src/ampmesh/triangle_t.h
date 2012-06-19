#ifndef TRIANGLE_T
#define TRIANGLE_T

#include <cassert>
#include <cmath>
#include <numeric>
#include <vector>

const double epsilon = 1.0e-12;

void scale_points(unsigned int direction, double scaling_factor, unsigned int n_points, double* points);
void scale_points(const std::vector<double> &scaling_factors, unsigned int n_points, double* points);
void translate_points(unsigned int direction, double distance, unsigned int n_points, double* points);
void translate_points(std::vector<double> translation_vector, unsigned int n_points, double* points);
void rotate_points(unsigned int rotation_axis, double rotation_angle, unsigned int n_points, double* points);
std::vector<double> compute_cross_product(const std::vector<double> &u, const std::vector<double> &v);
double compute_scalar_product(const std::vector<double> &u, const std::vector<double> &v);
std::vector<double> make_vector_from_two_points(const std::vector<double> &start_point, const std::vector<double> &end_point);

class triangle_t {
public:
  triangle_t();
  triangle_t(const std::vector<double> &A, const std::vector<double> &B, const std::vector<double> &C);
  void set_support_points(const std::vector<double> &A, const std::vector<double> &B, const std::vector<double> &C);
  bool above_point(const std::vector<double> &p);

private:
  std::vector<double> _A, _B, _C, _n_ABC;
  void compute_n_ABC();
};

#endif // TRIANGLE_T
