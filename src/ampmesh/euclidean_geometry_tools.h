#ifndef EUCLIDEAN_GEOMETRY_TOOLS
#define EUCLIDEAN_GEOMETRY_TOOLS

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
double compute_vector_norm(double const * vector);
void normalize_vector(double * vector);

#endif // EUCLIDEAN_GEOMETRY_TOOLS
