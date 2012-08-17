#ifndef EUCLIDEAN_GEOMETRY_TOOLS
#define EUCLIDEAN_GEOMETRY_TOOLS

void scale_points(unsigned int direction, double scaling_factor, unsigned int n_points, double * points);
void scale_points(double const * scaling_factors, unsigned int n_points, double * points);
void translate_points(unsigned int direction, double distance, unsigned int n_points, double * points);
void translate_points(double const * translation_vector, unsigned int n_points, double * points);
void rotate_points(unsigned int rotation_axis, double rotation_angle, unsigned int n_points, double * points);

void compute_cross_product(double const * u, double const * v, double * w);
double compute_scalar_product(double const * u, double const * const v);
void make_vector_from_two_points(double const * start_point, double const * end_point, double * vector);
double compute_vector_norm(double const * vector);
void normalize_vector(double * vector);
double compute_distance_between_two_points(double const * start_point, double const * end_point);

#endif // EUCLIDEAN_GEOMETRY_TOOLS
