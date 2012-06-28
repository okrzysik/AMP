#ifndef EDGE_T 
#define EDGE_T

#include <vector>

class edge_t {
public:
  edge_t(double const * A, double const * B, double const * ABC);
  void set_support_points(double const * A, double const * B);
  void set_containing_plane(double const * ABC);
  double const * get_support_point_ptr(unsigned int i) const;
  double const * get_normal();
  double const * get_direction();
  double const * get_center();
  bool above_point(double const * p, double tolerance = 1.0e-12);
  bool project_point(double const * point_in_containing_plane, double * projectioni, double tolerance = 1.0e-12);

private:
  std::vector<double const *> support_points_ptr;
  double const * containing_plane_ptr;
  std::vector<double> normal, direction, center, tmp;
  bool normal_updated, direction_updated, center_updated;
  void compute_normal();
  void compute_direction();
  void compute_center();
  double compute_distance_to_containing_line(double const * point_in_containing_plane);

};

#endif // EDGE_T 
