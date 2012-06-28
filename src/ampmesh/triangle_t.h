#ifndef TRIANGLE_T
#define TRIANGLE_T

#include <ampmesh/edge_t.h>

#include <vector>

class triangle_t {
public:
  triangle_t(double const * A, double const * B, double const * C);
  void set_support_points(double const * A, double const * B, double const * C);
  void set_support_points(double const * * ptr);
  bool above_point(double const * p, double tolerance = 1.0e-12);
  double const * get_support_point_ptr(unsigned int i) const;
  double const * get_normal();
  double const * get_centroid();
  edge_t * get_edge(unsigned int i);
  bool contains_point(double const * p, double tolerance = 1.0e-12);

private:
  std::vector<double const *> support_points_ptr;
  std::vector<double> normal, centroid, tmp;
  bool normal_updated, centroid_updated;
  void compute_normal();
  void compute_centroid();

  std::vector<edge_t> edges;
  bool edges_updated;
  void build_edges();

};

#endif // TRIANGLE_T
