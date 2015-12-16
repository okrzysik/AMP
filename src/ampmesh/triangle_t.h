
#ifndef TRIANGLE_T
#define TRIANGLE_T

#include <ampmesh/edge_t.h>

#include <vector>

class triangle_t {
public:
    triangle_t( double const *A, double const *B, double const *C );
    ~triangle_t();
    void set_support_points( double const *A, double const *B, double const *C );
    double const *get_support_point_ptr( unsigned int i ) const;
    double const *get_normal();
    double const *get_centroid();
    edge_t *get_edge( unsigned int i );
    bool above_point( double const *point, double tolerance = 1.0e-12 );
    bool contains_point( double const *point, double tolerance = 1.0e-12 );
    int project_point( double const *point, double *projection, double tolerance = 1.0e-12 );

private:
    edge_t *edges_ptr[3];
    double const *support_points_ptr[3];
    std::vector<double> normal;
    std::vector<double> centroid;
    std::vector<double> tmp;
    bool normal_updated;
    bool centroid_updated;
    bool edges_updated;

    void build_edges();
    void clear_edges();
    void compute_normal();
    void compute_centroid();
    double compute_distance_to_containing_plane( double const *point );
};

int project_point_onto_collection_of_triangles( unsigned int n_triangles,
                                                triangle_t **triangle_ptr,
                                                double const *point,
                                                double *projection,
                                                unsigned int &position,
                                                double &distance,
                                                double tolerance = 1.0e-12 );

#endif // TRIANGLE_T
