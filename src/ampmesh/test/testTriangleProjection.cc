#include <utils/AMPManager.h>
#include <utils/UnitTest.h>

#include <triangle_t.h>
#include <euclidean_geometry_tools.h>

void test_above_point(triangle_t * t_ptr, unsigned int n_random_candidate_points = 10000) {
  double const * centroid = t_ptr->get_centroid();
  double const * normal = t_ptr->get_normal();
  double tmp[9];
  make_vector_from_two_points(t_ptr->get_support_point_ptr(0), t_ptr->get_support_point_ptr(1), tmp+0);
  make_vector_from_two_points(t_ptr->get_support_point_ptr(1), t_ptr->get_support_point_ptr(2), tmp+3);
  make_vector_from_two_points(t_ptr->get_support_point_ptr(2), t_ptr->get_support_point_ptr(0), tmp+6);
  double random_candidate_point[3], random_motion[4]; 
  bool triangle_above_random_candidate_point;
  double tolerance = 1.0e-12;
  for (unsigned int i = 0; i < n_random_candidate_points; ++i) {
    for (unsigned int j = 0; j < 4; ++j) { random_motion[j] = -1.0+2.0*rand()/RAND_MAX; }
    triangle_above_random_candidate_point = (random_motion[0] < -tolerance);
    for (unsigned int j = 0; j < 3; ++j) { random_candidate_point[j] = centroid[j] + random_motion[0]*normal[j] + random_motion[1]*tmp[0+j] + random_motion[2]*tmp[3+j] + random_motion[3]*tmp[6+j]; }
    AMP_ASSERT(t_ptr->above_point(random_candidate_point, tolerance) == triangle_above_random_candidate_point);
  } // end for i
}

void test_contains_point(triangle_t * t_ptr, unsigned int n_random_candidate_points = 10000) {
  double const * triangle_centroid = t_ptr->get_centroid();
  double const * triangle_normal = t_ptr->get_normal();
  for (unsigned int i = 0; i < 3; ++i) {
    edge_t * e_ptr = t_ptr->get_edge(i);
    double tmp[3];
    make_vector_from_two_points(e_ptr->get_support_point_ptr(0), e_ptr->get_support_point_ptr(1), tmp);
    double const * edge_normal = e_ptr->get_normal();
    double const * edge_center = e_ptr->get_center();
    AMP_ASSERT(fabs(compute_scalar_product(tmp, edge_normal)) < 1.0e-14);
    AMP_ASSERT(fabs(compute_scalar_product(triangle_normal, edge_normal)) < 1.0e-14);
    make_vector_from_two_points(triangle_centroid, edge_center, tmp);
    double random_candidate_point[3], random_motion[2];
    bool triangle_contains_random_candidate_point;
    double tolerance = 1.0e-12;
    for (unsigned int j = 0; j < n_random_candidate_points; ++j) {
      for (unsigned int k = 0; k < 2; ++k) { random_motion[k] = -1.0+2.0*rand()/RAND_MAX; }
      for (unsigned int k = 0; k < 3; ++k) { random_candidate_point[k] = edge_center[k] + random_motion[0] * tmp[k] + random_motion[1] * triangle_normal[k]; }
      triangle_contains_random_candidate_point = (random_motion[0] < tolerance);
      AMP_ASSERT(e_ptr->above_point(random_candidate_point, tolerance) == triangle_contains_random_candidate_point);
      AMP_ASSERT(t_ptr->contains_point(random_candidate_point, tolerance) == triangle_contains_random_candidate_point);
    } // end for j
  } // end for i

  for (unsigned int i = 0; i < 3; ++i) {
    double tmp[3];
    double const * triangle_support_point = t_ptr->get_support_point_ptr(i);
    make_vector_from_two_points(triangle_centroid, triangle_support_point, tmp);
    double random_candidate_point[3], random_motion[2];
    bool triangle_contains_random_candidate_point;
    double tolerance = 1.0e-12;
    for (unsigned int j = 0; j < n_random_candidate_points; ++j) {
      for (unsigned int k = 0; k < 2; ++k) { random_motion[k] = -1.0+2.0*rand()/RAND_MAX; }
      for (unsigned int k = 0; k < 3; ++k) { random_candidate_point[k] = triangle_support_point[k] + random_motion[0] * tmp[k] + random_motion[1] * triangle_normal[k]; }
      triangle_contains_random_candidate_point = (random_motion[0] < tolerance);
      AMP_ASSERT(t_ptr->contains_point(random_candidate_point, tolerance) == triangle_contains_random_candidate_point);
    } // end for j
  } // end for i
}

void myTest(AMP::UnitTest *ut, std::string exeName) {

  double points[9] = {
    0.0, 0.0, 0.0, // 0
    1.0, 0.0, 0.0, // 1
    0.0, 1.0, 0.0, // 2
  }; 


  double scaling_factors[3] = { 4.0, 2.0, 1.0 };
  scale_points(std::vector<double>(scaling_factors, scaling_factors+3), 3, points);

  double translation_vector[3] = { 3.0, 1.0, 5.0 };
  translate_points(std::vector<double>(translation_vector, translation_vector+3), 3, points);

  rotate_points(2, M_PI/3.0, 3, points);

  rotate_points(0, 0.75*M_PI, 3, points);

  srand(0);
  for (unsigned int i = 0; i < 9; ++i) { points[i] += -0.1 + 0.2*rand()/RAND_MAX; }

  triangle_t triangle(points, points+3, points+6);

  test_above_point(&triangle);

  test_contains_point(&triangle);

  ut->passes(exeName);
}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::string exeName = "testTriangleProjection";

  try {
    myTest(&ut, exeName);
  } catch (std::exception &err) {
    std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
    ut.failure("ERROR: While testing");
  } catch( ... ) {
    std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
    ut.failure("ERROR: While testing");
  }

  ut.report();
  int num_failed = ut.NumFailGlobal();

  AMP::AMPManager::shutdown();
  return num_failed;
}  
