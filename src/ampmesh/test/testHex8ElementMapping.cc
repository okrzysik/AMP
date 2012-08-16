#include <utils/AMPManager.h>
#include <utils/UnitTest.h>

#include <ampmesh/hex8_element_t.h>
#include <ampmesh/latex_visualization_tools.h>
#include <ampmesh/euclidean_geometry_tools.h>

void test_mapping(hex8_element_t *volume_element, unsigned int n_random_candidate_points = 10000, double abs_tol = 1.0e-12, double rel_tol = 1.0e-12) {
  std::vector<double> random_candidate_point(3), candidate_point_global_coordinates(3), candidate_point_local_coordinates(3), mapping_error(3);
  double mapping_error_norm, tol;
  for (unsigned int i = 0; i < n_random_candidate_points; ++i) {
    for (unsigned int j = 0; j < 3; ++j) { random_candidate_point[j] = -1.0+2.0*rand()/RAND_MAX; }
    volume_element->map_local_to_global(&(random_candidate_point[0]), &(candidate_point_global_coordinates[0]));
    volume_element->map_global_to_local(&(candidate_point_global_coordinates[0]), &(candidate_point_local_coordinates[0]));

    for (unsigned int i = 0; i < 3; ++i) { mapping_error[i] = candidate_point_local_coordinates[i] - random_candidate_point[i]; }
    mapping_error_norm = compute_vector_norm(&(mapping_error[0]));
    tol = abs_tol + rel_tol * compute_vector_norm(&(random_candidate_point[0]));
    AMP_ASSERT(mapping_error_norm < tol);
  } // end for i
}

void myTest(AMP::UnitTest *ut, std::string exeName) {
  const double pi = 3.141592653589793;
  double points[24] = {
    -1.0, -1.0, -1.0, // 0
    +1.0, -1.0, -1.0, // 1
    +1.0, +1.0, -1.0, // 2
    -1.0, +1.0, -1.0, // 3
    -1.0, -1.0, +1.0, // 4
    +1.0, -1.0, +1.0, // 5
    +1.0, +1.0, +1.0, // 6
    -1.0, +1.0, +1.0  // 7
  }; 
  // shifting the points from [-1, 1]^3 to [0, 1]^3
  for (unsigned int i = 0; i < 24; ++i) { points[i] = 0.5*(points[i]+1.0); }

  hex8_element_t volume_element(points);

/*
  double scaling_factors[3] = { 4.0, 2.0, 1.0 };
  scale_points(std::vector<double>(scaling_factors, scaling_factors+3), 8, points);
*/

  double translation_vector[3] = { 3.0, 1.0, 5.0 };
  translate_points(std::vector<double>(translation_vector, translation_vector+3), 8, points);

  rotate_points(2, pi/3.0, 8, points);

  rotate_points(0, 0.75*pi, 8, points);

  srand(0);
  for (unsigned int i = 0; i < 24; ++i) { points[i] += -0.1 + 0.2*rand()/RAND_MAX; }

  volume_element.set_support_points(points);

//  test_mapping(&volume_element, 1000000);
  test_mapping(&volume_element);
//  std::cout<<"[test mapping] newton count = "<<volume_element.newton_count<<"\n";

  translate_points(0, 2.0, 8, points);
  translate_points(1, 9.0, 8, points);
  translate_points(2, 2.0, 8, points);
  volume_element.set_support_points(points);

  for (unsigned int i = 0; i < 8; ++i) { draw_point(volume_element.get_support_point(i), "red"); }
  double point_of_view[3] = { 1.0, 1.0, 1.0 };

  std::cout<<"% volume element\n";
  draw_hex8_element(&volume_element, point_of_view);

  std::cout<<"% bounding polyhedron\n";
  draw_bounding_polyhedron(&volume_element, point_of_view);

  std::cout<<"% bounding box\n";
  draw_bounding_box(&volume_element, point_of_view);

  ut->passes(exeName);
}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::string exeName = "testHex8ElementMapping";

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
