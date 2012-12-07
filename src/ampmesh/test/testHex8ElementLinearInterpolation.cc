#include <utils/AMPManager.h>
#include <utils/UnitTest.h>

#include <ampmesh/hex8_element_t.h>
#include <ampmesh/euclidean_geometry_tools.h>

double my_function(double const *xyz) {
  double x = xyz[0], y = xyz[1], z = xyz[2];
  return (1.0 + 6.0 * x) * (2.0 - 5.0 * y) * (3.0 + 4.0 * z);
}

double my_function_no_cross_terms(double const *xyz) {
  double x = xyz[0], y = xyz[1], z = xyz[2];
  return 1.0 + 6.0 * x - 5.0 * y + 4.0 * z;
}

void test_normal(hex8_element_t *volume_element, unsigned int n_random_candidate_points = 20) {
  double normal_vector[3];
  double local_coordinates_on_face[2], local_coordinates[3], global_coordinates[3];
  double local_coordinates_on_face_check[2];
  double const * face_support_points_ptr[4];
  for (unsigned int i = 0; i < n_random_candidate_points; ++i) {
    for (unsigned int d = 0; d < 2; ++d) { local_coordinates_on_face[d] = -1.0+2.0*rand()/RAND_MAX; }
    for (unsigned int f = 0; f < 6; ++f) {
      for (unsigned int v = 0; v < 4; ++v) { face_support_points_ptr[v] = volume_element->get_support_point(volume_element->get_face(f)[v]); }
      volume_element->map_face_to_local(f, local_coordinates_on_face, local_coordinates);
      volume_element->map_local_to_face(f, local_coordinates, local_coordinates_on_face_check);
      AMP_ASSERT(std::equal(local_coordinates_on_face, local_coordinates_on_face+2, local_coordinates_on_face_check));
      volume_element->map_local_to_global(local_coordinates, global_coordinates);
      volume_element->compute_normal_to_face(f, local_coordinates, global_coordinates, normal_vector);
      std::cout<<f<<" { ";
      for (unsigned int d = 0; d < 3; ++d) { std::cout<<normal_vector[d]<<" "; }
      std::cout<<"}  ";
      volume_element->compute_normal_to_face(f, local_coordinates, normal_vector);
      std::cout<<f<<" { ";
      for (unsigned int d = 0; d < 3; ++d) { std::cout<<normal_vector[d]<<" "; }
      std::cout<<"}  ";
      volume_element->get_normal_to_face(face_support_points_ptr, local_coordinates_on_face, normal_vector);
      std::cout<<f<<" { ";
      for (unsigned int d = 0; d < 3; ++d) { std::cout<<normal_vector[d]<<" "; }
      std::cout<<"}\n";
    } // end for f
  } // end for i
}

inline bool soft_equal_to(double x, double y) { return abs(x-y) < 1.0e-15; }

void test_recovering_local_coordinates_on_face_from_basis_functions_values(unsigned int n_random_candidate_points = 1000) {
  double x[2], x_prime[2], phi[4];
  for (unsigned int i = 0; i < n_random_candidate_points; ++i) {
    for (unsigned j = 0; j < 2; ++j) {
      x[j] = -1.0+2.0*rand()/RAND_MAX;
    } // end for j
    hex8_element_t::get_basis_functions_values_on_face(x, phi);
    hex8_element_t::get_local_coordinates_on_face(phi, x_prime);
    assert(std::equal(x, x+2, x_prime, soft_equal_to));
  } // end for i
}

unsigned int perform_battery_of_tests(hex8_element_t *volume_element, double (*my_function_to_call)(double const *), unsigned int n_random_candidate_points = 1000, double tol_abs = 1.0e-12, double tol_rel = 1.0e-12) {
  std::vector<double> my_function_at_support_points(8);
  for (unsigned int i = 0; i < 8; ++i) { 
    my_function_at_support_points[i] = (*my_function_to_call)(volume_element->get_support_point(i));
  } // end for i

  std::vector<double> candidate_point_local_coordinates(3), candidate_point_global_coordinates(3);
  std::vector<double> basis_functions_values(8);
  double interpolated_value, my_function_at_candidate_point, interpolation_error, tol;
  unsigned int count_tests_failing = 0;
  for (unsigned int i = 0; i < n_random_candidate_points; ++i) {
    for (unsigned int j = 0; j < 3; ++j) { candidate_point_local_coordinates[j] = -1.0+2.0*rand()/RAND_MAX; }

    hex8_element_t::get_basis_functions_values(&(candidate_point_local_coordinates[0]), &(basis_functions_values[0]));
    interpolated_value = 0.0;
    for (unsigned int j = 0; j < 8; ++j) {
      interpolated_value += my_function_at_support_points[j] * basis_functions_values[j];
    } // end for j

    volume_element->map_local_to_global(&(candidate_point_local_coordinates[0]), &(candidate_point_global_coordinates[0]));
    my_function_at_candidate_point = (*my_function_to_call)(&(candidate_point_global_coordinates[0]));
    interpolation_error = fabs(interpolated_value-my_function_at_candidate_point);
    tol = tol_abs+tol_rel*fabs(interpolated_value);
    if (interpolation_error > tol) { ++count_tests_failing; }
  } // end for i
  return count_tests_failing;
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

  hex8_element_t volume_element(points);
  AMP_ASSERT(perform_battery_of_tests(&volume_element, my_function) == 0);
  AMP_ASSERT(perform_battery_of_tests(&volume_element, my_function_no_cross_terms) == 0);
  test_normal(&volume_element);
  srand(0);

  double scaling_factors[3] = { 4.0, 2.0, 1.0 };
  scale_points(scaling_factors, 8, points);
  volume_element.set_support_points(points);
  AMP_ASSERT(perform_battery_of_tests(&volume_element, my_function) == 0);
  AMP_ASSERT(perform_battery_of_tests(&volume_element, my_function_no_cross_terms) == 0);

  double translation_vector[3] = { 3.0, 1.0, 5.0 };
  translate_points(translation_vector, 8, points);
  volume_element.set_support_points(points);
  AMP_ASSERT(perform_battery_of_tests(&volume_element, my_function) == 0);
  AMP_ASSERT(perform_battery_of_tests(&volume_element, my_function_no_cross_terms) == 0);

  rotate_points(2, pi/2.0, 8, points);
  volume_element.set_support_points(points);
  AMP_ASSERT(perform_battery_of_tests(&volume_element, my_function) == 0);
  AMP_ASSERT(perform_battery_of_tests(&volume_element, my_function_no_cross_terms) == 0);

  rotate_points(0, -0.75*pi, 8, points);
  volume_element.set_support_points(points);
  AMP_ASSERT(perform_battery_of_tests(&volume_element, my_function) > 0);
  AMP_ASSERT(perform_battery_of_tests(&volume_element, my_function_no_cross_terms) == 0);

  rotate_points(0, 0.75*pi, 8, points);
  volume_element.set_support_points(points);
  AMP_ASSERT(perform_battery_of_tests(&volume_element, my_function) == 0);
  AMP_ASSERT(perform_battery_of_tests(&volume_element, my_function_no_cross_terms) == 0);

  for (unsigned int i = 0; i < 24; ++i) { points[i] += -0.1 + 0.2*rand()/RAND_MAX; }
  volume_element.set_support_points(points);
  AMP_ASSERT(perform_battery_of_tests(&volume_element, my_function) > 0);
  AMP_ASSERT(perform_battery_of_tests(&volume_element, my_function_no_cross_terms) == 0);

  test_recovering_local_coordinates_on_face_from_basis_functions_values();

  ut->passes(exeName);
}


int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::string exeName = "testHex8ElementLinearInterpolation";

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



