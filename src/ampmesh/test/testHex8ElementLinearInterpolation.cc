#include "utils/AMPManager.h"
#include "utils/UnitTest.h"

#include "hex8_element_t.h"

double my_function(const std::vector<double> &xyz) {
  AMP_ASSERT(xyz.size() == 3);
  double x = xyz[0], y = xyz[1], z = xyz[2];
  return (1.0 + 6.0 * x) * (2.0 - 5.0 * y) * (3.0 + 4.0 * z);
}

double my_function_no_cross_terms(const std::vector<double> &xyz) {
  AMP_ASSERT(xyz.size() == 3);
  double x = xyz[0], y = xyz[1], z = xyz[2];
  return 1.0 + 6.0 * x - 5.0 * y + 4.0 * z;
}

unsigned int perform_battery_of_tests(hex8_element_t &volume_element, double (*my_function_to_call)(const std::vector<double> &)) {
  const unsigned int n_random_candidate_points = 100;
  const double tol_abs = 1.0e-13, tol_rel = 1.0e-13;
  std::vector<double> support_points = volume_element.get_support_points();
  std::vector<double> my_function_at_support_points(8);
  for (unsigned int i = 0; i < 8; ++i) { 
    my_function_at_support_points[i] = (*my_function_to_call)(std::vector<double>(&support_points[3*i], &support_points[3*i]+3));
  } // end for i

  unsigned int count_tests_failing = 0;
  for (unsigned int i = 0; i < n_random_candidate_points; ++i) {
    std::vector<double> candidate_point_local_coordinates(3);
    for (unsigned int j = 0; j < 3; ++j) { candidate_point_local_coordinates[j] = -1.0+2.0*rand()/RAND_MAX; }

    std::vector<double> basis_functions_values = get_basis_functions_values(candidate_point_local_coordinates);
    double interpolated_value = 0.0;
    for (unsigned int j = 0; j < 8; ++j) {
      interpolated_value += my_function_at_support_points[j] * basis_functions_values[j];
    } // end for j

    std::vector<double> candidate_point_global_coordinates = volume_element.map_local_to_global(candidate_point_local_coordinates);
    double my_function_at_candidate_point = (*my_function_to_call)(candidate_point_global_coordinates);
    double interpolation_error = fabs(interpolated_value-my_function_at_candidate_point);
    double tolerance = tol_abs+tol_rel*fabs(interpolated_value);
    if (interpolation_error > tolerance) { ++count_tests_failing; }
  } // end for i
  return count_tests_failing;
}

void myTest(AMP::UnitTest *ut, std::string exeName) {
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

  hex8_element_t volume_element(std::vector<double>(points, points+24));
  assert(perform_battery_of_tests(volume_element, my_function) == 0);
  assert(perform_battery_of_tests(volume_element, my_function_no_cross_terms) == 0);
  srand(0);

  double scaling_factors[3] = { 4.0, 2.0, 1.0 };
  scale_points(std::vector<double>(scaling_factors, scaling_factors+3), 8, points);
  volume_element.set_support_points(std::vector<double>(points, points+24));
  assert(perform_battery_of_tests(volume_element, my_function) == 0);
  assert(perform_battery_of_tests(volume_element, my_function_no_cross_terms) == 0);

  double translation_vector[3] = { 3.0, 1.0, 5.0 };
  translate_points(std::vector<double>(translation_vector, translation_vector+3), 8, points);
  volume_element.set_support_points(std::vector<double>(points, points+24));
  assert(perform_battery_of_tests(volume_element, my_function) == 0);
  assert(perform_battery_of_tests(volume_element, my_function_no_cross_terms) == 0);

  rotate_points(2, M_PI/2.0, 8, points);
  volume_element.set_support_points(std::vector<double>(points, points+24));
  assert(perform_battery_of_tests(volume_element, my_function) == 0);
  assert(perform_battery_of_tests(volume_element, my_function_no_cross_terms) == 0);

  rotate_points(0, -0.75*M_PI, 8, points);
  volume_element.set_support_points(std::vector<double>(points, points+24));
  assert(perform_battery_of_tests(volume_element, my_function) > 0);
  assert(perform_battery_of_tests(volume_element, my_function_no_cross_terms) == 0);

  rotate_points(0, 0.75*M_PI, 8, points);
  volume_element.set_support_points(std::vector<double>(points, points+24));
  assert(perform_battery_of_tests(volume_element, my_function) == 0);
  assert(perform_battery_of_tests(volume_element, my_function_no_cross_terms) == 0);

  for (unsigned int i = 0; i < 24; ++i) { points[i] += -0.1 + 0.2*rand()/RAND_MAX; }
  volume_element.set_support_points(std::vector<double>(points, points+24));
  assert(perform_battery_of_tests(volume_element, my_function) > 0);
  assert(perform_battery_of_tests(volume_element, my_function_no_cross_terms) == 0);

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



