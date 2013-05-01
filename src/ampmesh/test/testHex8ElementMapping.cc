#include <utils/AMPManager.h>
#include <utils/UnitTest.h>

#include <ampmesh/hex8_element_t.h>
#include <ampmesh/latex_visualization_tools.h>
#include <ampmesh/euclidean_geometry_tools.h>

class soft_equal_to {
public:
  soft_equal_to(double const tol = 1.0e-15) : tolerance(tol) { };
  inline bool operator()(double const x, double const y) { return std::abs(x-y) < tolerance; };
private:
  double tolerance;
};

//inline bool soft_equal_to(double const x, double const y, double const tolerance = 1.0e-15) { return std::abs(x-y) < tolerance; }

unsigned int test_mapping_global_to_local(hex8_element_t *volume_element, unsigned int n_random_candidate_points = 10000, double abs_tol = 1.0e-12, double rel_tol = 1.0e-12) {
  double local_coordinates[3], global_coordinates[3], computed_local_coordinates[3], tol, error_vector[3];
  unsigned int count_tests_failing = 0;
  for (unsigned int i = 0; i < n_random_candidate_points; ++i) {
    for (unsigned int j = 0; j < 3; ++j) { local_coordinates[j] = -1.0+2.0*rand()/RAND_MAX; }
    volume_element->map_local_to_global(local_coordinates, global_coordinates);
    volume_element->map_global_to_local(global_coordinates, computed_local_coordinates);
    make_vector_from_two_points(computed_local_coordinates, local_coordinates, error_vector);
    tol = abs_tol + rel_tol * compute_vector_norm(local_coordinates);
    if (compute_vector_norm(error_vector) > tol) { ++count_tests_failing; }
  } // end for i
  return count_tests_failing;
}

void test_mapping_face_to_local(hex8_element_t *volume_element, unsigned int n_random_candidate_points = 10000, double tol = 1.0e-16) {
  double local_coordinates_on_face[3], local_coordinates[3], computed_local_coordinates_on_face[3], computed_shift[3], error_vector[3];
  local_coordinates_on_face[2] = 0.0;
  computed_local_coordinates_on_face[2] = 0.0;
  for (unsigned int i = 0; i < n_random_candidate_points; ++i) {
    for (unsigned int j = 0; j < 2; ++j) { local_coordinates_on_face[j] = -1.0+2.0*rand()/RAND_MAX; }
    for (unsigned int f = 0; f < 6; ++f) {
      volume_element->map_face_to_local(f, local_coordinates_on_face, local_coordinates);
      volume_element->project_on_face(f, local_coordinates, computed_local_coordinates_on_face, computed_shift);
      AMP_ASSERT(compute_vector_norm(computed_shift) < tol);
      make_vector_from_two_points(computed_local_coordinates_on_face, local_coordinates_on_face, error_vector);
      AMP_ASSERT(compute_vector_norm(error_vector) < tol); 
    } // end for f
  } // end for i
}

void test_basis_functions_values_on_face(hex8_element_t *volume_element, unsigned int n_random_candidate_points = 10000, double tol = 1.0e-16) {
  double local_coordinates_on_face[2], local_coordinates[3], basis_functions_values_on_face[4], basis_functions_values[8]; 
  unsigned int const * face_ordering; 
  for (unsigned int i = 0; i < n_random_candidate_points; ++i) {
    for (unsigned int j = 0; j < 2; ++j) { local_coordinates_on_face[j] = -1.0+2.0*rand()/RAND_MAX; }
    for (unsigned int f = 0; f < 6; ++f) {
      face_ordering = hex8_element_t::get_face(f);
      volume_element->map_face_to_local(f, local_coordinates_on_face, local_coordinates);
      hex8_element_t::get_basis_functions_values(local_coordinates, basis_functions_values);
      hex8_element_t::get_basis_functions_values_on_face(local_coordinates_on_face, basis_functions_values_on_face);
      for (unsigned int v = 0; v < 4; ++v) {
        AMP_ASSERT(std::abs(basis_functions_values_on_face[v] - basis_functions_values[face_ordering[v]]) < tol);
      } // end for v
    } // end for f
  } // end for i
}

void test_mapping_basis_functions_values_to_local_coordinates_on_face(unsigned int n_random_candidate_points = 1000, double tol = 1.0e-16) {
  double local_coordinates_on_face[2], basis_functions_values_on_face[4], computed_local_coordinates_on_face[2];
  for (unsigned int i = 0; i < n_random_candidate_points; ++i) {
    for (unsigned int j = 0; j < 2; ++j) { local_coordinates_on_face[j] = -1.0+2.0*rand()/RAND_MAX; }
    hex8_element_t::get_basis_functions_values_on_face(local_coordinates_on_face, basis_functions_values_on_face);
    hex8_element_t::get_local_coordinates_on_face(basis_functions_values_on_face, computed_local_coordinates_on_face);
    assert(std::equal(local_coordinates_on_face, local_coordinates_on_face+2, computed_local_coordinates_on_face, soft_equal_to()));
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
  volume_element.set_support_points(points);
  AMP_ASSERT(test_mapping_global_to_local(&volume_element) == 0);

  test_mapping_face_to_local(&volume_element);

  test_basis_functions_values_on_face(&volume_element);

  test_mapping_basis_functions_values_to_local_coordinates_on_face();

///*
  double scaling_factors[3] = { 4.0, 2.0, 1.0 };
  scale_points(scaling_factors, 8, points);
  volume_element.set_support_points(points);
  AMP_ASSERT(test_mapping_global_to_local(&volume_element) == 0);
//*/

  double translation_vector[3] = { 3.0, 1.0, 5.0 };
  translate_points(translation_vector, 8, points);
  volume_element.set_support_points(points);
  AMP_ASSERT(test_mapping_global_to_local(&volume_element) == 0);

  rotate_points(2, pi/3.0, 8, points);
  volume_element.set_support_points(points);
  AMP_ASSERT(test_mapping_global_to_local(&volume_element) == 0);

  rotate_points(0, 0.75*pi, 8, points);
  volume_element.set_support_points(points);
  AMP_ASSERT(test_mapping_global_to_local(&volume_element) == 0);

  srand(0);
  for (unsigned int i = 0; i < 24; ++i) { points[i] += -0.1 + 0.2*rand()/RAND_MAX; }

  volume_element.set_support_points(points);

//  test_mapping(&volume_element, 1000000);
  AMP_ASSERT(test_mapping_global_to_local(&volume_element) == 0);
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
