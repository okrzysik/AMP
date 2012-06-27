#include "utils/AMPManager.h"
#include "utils/UnitTest.h"

#include "hex8_element_t.h"

void old_test_project_on_face() {
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

//  std::cout<<"volume element support points=\n";
//  for (unsigned int i = 0; i < 24; ++i) { std::cout<<points[i]<<"  "; if (i%3 == 2) {std::cout<<"\n";} }

  hex8_element_t volume_element(std::vector<double>(points, points+24));

  unsigned int n_candidates = 64;
  double candidates[] = {
// candidates coincide with support points
    0.0, 0.0, 0.0, // 0 -> 0 [-1, -1]
    1.0, 0.0, 0.0, // 1 -> 0 [-1, 1]
    1.0, 1.0, 0.0, // 2 -> 0 [1, 1]
    0.0, 1.0, 0.0, // 3 -> 0 [1, -1]
    0.0, 0.0, 1.0, // 4 -> 1 [-1, 1]
    1.0, 0.0, 1.0, // 5 -> 1 [1, 1]
    1.0, 1.0, 1.0, // 6 -> 2 [1, 1]
    0.0, 1.0, 1.0, // 7 -> 3 [1, 1]
// outside the volume element
    1.5, 0.0, 0.0, // 8 -> 99 [0, 0]
    0.0, 1.5, 0.0, // 9 -> 99 [0, 0]
    0.0, 0.0, 1.5, // 10 -> 99 [0, 0]
    1.5, 1.5, 0.0, // 11 -> 99 [0, 0]
    1.5, 0.0, 1.5, // 12 -> 99 [0, 0]
    0.0, 1.5, 1.5, // 13 -> 99 [0, 0]
    1.5, 1.5, 1.5, // 14 -> 99 [0, 0]
    -0.5, 0.0, 0.0, // 15 -> 99 [0, 0]
    0.0, -0.5, 0.0, // 16 -> 99 [0, 0]
    0.0, 0.0, -0.5, // 17 -> 99 [0, 0]
    -0.5, -0.5, 0.0, // 18 -> 99 [0, 0]
    -0.5, 0.0, -0.5, // 19 -> 99 [0, 0]
    0.0, -0.5, -0.5, // 20 -> 99 [0, 0]
    -0.5, -0.5, -0.5, // 21 -> 99 [0, 0]
// at center on faces
    0.5, 0.5, 0.0, // 22 -> 0 [0, 0]
    0.5, 0.0, 0.5, // 23 -> 1 [0, 0]
    1.0, 0.5, 0.5, // 24 -> 2 [0, 0]
    0.5, 1.0, 0.5, // 25 -> 3 [0, 0]
    0.0, 0.5, 0.5, // 26 -> 4 [0, 0]
    0.5, 0.5, 1.0, // 27 -> 5 [0, 0]
// at center on edges    
    0.0, 0.5, 0.0, // 28 -> 0 [0, -1] 
    0.5, 1.0, 0.0, // 29 -> 0 [1, 0] 
    1.0, 0.5, 0.0, // 30 -> 0 [0, 1]
    0.5, 0.0, 0.0, // 31 -> 0 [-1, 0] 
    1.0, 0.0, 0.5, // 32 -> 1 [1, 0] 
    0.5, 0.0, 1.0, // 33 -> 1 [0, 1] 
    0.0, 0.0, 0.5, // 34 -> 1 [-1, 0] 
    1.0, 1.0, 0.5, // 35 -> 2 [1, 0] 
    1.0, 0.5, 1.0, // 36 -> 2 [0, 1] 
    0.0, 1.0, 0.5, // 37 -> 3 [1, 0] 
    0.5, 1.0, 1.0, // 38 -> 3 [0, 1] 
    0.0, 0.5, 1.0, // 39 -> 4 [0, 1] 
// off center on faces 
    0.25, 0.25, 0.0, // 40 -> 0 [-0.5, -0.5]
    0.75, 0.25, 0.0, // 41 -> 0 [-0.5, 0.5]
    0.75, 0.75, 0.0, // 42 -> 0 [0.5, 0.5]
    0.25, 0.75, 0.0, // 43 -> 0 [0.5, -0.5]
    0.25, 0.0, 0.25, // 44 -> 1 [-0.5, -0.5]
    0.75, 0.0, 0.25, // 45 -> 1 [0.5, -0.5]
    0.75, 0.0, 0.75, // 46 -> 1 [0.5, 0.5]
    0.25, 0.0, 0.75, // 47 -> 1 [-0.5, 0.5]
    1.0, 0.25, 0.25, // 48 -> 2 [-0.5, -0.5]
    1.0, 0.75, 0.25, // 49 -> 2 [0.5, -0.5]
    1.0, 0.75, 0.75, // 50 -> 2 [0.5, 0.5]
    1.0, 0.25, 0.75, // 51 -> 2 [-0.5, 0.5]
    0.25, 1.0, 0.25, // 52 -> 3 [0.5, -0.5]
    0.75, 1.0, 0.25, // 53 -> 3 [-0.5, -0.5]
    0.75, 1.0, 0.75, // 54 -> 3 [-0.5, 0.5]
    0.25, 1.0, 0.75, // 55 -> 3 [0.5, 0.5] 
    0.0, 0.25, 0.25, // 56 -> 4 [0.5, -0.5]
    0.0, 0.75, 0.25, // 57 -> 4 [-0.5, -0.5]
    0.0, 0.75, 0.75, // 58 -> 4 [-0.5, 0.5]
    0.0, 0.25, 0.75, // 59 -> 4 [0.5, 0.5]
    0.25, 0.25, 1.0, // 60 -> 5 [-0.5, -0.5]
    0.75, 0.25, 1.0, // 61 -> 5 [0.5, -0.5]
    0.75, 0.75, 1.0, // 62 -> 5 [0.5, 0.5]
    0.25, 0.75, 1.0, // 63 -> 5 [-0.5, 0.5]
// off center on edges    
// ...    
  };

  unsigned int face[] = { 0, 0, 0, 0, 1, 1, 2, 3, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, };
  double coordinates[] = { -1, -1, -1, 1, 1, 1, 1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 1, -1, 0, 1, 0, 0, 1, -1, 0, 1, 0, 0, 1, 1, 0, 0, 1,  0, 1, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5, -0.5, -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5,  0.5, -0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, };

  for (unsigned int i = 0; i < n_candidates; ++i) {
    std::pair<unsigned int, std::vector<double> > out = volume_element.project_on_face(candidates[3*i+0], candidates[3*i+1], candidates[3*i+2]);
    AMP_ASSERT(out.first == face[i]);
    AMP_ASSERT(out.second[0] == coordinates[2*i+0]);
    AMP_ASSERT(out.second[1] == coordinates[2*i+1]);
  } // end for i
}

void test_mapping(hex8_element_t *volume_element, unsigned int n_random_candidate_points = 10000, double abs_tol = 1.0e-12, double rel_tol = 1.0e-12) {
  std::vector<double> random_candidate_point(3), candidate_point_global_coordinates(3), candidate_point_local_coordinates(3), mapping_error(3);
  double mapping_error_norm, tol;
  for (unsigned int i = 0; i < n_random_candidate_points; ++i) {
    for (unsigned int j = 0; j < 3; ++j) { random_candidate_point[j] = -1.0+2.0*rand()/RAND_MAX; }
    volume_element->map_local_to_global(&(random_candidate_point[0]), &(candidate_point_global_coordinates[0]));
    volume_element->map_global_to_local(&(candidate_point_global_coordinates[0]), &(candidate_point_local_coordinates[0]));

    for (unsigned int i = 0; i < 3; ++i) { mapping_error[i] = candidate_point_local_coordinates[i] - random_candidate_point[i]; }
    mapping_error_norm = sqrt(std::inner_product(&(mapping_error[0]), &(mapping_error[0])+3, &(mapping_error[0]), 0.0));
    tol = abs_tol + rel_tol * sqrt(std::inner_product(&(random_candidate_point[0]), &(random_candidate_point[0])+3, &(random_candidate_point[0]), 0.0));
    AMP_ASSERT(mapping_error_norm < tol);
  } // end for i

}

void write_point(double const * p, std::ostream & os = std::cout) {
  os<<"("<<p[0]<<","<<p[1]<<","<<p[2]<<")";
}

void write_cycle(unsigned int n, double const * * c, std::ostream & os = std::cout) {
  for (unsigned int i = 0; i < n; ++i) {
    write_point(c[i], os);
    os<<" -- ";
  } // end for i
  os<<"cycle ;\n";
}

void write_face(double const * * f, std::ostream & os = std::cout) {
  write_cycle(4, f, os);
}

void write_triangle(double const * * t, std::ostream & os = std::cout) {
  write_cycle(3, t, os);
}

void draw_face(hex8_element_t * e, unsigned int f, const std::string & option = "", std::ostream & os = std::cout) {
  assert(f < 6);
  std::vector<double const *> sp_ptr(4);
  unsigned int const * f_ptr = e->get_face(f);
  for (unsigned int p = 0; p < 4; ++p) { sp_ptr[p] = e->get_support_point(f_ptr[p]); } 
  os<<"\\draw["<<option<<"]\n";
  write_face(&(sp_ptr[0]), os);
}

void draw_triangles_on_face(hex8_element_t * e, unsigned int f, unsigned int t, const std::string & option = "", std::ostream & os = std::cout) {
  assert(f < 6);
  assert(t < 2);
  triangle_t * t_ptr = e->get_bounding_polyhedron();
  std::vector<double const *> sp_ptr(3);
  for (unsigned int i = 0; i < 3; ++i) { sp_ptr[i] = (t_ptr+2*f+t)->get_support_point_ptr(i); }
  os<<"\\draw["<<option<<"]\n";
  write_triangle(&(sp_ptr[0]), os);
}

void draw_bounding_box(hex8_element_t * e, double const * point_of_view, std::ostream & os = std::cout) {
  double const *bb_ptr = e->get_bounding_box();
  std::string plane[6] = { "yx", "zx", "zy", "zx", "zy", "yx" };
  std::string fixed_coord[6] = { "z", "y", "x", "y", "x", "z" };
  unsigned int fixed_coord_index[6] = { 2, 1, 3, 4, 0, 5 };

  unsigned int first_point_second_coord_index[6] = { 0, 0, 1, 0, 1, 0 }; 
  unsigned int first_point_first_coord_index[6] = { 1, 2, 2, 2, 2, 1 }; 
  unsigned int second_point_second_coord_index[6] = { 3, 3, 4, 3, 4, 3 }; 
  unsigned int second_point_first_coord_index[6] = { 4, 5, 5, 5, 5, 4 }; 

  double normals[18] = {
    0.0, 0.0, -1.0,
    0.0, -1.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    -1.0, 0.0, 0.0,
    0.0, 0.0, 1.0
  };

  os<<"\\tikzset{facestyle/.style={fill=none,draw=black,line join=round}}\n";
  for (unsigned int f = 0; f < 6; ++f) {
    os<<"\\begin{scope}[canvas is "<<plane[f]<<" plane at "<<fixed_coord[f]<<"="<<bb_ptr[fixed_coord_index[f]]<<"]\n";
    os<<"\\path[facestyle"<<((compute_scalar_product(point_of_view, normals+3*f) > 0.0) ? "" : ",dotted")<<"] ";
    os<<"("<<bb_ptr[first_point_first_coord_index[f]]<<","<<bb_ptr[first_point_second_coord_index[f]]<<")";
    os<<" rectangle ";
    os<<"("<<bb_ptr[second_point_first_coord_index[f]]<<","<<bb_ptr[second_point_second_coord_index[f]]<<")";
    os<<" ;\n";
    os<<"\\end{scope}\n";
  } // end for f
}  

std::string rubiks_cube_color_arrangement[6] = { "orange", "green", "white", "blue", "yellow", "red" };

void draw_bounding_polyhedron(hex8_element_t * e, double const * point_of_view, std::ostream & os = std::cout) {
  os<<"\\tikzset{facestyle/.style={opacity=0.4,line join=round}}\n";
  std::vector<std::string> options(12, "facestyle,");
  triangle_t * t_ptr = e->get_bounding_polyhedron();
  for (unsigned int f = 0; f < 6; ++f) { 
    for (unsigned int t = 0; t < 2; ++t) {
      if (compute_scalar_product(point_of_view, (t_ptr+2*f+t)->get_normal()) > 0.0) {
        options[2*f+t] += "fill=" + rubiks_cube_color_arrangement[f];
      } else {
        options[2*f+t] += "fill=none,dotted";
      } // end if
      draw_triangles_on_face(e, f, t, options[2*f+t]);
    } // end for t
  } // end for f
}

void draw_hex8_element(hex8_element_t * e, double const * point_of_view, std::ostream & os = std::cout) {
  os<<"\\tikzset{facestyle/.style={opacity=0.4,line join=round}}\n";
  std::vector<std::string> options(6, "facestyle,");
  triangle_t * t_ptr = e->get_bounding_polyhedron();
  for (unsigned int f = 0; f < 6; ++f) { 
    if (compute_scalar_product(point_of_view, (t_ptr+2*f)->get_normal()) > 0.0) {
      options[f] += "fill=" + rubiks_cube_color_arrangement[f];
    } else {
      options[f] += "fill=none,dotted";
    } // end if
    draw_face(e, f, options[f]); 
  } 
}

void draw_point(double const * p, std::string option = "", std::ostream & os = std::cout) {
  os<<"\\node["<<option<<"] at ";
  write_point(p, os);
  os<<" {+} ;\n";
}

void myTest(AMP::UnitTest *ut, std::string exeName) {
  old_test_project_on_face();

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

  hex8_element_t volume_element(std::vector<double>(points, points+24));

/*
  double scaling_factors[3] = { 4.0, 2.0, 1.0 };
  scale_points(std::vector<double>(scaling_factors, scaling_factors+3), 8, points);
*/

  double translation_vector[3] = { 3.0, 1.0, 5.0 };
  translate_points(std::vector<double>(translation_vector, translation_vector+3), 8, points);

  rotate_points(2, M_PI/3.0, 8, points);

  rotate_points(0, 0.75*M_PI, 8, points);

  srand(0);
  for (unsigned int i = 0; i < 24; ++i) { points[i] += -0.1 + 0.2*rand()/RAND_MAX; }

  volume_element.set_support_points(points);

  test_mapping(&volume_element);

  translate_points(1, 5.0, 8, points);
  translate_points(0, 2.0, 8, points);
  volume_element.set_support_points(points);

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
