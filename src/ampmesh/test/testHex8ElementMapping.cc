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

void draw_face(hex8_element_t * e, unsigned int f, std::string option = "", std::ostream & os = std::cout) {
  std::vector<double const *> sp_ptr(4);
  unsigned int const * f_ptr = e->get_face(f);
  for (unsigned int p = 0; p < 4; ++p) { sp_ptr[p] = e->get_support_point(f_ptr[p]); } 
  os<<"\\draw[line join=round"<<option<<"]\n";
  write_face(&(sp_ptr[0]), os);
}

void draw_triangles_on_face(hex8_element_t * e, unsigned int f, std::string option = "", std::ostream & os = std::cout) {
  triangle_t const * t_ptr = e->get_bounding_polyhedron();
  for (unsigned int i = 0; i < f; ++i) { ++t_ptr; ++t_ptr; }
  std::vector<double const *> sp_ptr(3);
  ++t_ptr;
  for (unsigned int i = 0; i < 3; ++i) { sp_ptr[i] = t_ptr->get_support_point_ptr(i); }
  os<<"\\draw[line join=round"<<option<<"]\n";
  write_triangle(&(sp_ptr[0]), os);
  --t_ptr;
  for (unsigned int i = 0; i < 3; ++i) { sp_ptr[i] = t_ptr->get_support_point_ptr(i); }
  os<<"\\draw[line join=round"<<option<<"]\n";
  write_triangle(&(sp_ptr[0]), os);
}

void draw_bounding_box(hex8_element_t * e, std::ostream & os = std::cout) {
  double const *bb_ptr = e->get_bounding_box();
  std::string path_front = "\\path[facestyle]";
  std::string path_behind = "\\path[facestyle,dotted]";
  os<<"\\tikzset{facestyle/.style={fill=none,draw=black,line join=round}}\n";
  os<<"\\begin{scope}[canvas is zy plane at x="<<bb_ptr[0]<<"]\n";
  os<<path_behind<<" ("<<bb_ptr[2]<<","<<bb_ptr[1]<<") rectangle ("<<bb_ptr[5]<<","<<bb_ptr[4]<<") ;\n";
  os<<"\\end{scope}\n";
  os<<"\\begin{scope}[canvas is zx plane at y="<<bb_ptr[1]<<"]\n";
  os<<path_behind<<" ("<<bb_ptr[2]<<","<<bb_ptr[0]<<") rectangle ("<<bb_ptr[5]<<","<<bb_ptr[3]<<") ;\n";
  os<<"\\end{scope}\n";
  os<<"\\begin{scope}[canvas is zx plane at y="<<bb_ptr[4]<<"]\n";
  os<<path_front<<" ("<<bb_ptr[2]<<","<<bb_ptr[0]<<") rectangle ("<<bb_ptr[5]<<","<<bb_ptr[3]<<") ;\n";
  os<<"\\end{scope}\n";
  os<<"\\begin{scope}[canvas is zy plane at x="<<bb_ptr[3]<<"]\n";
  os<<path_front<<" ("<<bb_ptr[2]<<","<<bb_ptr[1]<<") rectangle ("<<bb_ptr[5]<<","<<bb_ptr[4]<<") ;\n";
  os<<"\\end{scope}\n";
  os<<"\\begin{scope}[canvas is yx plane at z="<<bb_ptr[5]<<"]\n";
  os<<path_front<<" ("<<bb_ptr[1]<<","<<bb_ptr[0]<<") rectangle ("<<bb_ptr[4]<<","<<bb_ptr[3]<<") ;\n";
  os<<"\\end{scope}\n";
}

/*void draw_bounding_polyhedron(hex8_element_t * e, std::ostream & os = std::cout) {
  for (unsigned int f = 0; f < 6; ++f) { 
    draw_triangle_on_face(e, f, 0, options[f], 
    os<<"\n";
  } // end for f
}*/

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

  // rubik's cube color arrangement
  std::string options[] = { ",fill=orange", ",fill=green", ",fill=white", ",fill=blue", ",fill=yellow", ",fill=red" };

  // draw_hex8_element(&volume_element)
  std::cout<<"% volume element\n";
  for (unsigned int f = 0; f < 6; ++f) { draw_face(&volume_element, f, options[f]); } 

  // draw_bounding_polyhedron(&volume_element);
  std::cout<<"% bounding polyhedron\n";
  for (unsigned int f = 0; f < 6; ++f) { 
    draw_triangles_on_face(&volume_element, f, options[f]);
  } // end for f

  std::cout<<"% bounding box\n";
  draw_bounding_box(&volume_element);

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
