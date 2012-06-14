#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
//#include "utils/Utilities.h"
//#include "utils/Database.h"
//#include "utils/InputDatabase.h"
//#include "utils/InputManager.h"
//#include "utils/AMP_MPI.h"
//#include "utils/PIO.h"

//#include "ampmesh/Mesh.h"

#include "hex8_element_t.h"

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
    assert(out.first == face[i]);
    assert(out.second[0] == coordinates[2*i+0]);
    assert(out.second[1] == coordinates[2*i+1]);
  } // end for i


  const unsigned int n_random_candidates = 10000;
  srand(0);
  // scaling 
  for (unsigned int i = 0; i < 8; ++i) { 
    points[3*i+0] *= 4.0; 
    points[3*i+1] *= 2.0; 
    points[3*i+2] *= 1.0; 
  } // end for i
  // translating and rotating the support points 
  for (unsigned int i = 0; i < 8; ++i) { 
    points[3*i+0] += 3.0; 
    points[3*i+1] += 1.0; 
    points[3*i+2] += 5.0; 
  } // end for i
  double theta = M_PI/3.0;
  for (unsigned int i = 0; i < 8; ++i) { 
    double tmp[3];
    tmp[0] = cos(theta)*points[3*i+0]-sin(theta)*points[3*i+1];
    tmp[1] = sin(theta)*points[3*i+0]+cos(theta)*points[3*i+1];
    tmp[2] = points[3*i+2];
    std::copy(tmp, tmp+3, points+3*i);
  } // end for i
  theta = -0.75*M_PI;
  for (unsigned int i = 0; i < 8; ++i) { 
    double tmp[3];
    tmp[0] = points[3*i+0];
    tmp[1] = cos(theta)*points[3*i+1]-sin(theta)*points[3*i+2];
    tmp[2] = sin(theta)*points[3*i+1]+cos(theta)*points[3*i+2];
    std::copy(tmp, tmp+3, points+3*i);
  } // end for i
  // moving them around
  for (unsigned int i = 0; i < 24; ++i) { points[i] += -0.1 + 0.2*rand()/RAND_MAX; }
//  for (unsigned int i = 0; i < 24; ++i) { std::cout<<points[i]<<"  "; if (i%3 == 2) {std::cout<<"\n";} }
  volume_element.set_support_points(std::vector<double>(points, points+24));
  std::vector<double> bounding_box = volume_element.get_bounding_box();

  unsigned int count = 0;
  for (unsigned int i = 0; i < n_random_candidates; ++i) {
    std::vector<double> candidate(3);
//    for (unsigned int j = 0; j < 3; ++j) { candidate[j] = (j==2 ? -10.0 : -5.0) + 10.0*rand()/RAND_MAX; }
    for (unsigned int j = 0; j < 3; ++j) { candidate[j] = bounding_box[j+0]+(bounding_box[j+3]-bounding_box[j+0])*rand()/RAND_MAX; }
    std::vector<double> tmp = volume_element.map_global_to_local(candidate);
//    for (unsigned int j = 0; j < 3; ++j) { assert(candidate[j] == 0.5*(tmp[j]+1.0)); }
    std::cout<<i<<"  ";
    for (unsigned int j = 0; j < 3; ++j) { std::cout<<i<<"  "<<candidate[j]<<"  "<<0.5*(tmp[j]+1.0)<<"\n"; }
    std::cout<<"\n";
    volume_element.do_mapping_verification_test(candidate);
    if (volume_element.project_on_face(candidate).first != 99) { ++count; }
  } // end for i
//  std::cout<<"count="<<count<<"\n";
  assert(count != 0);

  ut->passes(exeName);
}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::string exeName = "testHex8ElementContactSearch";

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
