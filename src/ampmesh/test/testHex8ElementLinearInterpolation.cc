#include "utils/AMPManager.h"
#include "utils/UnitTest.h"

#include "hex8_element_t.h"

void scale_points(unsigned int direction, double scaling_factor, unsigned int n_points, double* points) {
  assert(direction < 3);
  for (unsigned int i = 0; i < n_points; ++i) { 
    points[3*i+direction] *= scaling_factor; 
  } // end for i
}

void scale_points(const std::vector<double> &scaling_factors, unsigned int n_points, double* points) {
  assert(scaling_factors.size() == 3);
  for (unsigned int i = 0; i < 3; ++i) { 
    scale_points(i, scaling_factors[i], n_points, points);
  } // end for i
}

void translate_points(unsigned int direction, double distance, unsigned int n_points, double* points) {
  assert(direction < 3);
  for (unsigned int i = 0; i < n_points; ++i) { 
    points[3*i+direction] += distance; 
  } // end for i
}

void translate_points(std::vector<double> translation_vector, unsigned int n_points, double* points) {
  assert(translation_vector.size() == 3);
  for (unsigned int i = 0; i < 3; ++i) { 
   translate_points(i, translation_vector[i], n_points, points); 
  } // end for i
}

void rotate_points(unsigned int rotation_axis, double rotation_angle, unsigned int n_points, double* points) {
  assert(rotation_axis < 3);
  unsigned int non_fixed_directions[2];
  unsigned int i = 0;
  for (unsigned int j = 0; j < 3; ++j) { 
    if (j != rotation_axis) {
      non_fixed_directions[i++] = j;
    } // end if
  } // end for j
  double tmp[3];
  for (unsigned int j = 0; j < n_points; ++j) {
    tmp[non_fixed_directions[0]] = cos(rotation_angle)*points[3*j+non_fixed_directions[0]]-sin(rotation_angle)*points[3*j+non_fixed_directions[1]];
    tmp[non_fixed_directions[1]] = sin(rotation_angle)*points[3*j+non_fixed_directions[0]]+cos(rotation_angle)*points[3*j+non_fixed_directions[1]];
    tmp[rotation_axis] = points[3*j+rotation_axis];
    std::copy(tmp, tmp+3, points+3*j);
  } // end for j
}

double my_function(const std::vector<double> &xyz) {
  AMP_ASSERT(xyz.size() == 3);
  double x = xyz[0], y = xyz[1], z = xyz[2];
  return (1.0 + 6.0 * x) * (2.0 - 5.0 * y) * (3.0 + 4.0 * z);
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

  double scaling_factors[3] = { 4.0, 2.0, 1.0 };
  scale_points(std::vector<double>(scaling_factors, scaling_factors+3), 8, points);

  double translation_vector[3] = { 3.0, 1.0, 5.0 };
  translate_points(std::vector<double>(translation_vector, translation_vector+3), 8, points);

//  rotate_points(2, M_PI/3.0, 8, points);
//  rotate_points(0, -0.74*M_PI, 8, points);

  srand(0);
/*  // moving them around
  for (unsigned int i = 0; i < 24; ++i) { points[i] += -0.1 + 0.2*rand()/RAND_MAX; }
//  for (unsigned int i = 0; i < 24; ++i) { std::cout<<points[i]<<"  "; if (i%3 == 2) {std::cout<<"\n";} }
  volume_element.set_support_points(std::vector<double>(points, points+24));*/
  hex8_element_t volume_element(std::vector<double>(points, points+24));

  std::vector<double> my_function_at_support_points(8);
  for (unsigned int i = 0; i < 8; ++i) { my_function_at_support_points[i] = my_function(std::vector<double>(points+3*i, points+3*i+3)); }
  std::vector<double> bounding_box = volume_element.get_bounding_box();

  unsigned int count = 0;
  const unsigned int n_random_candidates = 10;
  for (unsigned int i = 0; i < n_random_candidates; ++i) {
    std::vector<double> candidate_global_coordinates(3);
    for (unsigned int j = 0; j < 3; ++j) { candidate_global_coordinates[j] = bounding_box[j+0]+(bounding_box[j+3]-bounding_box[j+0])*rand()/RAND_MAX; }

    std::vector<double> candidate_local_coordinates = volume_element.map_global_to_local(candidate_global_coordinates);
    bool coordinates_are_local = true;
    if (volume_element.contains_point(candidate_local_coordinates, coordinates_are_local)) { ++count; }

    std::vector<double> basis_functions_values = get_basis_functions_values(candidate_local_coordinates);
    double interpolated_value = 0.0;
    for (unsigned int j = 0; j < 8; ++j) {
      interpolated_value += my_function_at_support_points[j] * basis_functions_values[j];
    } // end for j

    double my_function_at_candidate = my_function(candidate_global_coordinates);
    double interpolation_error = fabs(interpolated_value-my_function_at_candidate);
    std::cout<<std::setprecision(15)<<"interpolated_value="<<interpolated_value<<"  my_function_at_candidate="<<my_function_at_candidate<<"  interpolation_error="<<interpolation_error<<"\n";
//    assert(std::min(fabs(interpolated_value-my_function_at_candidate), fabs((interpolated_value-my_function_at_candidate)/interpolated_value)) < 1.0e-13);
    
  } // end for i
//  std::cout<<"count="<<count<<"\n";
  assert(count != 0);
  
/*        hex8_element_t volume_element(dummy);
          double value = 0.0;
          for (unsigned int j = 0; j < 8; ++j) {
            std::vector<size_t> globalID;
            dofManager->getDOFs(support_points[j].globalID(), globalID);
            AMP_ASSERT(globalID.size() == 1);
            value += vectorField->getValueByGlobalID(globalID.front()) * basis_functions_values[j];
          } // end for j
*/
/*  DendroSearch dendroSearch(globalComm, meshAdapter);
  std::vector<double> interpolatedData; 
  std::vector<bool> interpolationWasDone;
  dendroSearch.interpolate(dummyVector, pts, interpolatedData, interpolationWasDone);
  AMP_ASSERT(interpolatedData.size() == numLocalPts);
  AMP_ASSERT(interpolationWasDone.size() == numLocalPts);

  int localNotFound = static_cast<int>(std::count(interpolationWasDone.begin(), interpolationWasDone.end(), false));
  int globalNotFound = -1;
  MPI_Allreduce(&localNotFound, &globalNotFound, 1, MPI_INT, MPI_SUM, globalComm.getCommunicator());
  if(!rank) {
    std::cout<<globalNotFound<<" points total weren't found"<<std::endl;
  }


  std::vector<double> interpolationError(numLocalPts, 0.0);
  for (unsigned int i = 0; i < numLocalPts; ++i) {
    if (interpolationWasDone[i]) {
      interpolationError[i] = fabs(interpolatedData[i] - dummyFunction(std::vector<double>(&(pts[3*i]), &(pts[3*i+3]))));
    }
  } // end for i
  double localErrorSquaredNorm = std::inner_product(interpolationError.begin(), interpolationError.end(), interpolationError.begin(), 0.0);
  double localMaxError = *std::max_element(interpolationError.begin(), interpolationError.end());
  if(!rank) {
    std::cout<<"Finished computing the local squared norm of the interpolation error."<<std::endl;
  }
  globalComm.barrier();


  double globalMaxError = globalComm.maxReduce<double>(localMaxError);
  double globalErrorSquaredNorm = -1.0;
  MPI_Allreduce(&localErrorSquaredNorm, &globalErrorSquaredNorm, 1, MPI_DOUBLE, MPI_SUM, globalComm.getCommunicator());
  if(!rank) {
    std::cout<<"Global error norm is "<<sqrt(globalErrorSquaredNorm)<<std::endl;
    std::cout<<"Global max error is "<<globalMaxError<<std::endl;
  }

//  AMP_ASSERT(sqrt(globalErrorSquaredNorm) < 1.0e-15);

*/

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



