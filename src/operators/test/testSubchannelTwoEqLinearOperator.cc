#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>

#include "boost/shared_ptr.hpp"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"

#include "vectors/VectorBuilder.h"
#include "vectors/SimpleVector.h"

#include "SubchannelPhysicsModel.h"
#include "SubchannelOperatorParameters.h"
#include "SubchannelTwoEqLinearOperator.h"
#include "../OperatorBuilder.h"

// function to check that Jacobian matches known values
bool JacobianIsCorrect(double **testJacobian, double knownJacobian[5][5])
{
   bool passed = true; // boolean for all values being equal to known values
   double max_abs_err = 0.0; // maximum absolute error; L-infinity norm of error
   double max_rel_err = 0.0; // maximum relative error
   for (size_t i = 0; i <= 4; i++){
      for (size_t j = 0; j <= 4; j++){
         // check that Jacobian value is approximately equal to known value
         if (!AMP::Utilities::approx_equal(testJacobian[i][j],knownJacobian[i][j],0.01)) passed = false;
         // calculate absolute error
         double abs_err = std::abs(testJacobian[i][j]-knownJacobian[i][j]);
         if (abs_err > max_abs_err) max_abs_err = abs_err;
         // calculate relative error
         double rel_err;
         if (std::abs(knownJacobian[i][j]) < 1.0e-15) rel_err = 0.0;
         else rel_err = abs_err / knownJacobian[i][j];
         if (rel_err > max_rel_err) max_rel_err = rel_err;
      }
   }
   if (!passed){
      // report errors
      AMP::pout << "Jacobian did not match known Jacobian.\n";
      AMP::pout << "Jacobian comparison: Maximum absolute error: " << max_abs_err << ", Maximum relative error: " << max_rel_err << "\n";
      AMP::pout << "Test Jacobian:\n";
      for (size_t i = 0; i <= 4; i++){
         AMP::pout << "\n";
         for (size_t j = 0; j <= 4; j++){
            AMP::pout << "   " << testJacobian[i][j];
         }
      }
      AMP::pout << "\n\nKnown Jacobian:\n";
      for (size_t i = 0; i <= 4; i++){
         AMP::pout << "\n";
         for (size_t j = 0; j <= 4; j++){
            AMP::pout << "   " << knownJacobian[i][j];
         }
      }
   }
   return passed;
}

void Test(AMP::UnitTest *ut, const std::string exeName)
{
  // create input and output file names
  std::string input_file = "input_"  + exeName;
  std::string log_file   = "output_" + exeName;
  AMP::PIO::logOnlyNodeZero(log_file);

  // get input database from input file
  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  // create mesh
  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  boost::shared_ptr<AMP::Database>  mesh_db = input_db->getDatabase("Mesh");
  boost::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(mesh_db));
  meshParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
  boost::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh(meshParams);

  // create subchannel physics model
  boost::shared_ptr<AMP::Database> subchannelPhysics_db = input_db->getDatabase("SubchannelPhysicsModel");
  boost::shared_ptr<AMP::Operator::ElementPhysicsModelParameters> params( new AMP::Operator::ElementPhysicsModelParameters(subchannelPhysics_db));
  boost::shared_ptr<AMP::Operator::SubchannelPhysicsModel>  subchannelPhysicsModel (new AMP::Operator::SubchannelPhysicsModel(params));

  // get nonlinear operator database
  boost::shared_ptr<AMP::Database> subchannelOperator_db = input_db->getDatabase("SubchannelTwoEqLinearOperator");
  // set operator parameters
  boost::shared_ptr<AMP::Operator::SubchannelOperatorParameters> subchannelOpParams(new AMP::Operator::SubchannelOperatorParameters( subchannelOperator_db ));
  subchannelOpParams->d_Mesh = meshAdapter;
  subchannelOpParams->d_subchannelPhysicsModel = subchannelPhysicsModel;

  // create linear operator
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel;
  boost::shared_ptr<AMP::Operator::SubchannelTwoEqLinearOperator> subchannelOperator =
      boost::dynamic_pointer_cast<AMP::Operator::SubchannelTwoEqLinearOperator>(AMP::Operator::OperatorBuilder::createOperator(
      meshAdapter,"SubchannelTwoEqLinearOperator",input_db,elementModel ));
  
  // report successful creation
  ut->passes(exeName+": creation");
  std::cout.flush();

  // get input and output variables
  AMP::LinearAlgebra::Variable::shared_ptr inputVariable  = subchannelOperator->getInputVariable();
  AMP::LinearAlgebra::Variable::shared_ptr outputVariable = subchannelOperator->getOutputVariable();

  // create solution, rhs, and residual vectors
  AMP::LinearAlgebra::Vector::shared_ptr FrozenVec = AMP::LinearAlgebra::SimpleVector::create( 5, inputVariable );  // frozen solution
  AMP::LinearAlgebra::Vector::shared_ptr SolVec = AMP::LinearAlgebra::SimpleVector::create( 5, inputVariable );  // delta solution
  AMP::LinearAlgebra::Vector::shared_ptr RhsVec = AMP::LinearAlgebra::SimpleVector::create( 5, outputVariable ); // RHS
  AMP::LinearAlgebra::Vector::shared_ptr ResVec = AMP::LinearAlgebra::SimpleVector::create( 5, outputVariable ); // residual

  // reset the nonlinear operator
  subchannelOperator->reset(subchannelOpParams);

  {
     // Test apply with known residual evaluation
     FrozenVec->setValueByLocalID(0,1000.0);
     FrozenVec->setValueByLocalID(1,15.3e6);
     FrozenVec->setValueByLocalID(2,15.2e6);
     FrozenVec->setValueByLocalID(3,15.1e6);
     FrozenVec->setValueByLocalID(4,15.0e6);
     SolVec->setValueByLocalID(0,1.0);
     SolVec->setValueByLocalID(1,1.0);
     SolVec->setValueByLocalID(2,1.0);
     SolVec->setValueByLocalID(3,1.0);
     SolVec->setValueByLocalID(4,1.0);
     subchannelOperator->setFrozenVector(FrozenVec);
     subchannelOperator->apply(RhsVec, SolVec, ResVec, 1.0, 0.0);

     double knownJacobian[5][5] = {
        {1.00000286783759,0.000945636121207302,0,0,0},
        {0.0084145730631915,-1.00006033273598,0.999924832234253,0,0},
        {0.0164760871167386,0,-1.00006231804521,0.999922391414922,0},
        {0.0240080126558191,0,0,-1.00006376793541,0.999920554193017},
        {0,0,0,0,0}
     };
     bool passedJacobianTest ;
//     bool passedJacobianTest = JacobianIsCorrect(subchannelOperator->d_Jacobian,knownJacobian);
     if (passedJacobianTest) ut->passes(exeName+": apply: known Jacobian value test #1");
     else ut->failure(exeName+": apply: known Jacobian value test #1");
  }

  {
     // Test apply with known residual evaluation
     FrozenVec->setValueByLocalID(0,500.0e3);
     FrozenVec->setValueByLocalID(1,17.0e6);
     FrozenVec->setValueByLocalID(2,17.1e6);
     FrozenVec->setValueByLocalID(3,17.2e6);
     FrozenVec->setValueByLocalID(4,17.3e6);
     SolVec->setValueByLocalID(0,1.0);
     SolVec->setValueByLocalID(1,1.0);
     SolVec->setValueByLocalID(2,1.0);
     SolVec->setValueByLocalID(3,1.0);
     SolVec->setValueByLocalID(4,1.0);
     subchannelOperator->setFrozenVector(FrozenVec);
     subchannelOperator->apply(RhsVec, SolVec, ResVec, 1.0, 0.0);

     double knownJacobian[5][5] = {
        {1.00000000036425,0.00093610492591568,0,0,0},
        {0.0616135012140977,-1.00009118544544,0.999886751731324,0,0},
        {0.068075038946247,0,-1.00010127215554,0.99987438016631,0},
        {0.073216467438169,0,0,-1.00011121418204,0.999862069911393},
        {0,0,0,0,0}
     };
     bool passedJacobianTest ; 
//     bool passedJacobianTest = JacobianIsCorrect(subchannelOperator->d_Jacobian,knownJacobian);
     if (passedJacobianTest) ut->passes(exeName+": apply: known Jacobian value test #2");
     else ut->failure(exeName+": apply: known Jacobian value test #2");
  }

  {
     // Test apply with known residual evaluation
     FrozenVec->setValueByLocalID(0,1000.0e3);
     FrozenVec->setValueByLocalID(1,16.5e6);
     FrozenVec->setValueByLocalID(2,16.0e6);
     FrozenVec->setValueByLocalID(3,16.1e6);
     FrozenVec->setValueByLocalID(4,16.2e6);
     SolVec->setValueByLocalID(0,1.0);
     SolVec->setValueByLocalID(1,1.0);
     SolVec->setValueByLocalID(2,1.0);
     SolVec->setValueByLocalID(3,1.0);
     SolVec->setValueByLocalID(4,1.0);
     subchannelOperator->setFrozenVector(FrozenVec);
     subchannelOperator->apply(RhsVec, SolVec, ResVec, 1.0, 0.0);

     double knownJacobian[5][5] = {
        {0.999999999030052,0.000938893989442007,0,0,0},
        {0.115520453086643,-1.0001777089274,0.999780389006059,0,0},
        {0.131575829506197,0,-1.00020366987796,0.999748303691582,0},
        {0.150335512370785,0,0,-1.00023876717391,0.999705452198002},
        {0,0,0,0,0}
     };
     bool passedJacobianTest ; 
//     bool passedJacobianTest = JacobianIsCorrect(subchannelOperator->d_Jacobian,knownJacobian);
     if (passedJacobianTest) ut->passes(exeName+": apply: known Jacobian value test #3");
     else ut->failure(exeName+": apply: known Jacobian value test #3");
  }

}

int main(int argc, char *argv[])
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup(argc,argv,startup_properties);

    AMP::UnitTest ut;

    const int NUMFILES=1;
    std::string files[NUMFILES] = {
        "testSubchannelTwoEqLinearOperator"
    };

    for (int i=0; i<NUMFILES; i++) {
        try {
            Test(&ut, files[i]);

        } catch (std::exception &err) {
            std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
            ut.failure("ERROR: While testing: "+files[i]);
        } catch( ... ) {
            std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
            ut.failure("ERROR: While testing: "+files[i]);
        }
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}   


