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

//#include "ampmesh/Mesh.h"
#include "vectors/VectorBuilder.h"
#include "vectors/SimpleVector.h"

#include "SubchannelPhysicsModel.h"
#include "SubchannelOperatorParameters.h"
#include "SubchannelTwoEqNonlinearOperator.h"
#include "../OperatorBuilder.h"
#include "discretization/simpleDOF_Manager.h"

#include "SubchannelHelpers.h"

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
  boost::shared_ptr<AMP::Mesh::Mesh> subchannelMesh = AMP::Mesh::Mesh::buildMesh(meshParams);
  AMP::Mesh::Mesh::shared_ptr xyFaceMesh;
  xyFaceMesh = subchannelMesh->Subset( getFaceIterator( subchannelMesh , 1 ) );

  // create subchannel physics model
  boost::shared_ptr<AMP::Database> subchannelPhysics_db = input_db->getDatabase("SubchannelPhysicsModel");
  boost::shared_ptr<AMP::Operator::ElementPhysicsModelParameters> params( new AMP::Operator::ElementPhysicsModelParameters(subchannelPhysics_db));
  boost::shared_ptr<AMP::Operator::SubchannelPhysicsModel>  subchannelPhysicsModel (new AMP::Operator::SubchannelPhysicsModel(params));

  // get nonlinear operator database
  boost::shared_ptr<AMP::Database> subchannelOperator_db = input_db->getDatabase("SubchannelTwoEqNonlinearOperator");
  // set operator parameters
  boost::shared_ptr<AMP::Operator::SubchannelOperatorParameters> subchannelOpParams(new AMP::Operator::SubchannelOperatorParameters( subchannelOperator_db ));
  subchannelOpParams->d_Mesh = xyFaceMesh ;
  subchannelOpParams->d_subchannelPhysicsModel = subchannelPhysicsModel;

  // create nonlinear operator
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel;
  boost::shared_ptr<AMP::Operator::SubchannelTwoEqNonlinearOperator> subchannelOperator =
      boost::dynamic_pointer_cast<AMP::Operator::SubchannelTwoEqNonlinearOperator>(AMP::Operator::OperatorBuilder::createOperator(
      xyFaceMesh ,"SubchannelTwoEqNonlinearOperator",input_db,elementModel ));
  
  // report successful creation
  ut->passes(exeName+": creation");
  std::cout.flush();

  int DofsPerFace =  2;
  AMP::Discretization::DOFManager::shared_ptr faceDOFManager = AMP::Discretization::simpleDOFManager::create( xyFaceMesh, AMP::Mesh::Face, 1, DofsPerFace, true);

  // get input and output variables
  AMP::LinearAlgebra::Variable::shared_ptr inputVariable  = subchannelOperator->getInputVariable();
  AMP::LinearAlgebra::Variable::shared_ptr outputVariable = subchannelOperator->getOutputVariable();

  // create solution, rhs, and residual vectors
  AMP::LinearAlgebra::Vector::shared_ptr SolVec = AMP::LinearAlgebra::createVector( faceDOFManager , inputVariable , true );
  AMP::LinearAlgebra::Vector::shared_ptr RhsVec = AMP::LinearAlgebra::createVector( faceDOFManager , outputVariable, true );
  AMP::LinearAlgebra::Vector::shared_ptr ResVec = AMP::LinearAlgebra::createVector( faceDOFManager , outputVariable, true );

  // reset the nonlinear operator
  subchannelOperator->reset(subchannelOpParams);
  std::vector<size_t> dofs;
  AMP::Mesh::MeshIterator face     = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
  AMP::Mesh::MeshIterator end_face = face.end();

  {
     // Test apply with known residual evaluation
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], 1000 );
    double j = 15.3;
    for( ; face != end_face; ++face,j=j-0.1){
      faceDOFManager->getDOFs( face->globalID(), dofs );
      SolVec->setValueByGlobalID(dofs[1], j*1.e6);
    }
    subchannelOperator->apply(RhsVec, SolVec, ResVec, 1.0, 0.0);
    bool passedKnownTest = true;
    double known[8] = {-1.3166e6,2.4471e5, 0.0,2.4563e5, 0.0,7.602e5, 0.0,0.0};
    face     = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
    int i=0;
    for( ; face != end_face; ++face,++i){
      faceDOFManager->getDOFs( face->globalID(), dofs );
      double val = ResVec->getValueByGlobalID(dofs[1]);
      if (!AMP::Utilities::approx_equal(val,known[i],0.01)) passedKnownTest = false;
    }
    if (passedKnownTest) ut->passes(exeName+": known value test #1");
    else ut->failure(exeName+": apply: known residual test #1");
  }

  face     = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
  {
     // Test apply with known residual evaluation
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], 500e3 );
    double j = 17.0;
    for( ; face != end_face; ++face,j=j+0.1){
      faceDOFManager->getDOFs( face->globalID(), dofs );
      SolVec->setValueByLocalID(dofs[1],j*1.0e6);
    } 
    subchannelOperator->apply(RhsVec, SolVec, ResVec, 1.0, 0.0);
    bool passedKnownTest = true;
    double known[8] = {-8.15720e5,4.62164e5, 0.0,4.66929e5, 0.0,-1.31480e6, 0.0,0.0};
    face     = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
    int i=0;
    for( ; face != end_face; ++face,++i){
      faceDOFManager->getDOFs( face->globalID(), dofs );
      double val = ResVec->getValueByGlobalID(dofs[1]);
      if (!AMP::Utilities::approx_equal(val,known[i],0.01)) passedKnownTest = false;
    }
     if (passedKnownTest) ut->passes(exeName+": known value test #2");
     else ut->failure(exeName+": apply: known residual test #2");
  }

  face     = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
  {
     // Test apply with known residual evaluation
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], 1000e3 );
    double j = 16.5;
    for( ; face != end_face; ++face,j=j-0.1){
      faceDOFManager->getDOFs( face->globalID(), dofs );
      SolVec->setValueByLocalID(dofs[1],j*1.0e6);
    } 
    subchannelOperator->apply(RhsVec, SolVec, ResVec, 1.0, 0.0);
    bool passedKnownTest = true;
    double known[8] = {-3.161889e5,3.04554e5, 0.0,3.13654e5, 0.0,-3.63256e5, 0.0,0.0};
    face     = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
    int i=0;
    for( ; face != end_face; ++face,++i){
      faceDOFManager->getDOFs( face->globalID(), dofs );
      double val = ResVec->getValueByGlobalID(dofs[1]);
      if (!AMP::Utilities::approx_equal(val,known[i],0.01)) passedKnownTest = false;
    }
    if (passedKnownTest) ut->passes(exeName+": known value test #3");
    else ut->failure(exeName+": apply: known residual test #3");
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
        "testSubchannelTwoEqNonlinearOperator"
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


