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

#include "ampmesh/StructuredMeshHelper.h"

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
  xyFaceMesh = subchannelMesh->Subset( AMP::Mesh::StructuredMeshHelper::getXYFaceIterator( subchannelMesh , 0 ) );

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
    face = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], 1000.0e3 );
    double k = 10.0;
    double j = 16.4;
    for( ; face != end_face; ++face, j = j-0.1, k = k-1.0){
      faceDOFManager->getDOFs( face->globalID(), dofs );
      SolVec->setValueByGlobalID(dofs[0], k*1.e5);
      SolVec->setValueByGlobalID(dofs[1], j*1.e6);
    }
    subchannelOperator->apply(RhsVec, SolVec, ResVec, 1.0, 0.0);
    bool passedKnownTest = true;
    double known[8] = {
       -316282.816245409,
       48038.9385840241,
       -147423.339011925,
       44515.6919998378,
       -194846.67802385,
       41390.8375866283,
       -147423.339011925,
       586800
    };
    face     = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
    int i=0;
    for( ; face != end_face; ++face,++i){
      faceDOFManager->getDOFs( face->globalID(), dofs );
      double h_val = ResVec->getValueByGlobalID(dofs[0]);
      double p_val = ResVec->getValueByGlobalID(dofs[1]);
      if (!AMP::Utilities::approx_equal(h_val,known[2*i],0.01)){
         passedKnownTest = false;
         AMP::pout<<"Calculated: "<<h_val<<", Known: "<<known[2*i]<<"\n";
      }
      if (!AMP::Utilities::approx_equal(p_val,known[2*i+1],0.01)){
         passedKnownTest = false;
         AMP::pout<<"Calculated: "<<p_val<<", Known: "<<known[2*i+1]<<"\n";
      }
    }
    if (passedKnownTest) ut->passes(exeName+": known value test #1");
    else ut->failure(exeName+": known residual test #1");
  }

  {
    // Test apply with known residual evaluation
    face = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], 950.0e3 );
    SolVec->setValueByGlobalID(dofs[1], 15.0e6 );
    ++face;
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], 850.0e3 );
    SolVec->setValueByGlobalID(dofs[1], 15.1e6 );
    ++face;
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], 700.0e3 );
    SolVec->setValueByGlobalID(dofs[1], 15.25e6 );
    ++face;
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], 500.0e3 );
    SolVec->setValueByGlobalID(dofs[1], 15.26e6 );

    subchannelOperator->apply(RhsVec, SolVec, ResVec, 1.0, 0.0);
    bool passedKnownTest = true;
    double known[8] = {
       -367603.556722071,
       246371.994992077,
       -147423.339011925,
       292021.286966963,
       -244846.67802385,
       146982.197373965,
       -247423.339011925,
       -253200
    };
    face     = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
    int i=0;
    for( ; face != end_face; ++face,++i){
      faceDOFManager->getDOFs( face->globalID(), dofs );
      double h_val = ResVec->getValueByGlobalID(dofs[0]);
      double p_val = ResVec->getValueByGlobalID(dofs[1]);
      if (!AMP::Utilities::approx_equal(h_val,known[2*i],0.01)){
         passedKnownTest = false;
         AMP::pout<<"Calculated: "<<h_val<<", Known: "<<known[2*i]<<"\n";
      }
      if (!AMP::Utilities::approx_equal(p_val,known[2*i+1],0.01)){
         passedKnownTest = false;
         AMP::pout<<"Calculated: "<<p_val<<", Known: "<<known[2*i+1]<<"\n";
      }
    }
    if (passedKnownTest) ut->passes(exeName+": known value test #2");
    else ut->failure(exeName+": known residual test #2");
  }

  {
    // Test apply with known residual evaluation
    face = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], 700.0e3 );
    SolVec->setValueByGlobalID(dofs[1], 12.4e6 );
    ++face;
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], 900.0e3 );
    SolVec->setValueByGlobalID(dofs[1], 12.3e6 );
    ++face;
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], 800.0e3 );
    SolVec->setValueByGlobalID(dofs[1], 16.2e6 );
    ++face;
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], 650.0e3 );
    SolVec->setValueByGlobalID(dofs[1], 14.1e5 );

    subchannelOperator->apply(RhsVec, SolVec, ResVec, 1.0, 0.0);
    bool passedKnownTest = true;
    double known[8] = {
       -620085.843930095,
       44727.2351558046,
       152576.660988075,
       4044682.60608771,
       -194846.67802385,
       -14648547.7970646,
       -197423.339011925,
       -14103200
    };
    face     = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
    int i=0;
    for( ; face != end_face; ++face,++i){
      faceDOFManager->getDOFs( face->globalID(), dofs );
      double h_val = ResVec->getValueByGlobalID(dofs[0]);
      double p_val = ResVec->getValueByGlobalID(dofs[1]);
      if (!AMP::Utilities::approx_equal(h_val,known[2*i],0.01)){
         passedKnownTest = false;
         AMP::pout<<"Calculated: "<<h_val<<", Known: "<<known[2*i]<<"\n";
      }
      if (!AMP::Utilities::approx_equal(p_val,known[2*i+1],0.01)){
         passedKnownTest = false;
         AMP::pout<<"Calculated: "<<p_val<<", Known: "<<known[2*i+1]<<"\n";
      }
    }
    if (passedKnownTest) ut->passes(exeName+": known value test #3");
    else ut->failure(exeName+": known residual test #3");
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


