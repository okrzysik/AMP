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

#include "operators/subchannel/SubchannelPhysicsModel.h"
#include "operators/subchannel/SubchannelOperatorParameters.h"
#include "operators/subchannel/SubchannelTwoEqNonlinearOperator.h"
#include "operators/subchannel/SubchannelConstants.h"
#include "operators/OperatorBuilder.h"

#include "ampmesh/StructuredMeshHelper.h"
#include "discretization/simpleDOF_Manager.h"

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

  // get dof manager
  int DofsPerFace =  2;
  AMP::Discretization::DOFManager::shared_ptr faceDOFManager = AMP::Discretization::simpleDOFManager::create( subchannelMesh,
     AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,1), AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,0), DofsPerFace);

  // get input and output variables
  AMP::LinearAlgebra::Variable::shared_ptr inputVariable  (new AMP::LinearAlgebra::Variable("flow"));
  AMP::LinearAlgebra::Variable::shared_ptr outputVariable (new AMP::LinearAlgebra::Variable("flow"));

  // create solution, rhs, and residual vectors
  AMP::LinearAlgebra::Vector::shared_ptr SolVec = AMP::LinearAlgebra::createVector( faceDOFManager , inputVariable , true );
  AMP::LinearAlgebra::Vector::shared_ptr RhsVec = AMP::LinearAlgebra::createVector( faceDOFManager , outputVariable, true );
  AMP::LinearAlgebra::Vector::shared_ptr ResVec = AMP::LinearAlgebra::createVector( faceDOFManager , outputVariable, true );

  // create subchannel physics model
  boost::shared_ptr<AMP::Database> subchannelPhysics_db = input_db->getDatabase("SubchannelPhysicsModel");
  boost::shared_ptr<AMP::Operator::ElementPhysicsModelParameters> params( new AMP::Operator::ElementPhysicsModelParameters(subchannelPhysics_db));
  boost::shared_ptr<AMP::Operator::SubchannelPhysicsModel>  subchannelPhysicsModel (new AMP::Operator::SubchannelPhysicsModel(params));

  // get nonlinear operator database
  boost::shared_ptr<AMP::Database> subchannelOperator_db = input_db->getDatabase("SubchannelTwoEqNonlinearOperator");
  // set operator parameters
  boost::shared_ptr<AMP::Operator::SubchannelOperatorParameters> subchannelOpParams(new AMP::Operator::SubchannelOperatorParameters( subchannelOperator_db ));
  subchannelOpParams->d_Mesh = subchannelMesh ;
  subchannelOpParams->d_subchannelPhysicsModel = subchannelPhysicsModel;
  subchannelOpParams->d_dofMap = faceDOFManager ;

  // create nonlinear operator
  boost::shared_ptr<AMP::Operator::SubchannelTwoEqNonlinearOperator> subchannelOperator (new AMP::Operator::SubchannelTwoEqNonlinearOperator(subchannelOpParams));

  // report successful creation
  ut->passes(exeName+": creation");
  std::cout.flush();

  // reset the nonlinear operator
  subchannelOperator->reset(subchannelOpParams);
  std::vector<size_t> dofs;
  AMP::Mesh::MeshIterator face     = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
  AMP::Mesh::MeshIterator end_face = face.end();

  const double h_scale = AMP::Operator::Subchannel::scaleEnthalpy;
  const double P_scale = AMP::Operator::Subchannel::scalePressure;
  {
    // Test apply with known residual evaluation
    face = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], h_scale*1000.0e3 );
    SolVec->setValueByGlobalID(dofs[1], P_scale*16.4e6 );
    ++face;
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], h_scale*900.0e3 );
    SolVec->setValueByGlobalID(dofs[1], P_scale*16.3e6 );
    ++face;
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], h_scale*800.0e3 );
    SolVec->setValueByGlobalID(dofs[1], P_scale*16.2e6 );
    ++face;
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], h_scale*700.0e3 );
    SolVec->setValueByGlobalID(dofs[1], P_scale*16.1e6 );
    ++face;
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], h_scale*300.0e3 );
    SolVec->setValueByGlobalID(dofs[1], P_scale*13.5e6 );
    ++face;
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], h_scale*450.0e3 );
    SolVec->setValueByGlobalID(dofs[1], P_scale*9.0e6 );
    ++face;
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], h_scale*570.0e3 );
    SolVec->setValueByGlobalID(dofs[1], P_scale*12.0e5 );
    ++face;
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], h_scale*230.0e2 );
    SolVec->setValueByGlobalID(dofs[1], P_scale*4.0e6 );
    ++face;
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], h_scale*999.9e3 );
    SolVec->setValueByGlobalID(dofs[1], P_scale*14.0e6 );
    ++face;
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], h_scale*235.6e3 );
    SolVec->setValueByGlobalID(dofs[1], P_scale*12.5e6 );

    subchannelOperator->apply(RhsVec, SolVec, ResVec, 1.0, 0.0);
    bool passedKnownTest = true;
    double known[20] = {
       -316282.816245409,
       -49816.4072925864,
       -105719.954578781,
       -50981.1590877952,
       -116469.952796604,
       -52019.8509439426,
       -125233.43163654,
       -2555060.31658177,
       -430953.386215321,
       -4454165.84966422,
       117060.094406793,
       -7752956.92558845,
       89046.6137846786,
       2843912.28207917,
       -572233.43163654,
       10049463.6715359,
       960430.047203396,
       -1455636.07309586,
       -770019.954578781,
       -3013200
    };
    face     = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
    int i=0;
    for( ; face != end_face; ++face,++i){
      faceDOFManager->getDOFs( face->globalID(), dofs );
      double h_val = ResVec->getValueByGlobalID(dofs[0])/h_scale;
      double p_val = ResVec->getValueByGlobalID(dofs[1])/P_scale;
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
    SolVec->setValueByGlobalID(dofs[0], h_scale*950.0e3 );
    SolVec->setValueByGlobalID(dofs[1], P_scale*15.0e6 );
    ++face;
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], h_scale*850.0e3 );
    SolVec->setValueByGlobalID(dofs[1], P_scale*15.1e6 );
    ++face;
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], h_scale*700.0e3 );
    SolVec->setValueByGlobalID(dofs[1], P_scale*15.25e6 );
    ++face;
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], h_scale*500.0e3 );
    SolVec->setValueByGlobalID(dofs[1], P_scale*15.26e6 );
    ++face;
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], h_scale*324.6e3 );
    SolVec->setValueByGlobalID(dofs[1], P_scale*11.0e5 );
    ++face;
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], h_scale*457.7e3 );
    SolVec->setValueByGlobalID(dofs[1], P_scale*12.5e5 );
    ++face;
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], h_scale*134.6e2 );
    SolVec->setValueByGlobalID(dofs[1], P_scale*34.5e5 );
    ++face;
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], h_scale*457.6e3 );
    SolVec->setValueByGlobalID(dofs[1], P_scale*12.0e6 );
    ++face;
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], h_scale*325.7e3 );
    SolVec->setValueByGlobalID(dofs[1], P_scale*11.5e6 );
    ++face;
    faceDOFManager->getDOFs( face->globalID(), dofs );
    SolVec->setValueByGlobalID(dofs[0], h_scale*898.6e3 );
    SolVec->setValueByGlobalID(dofs[1], P_scale*15.7e6 );

    subchannelOperator->apply(RhsVec, SolVec, ResVec, 1.0, 0.0);
    bool passedKnownTest = true;
    double known[20] = {
       -367603.556722071,
       149631.268567802,
       -105719.954578781,
       198039.001082933,
       -166469.952796604,
       56265.1563678032,
       -225233.43163654,
       -14114645.1297306,
       -206353.386215321,
       196109.460143209,
       100160.094406793,
       2244037.06984132,
       -475193.386215321,
       8595470.36966505,
       418906.56836346,
       -454851.582252684,
       -148369.952796604,
       4249245.34387701,
       567180.045421219,
       186800
    };
    face     = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
    int i=0;
    for( ; face != end_face; ++face,++i){
      faceDOFManager->getDOFs( face->globalID(), dofs );
      double h_val = ResVec->getValueByGlobalID(dofs[0])/h_scale;
      double p_val = ResVec->getValueByGlobalID(dofs[1])/P_scale;
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
    ++face;
    SolVec->setValueByGlobalID(dofs[0], h_scale*700.0e3);
    SolVec->setValueByGlobalID(dofs[1], P_scale*12.4e6);
    faceDOFManager->getDOFs( face->globalID(), dofs );
    ++face;
    SolVec->setValueByGlobalID(dofs[0], h_scale*900.0e3);
    SolVec->setValueByGlobalID(dofs[1], P_scale*12.3e6);
    faceDOFManager->getDOFs( face->globalID(), dofs );
    ++face;
    SolVec->setValueByGlobalID(dofs[0], h_scale*800.0e3);
    SolVec->setValueByGlobalID(dofs[1], P_scale*16.2e6);
    faceDOFManager->getDOFs( face->globalID(), dofs );
    ++face;
    SolVec->setValueByGlobalID(dofs[0], h_scale*650.0e3);
    SolVec->setValueByGlobalID(dofs[1], P_scale*14.1e5);
    faceDOFManager->getDOFs( face->globalID(), dofs );
    ++face;
    SolVec->setValueByGlobalID(dofs[0], h_scale*367.4e3);
    SolVec->setValueByGlobalID(dofs[1], P_scale*31.5e5);
    faceDOFManager->getDOFs( face->globalID(), dofs );
    ++face;
    SolVec->setValueByGlobalID(dofs[0], h_scale*657.2e3);
    SolVec->setValueByGlobalID(dofs[1], P_scale*12.5e6);
    faceDOFManager->getDOFs( face->globalID(), dofs );
    ++face;
    SolVec->setValueByGlobalID(dofs[0], h_scale*788.5e3);
    SolVec->setValueByGlobalID(dofs[1], P_scale*12.7e6);
    faceDOFManager->getDOFs( face->globalID(), dofs );
    ++face;
    SolVec->setValueByGlobalID(dofs[0], h_scale*235.7e2);
    SolVec->setValueByGlobalID(dofs[1], P_scale*17.8e6);
    faceDOFManager->getDOFs( face->globalID(), dofs );
    ++face;
    SolVec->setValueByGlobalID(dofs[0], h_scale*673.1e3);
    SolVec->setValueByGlobalID(dofs[1], P_scale*13.6e6);
    faceDOFManager->getDOFs( face->globalID(), dofs );
    ++face;
    SolVec->setValueByGlobalID(dofs[0], h_scale*385.2e3);
    SolVec->setValueByGlobalID(dofs[1], P_scale*16.3e6);

    subchannelOperator->apply(RhsVec, SolVec, ResVec, 1.0, 0.0);
    bool passedKnownTest = true;
    double known[20] = {
       -620085.843930095,
       -49986.4688121105,
       194280.045421219,
       3949028.44119528,
       -116469.952796604,
       -14741988.3976949,
       -175233.43163654,
       1785589.98266971,
       -313553.386215321,
       9397345.85264865,
       256860.094406793,
       248813.913712698,
       100346.613784679,
       5143446.55072186,
       -790163.43163654,
       -4153289.69480504,
       633060.047203396,
       2745414.89283861,
       -293619.954578781,
       786800
    };
    face     = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
    int i=0;
    for( ; face != end_face; ++face,++i){
      faceDOFManager->getDOFs( face->globalID(), dofs );
      double h_val = ResVec->getValueByGlobalID(dofs[0])/h_scale;
      double p_val = ResVec->getValueByGlobalID(dofs[1])/P_scale;
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


