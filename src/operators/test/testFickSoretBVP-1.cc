#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>

#include "utils/shared_ptr.h"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"

#include "ampmesh/Mesh.h"
#include "vectors/Variable.h"
#include "vectors/MultiVariable.h"
#include "vectors/VectorBuilder.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"

#include "libmesh/libmesh.h"

#include "operators/diffusion/FickSoretNonlinearFEOperator.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"

#include "operators/boundary/DirichletMatrixCorrection.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/boundary/libmesh/NeumannVectorCorrection.h"

#include "../BVPOperatorParameters.h"
#include "../LinearBVPOperator.h"
#include "../NonlinearBVPOperator.h"
#include "../OperatorBuilder.h"

#include "applyTests.h"


void bvpTest1(AMP::UnitTest *ut, std::string exeName)
{
  // Tests FickSoret Dirchlet BVP operator for temperature

  // Initialization
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  // Input database
  AMP::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

//--------------------------------------------------
//   Create the Mesh.
//--------------------------------------------------
  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  AMP::shared_ptr<AMP::Database>  mesh_db = input_db->getDatabase("Mesh");
  AMP::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(new AMP::Mesh::MeshParameters(mesh_db));
  mgrParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
  AMP::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh(mgrParams);
//--------------------------------------------------

  // Create nonlinear FickSoret BVP operator and access volume nonlinear FickSoret operator
  AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
  AMP::shared_ptr<AMP::Operator::Operator> nlinBVPOperator =
    AMP::Operator::OperatorBuilder::createOperator(meshAdapter, "testFickSoretBVPOperator", input_db, elementPhysicsModel);
  AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nlinBVPOp =
          AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(nlinBVPOperator);
  AMP::shared_ptr<AMP::Operator::FickSoretNonlinearFEOperator> nlinOp =
         AMP::dynamic_pointer_cast<AMP::Operator::FickSoretNonlinearFEOperator>(nlinBVPOp->getVolumeOperator());
  AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> fickOp =
         AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(nlinOp->getFickOperator());
  AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> soretOp =
         AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(nlinOp->getSoretOperator());
  
  // use the linear BVP operator to create a Fick linear operator with bc's
  AMP::shared_ptr<AMP::Operator::Operator> linBVPOperator =
    AMP::Operator::OperatorBuilder::createOperator(meshAdapter, "testLinearFickBVPOperator", input_db, elementPhysicsModel);
  AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linBVPOp =
          AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(linBVPOperator);

  ut->passes(exeName+": creation");
  std::cout.flush();

  // set shift, scale for applyTests
  double shift[2];
  double scale[2];
  shift[0] = 0.;
  shift[1] = 0.;
  scale[0] = 1.;
  scale[1] = 1.;
  std::vector<double> trange(2), crange(2);
  AMP::shared_ptr<AMP::Operator::DiffusionTransportModel> fickTransportModel = fickOp->getTransportModel();
  std::vector<double> defaults(2,0);
  AMP::Materials::Material::shared_ptr fmat = fickTransportModel->getMaterial();
  // the Soret has a principal variable of temperature
  if (soretOp->getPrincipalVariableId() == AMP::Operator::Diffusion::TEMPERATURE) {
      std::string property="ThermalDiffusionCoefficient";
      if( (fmat->property(property))->is_argument("temperature") ) {
          trange = (fmat->property(property))->get_arg_range("temperature");  // Compile error
          scale[1] = trange[1]-trange[0];
          shift[1] = trange[0]+0.001*scale[1];
          scale[1] *= 0.999;
          defaults = (fmat->property(property))->get_defaults();  // compile error
      }
  }
  // the Fick has a principal variable of temperature
  if (fickOp->getPrincipalVariableId() == AMP::Operator::Diffusion::CONCENTRATION) {
      std::string property="FickCoefficient";
      if( (fmat->property(property))->is_argument("concentration") ) {
          crange = (fmat->property(property))->get_arg_range("concentration");  // Compile error
          scale[0] = crange[1]-crange[0];
          shift[0] = crange[0]+0.001*scale[0];
          scale[0] *= 0.999;
          defaults = (fmat->property(property))->get_defaults();  // compile error
      }
  }

  // Set up input and output vectors
  AMP::LinearAlgebra::Variable::shared_ptr cVar, tVar;
  cVar = AMP::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVariable>(fickOp->getInputVariable())->getVariable(AMP::Operator::Diffusion::CONCENTRATION);
  tVar = AMP::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVariable>(soretOp->getInputVariable())->getVariable(AMP::Operator::Diffusion::TEMPERATURE);
  AMP::shared_ptr<AMP::LinearAlgebra::MultiVariable> fsInpVar(new AMP::LinearAlgebra::MultiVariable("fsInput"));
  fsInpVar->add(cVar);
  fsInpVar->add(tVar);

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // Create a DOF manager for a nodal vector 
  int DOFsPerNode = 1;
  int nodalGhostWidth = 1;
  bool split = true;
  AMP::Discretization::DOFManager::shared_ptr nodalDofMap = AMP::Discretization::simpleDOFManager::create(meshAdapter, AMP::Mesh::Vertex, nodalGhostWidth, DOFsPerNode, split);
  //----------------------------------------------------------------------------------------------------------------------------------------------//

  // create solution, rhs, and residual vectors
  AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> multivector = AMP::LinearAlgebra::MultiVector::create ( "mulitvector", meshAdapter->getComm() );
  multivector->addVector( AMP::LinearAlgebra::createVector( nodalDofMap, cVar ) );
  multivector->addVector( AMP::LinearAlgebra::createVector( nodalDofMap, tVar ) );
  AMP::LinearAlgebra::Vector::shared_ptr solVec = multivector;
  AMP::LinearAlgebra::Vector::shared_ptr rhsVec = solVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr resVec = solVec->cloneVector();

  // set default values of input variables
  AMP::LinearAlgebra::Vector::shared_ptr inConcVec = solVec->subsetVectorForVariable(cVar);
  AMP::LinearAlgebra::Vector::shared_ptr inTempVec = solVec->subsetVectorForVariable(tVar);
  inConcVec->setToScalar(defaults[1]);  // compile error
  inTempVec->setToScalar(defaults[0]);  // compile error
  rhsVec->setToScalar(0.);

  AMP_INSIST(nlinOp->isValidInput(solVec), "input variable not set up correctly");

  // Test apply
  std::string msgPrefix=exeName+": apply";
  {
    applyTests(ut, msgPrefix, nlinBVPOperator, rhsVec, solVec, resVec, shift, scale, 2);
  }
  std::cout.flush();

  // Test linear reset from getJacobianParameters
  for(int i = 0; i < 3; i++) {
    inConcVec->setRandomValues();
    inTempVec->setRandomValues();
    adjust(solVec, shift, scale, 2);
    rhsVec->setRandomValues();
    resVec->setRandomValues();
    AMP::shared_ptr<AMP::Operator::OperatorParameters> jacparams = nlinBVPOp->getParameters("Jacobian", solVec);
    linBVPOp->reset(jacparams);
  }//end for i

  ut->passes(exeName+": getJacobianParameters");
  std::cout.flush();

}

int main(int argc, char *argv[])
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup(argc,argv,startup_properties);

    AMP::UnitTest ut;

  const int NUMFILES=2;
  std::string files[NUMFILES] = {
        "FickSoret-BVP-TUI-1", "FickSoret-BVP-UO2MSRZC09-1"
  };

    for (int i=0; i<NUMFILES; i++) 
        bvpTest1(&ut, files[i]);


    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}   


