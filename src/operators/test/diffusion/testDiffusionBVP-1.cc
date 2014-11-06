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
#include "vectors/VectorBuilder.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"

#include "libmesh/libmesh.h"

#include "operators/diffusion/DiffusionNonlinearElement.h"
#include "operators/diffusion/DiffusionLinearElement.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"

#include "operators/boundary/DirichletMatrixCorrection.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/boundary/libmesh/NeumannVectorCorrection.h"

#include "operators/BVPOperatorParameters.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/OperatorBuilder.h"

#include "../applyTests.h"



void bvpTest1(AMP::UnitTest *ut, std::string exeName)
{
  // Tests diffusion Dirchlet BVP operator for temperature

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

  // Create nonlinear Diffusion BVP operator and access volume nonlinear Diffusion operator
  AMP::shared_ptr<AMP::InputDatabase> nbvp_db =
          AMP::dynamic_pointer_cast<AMP::InputDatabase>( input_db->getDatabase("ThermalNonlinearBVPOperator"));
  AMP::shared_ptr<AMP::Operator::Operator> nlinBVPOperator =
    AMP::Operator::OperatorBuilder::createOperator(meshAdapter,"ThermalNonlinearBVPOperator",input_db);
  AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nlinBVPOp =
          AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(nlinBVPOperator);
  AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> nlinOp =
         AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(nlinBVPOp->getVolumeOperator());
  AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel = nlinOp->getTransportModel();
  
  // use the linear BVP operator to create a thermal linear operator with bc's
  AMP::shared_ptr<AMP::Operator::Operator> linBVPOperator =
    AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
						   "ThermalLinearBVPOperator",
						   input_db,
						   elementPhysicsModel);
  AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linBVPOp =
          AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(linBVPOperator);
  //AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> linOp =
  //        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(linBVPOp->getVolumeOperator());

  ut->passes(exeName+": creation");
  std::cout.flush();

  // Set up input and output vectors
  AMP::LinearAlgebra::Variable::shared_ptr bvpSolVar = nlinOp->getOutputVariable();
  AMP::LinearAlgebra::Variable::shared_ptr bvpRhsVar = nlinOp->getOutputVariable();
  AMP::LinearAlgebra::Variable::shared_ptr bvpResVar = nlinOp->getOutputVariable();

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // Create a DOF manager for a nodal vector 
  int DOFsPerNode = 1;
  int nodalGhostWidth = 1;
  bool split = true;
  AMP::Discretization::DOFManager::shared_ptr nodalDofMap = AMP::Discretization::simpleDOFManager::create(meshAdapter, AMP::Mesh::Vertex, nodalGhostWidth, DOFsPerNode, split);
  //----------------------------------------------------------------------------------------------------------------------------------------------//

  // create solution, rhs, and residual vectors
  AMP::LinearAlgebra::Vector::shared_ptr bvpSolVec = AMP::LinearAlgebra::createVector( nodalDofMap, bvpSolVar );
  AMP::LinearAlgebra::Vector::shared_ptr bvpRhsVec = AMP::LinearAlgebra::createVector( nodalDofMap, bvpRhsVar );
  AMP::LinearAlgebra::Vector::shared_ptr bvpResVec = AMP::LinearAlgebra::createVector( nodalDofMap, bvpResVar );

  bvpRhsVec->setToScalar(0.0);

  AMP::shared_ptr<AMP::Database> volOp_db =  input_db->getDatabase(nbvp_db->getString("VolumeOperator"));
  AMP::shared_ptr<AMP::Database> model_db = input_db->getDatabase(volOp_db->getString("LocalModel"));
  std::string property=model_db->getString("Property");

  // set shift, scale for applyTests
  double shift=0., scale=1.;
  std::vector<double> range(2);
  AMP::shared_ptr<AMP::Operator::DiffusionTransportModel> transportModel =
          AMP::dynamic_pointer_cast<AMP::Operator::DiffusionTransportModel>(elementPhysicsModel);
  AMP::Materials::Material::shared_ptr mat = transportModel->getMaterial();
  if (nlinOp->getPrincipalVariableId() == AMP::Operator::Diffusion::TEMPERATURE) {
      if( (mat->property(property))->is_argument("temperature") ) {
          range = (mat->property(property))->get_arg_range("temperature");  // Compile error
          scale = range[1]-range[0];
          shift = range[0]+0.001*scale;
          scale *= 0.999;
      }
  }
  if (nlinOp->getPrincipalVariableId() == AMP::Operator::Diffusion::CONCENTRATION) {
      if( (mat->property(property))->is_argument("concentration") ) {
          range = (mat->property(property))->get_arg_range("concentration");  // Compile error
          scale = range[1]-range[0];
          shift = range[0]+0.001*scale;
          scale *= 0.999;
      }
  }
  if (nlinOp->getPrincipalVariableId() == AMP::Operator::Diffusion::BURNUP) {
      AMP_INSIST(false, "do not know what to do");
  }

  // Test apply
  std::string msgPrefix=exeName+": apply";
  {
    applyTests(ut, msgPrefix, nlinBVPOperator, bvpRhsVec, bvpSolVec, bvpResVec, shift, scale);
  }
  std::cout.flush();

  // Test linear reset from getJacobianParameters
  for(int i = 0; i < 3; i++) {
    bvpSolVec->setRandomValues();
    adjust(bvpSolVec, shift, scale);
    bvpRhsVec->setRandomValues();
    bvpResVec->setRandomValues();
    AMP::shared_ptr<AMP::Operator::OperatorParameters> jacparams = nlinBVPOp->getJacobianParameters(bvpSolVec);
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

    const int NUMFILES=6;
    std::string files[NUMFILES] = {
        "Diffusion-TUI-Thermal-2", "Diffusion-TUI-Fick-2", "Diffusion-TUI-Soret-2",
        "Diffusion-UO2MSRZC09-Thermal-2", "Diffusion-UO2MSRZC09-Fick-2", "Diffusion-UO2MSRZC09-Soret-2"//,
        //"Diffusion-TUI-Fick-3", "Diffusion-UO2MSRZC09-Fick-3"
    };

    for (int i=0; i<NUMFILES; i++) {
        try {
            bvpTest1(&ut, files[i]);

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


