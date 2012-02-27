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

#include "ampmesh/Mesh.h"

#include "operators/diffusion/DiffusionTransportModel.h"
#include "operators/diffusion/DiffusionLinearElement.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionLinearFEOperatorParameters.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "operators/diffusion/DiffusionConstants.h"
#include "operators/OperatorParameters.h"
#include "operators/OperatorBuilder.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"


void linearTest1(AMP::UnitTest *ut, std::string exeName)
{
  // this tests creation from database and usage

  // Test create
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  std::string mesh_file = input_db->getString("Mesh");

  // Create the mesh parameter object
  boost::shared_ptr<AMP::MemoryDatabase> database(new AMP::MemoryDatabase("Mesh"));
  database->putInteger("dim",3);
  database->putString("MeshName","mesh");
  database->putString("MeshType","libMesh");
  database->putString("FileName",mesh_file);
  boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
  params->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));

  // Create the mesh
  AMP::Mesh::Mesh::shared_ptr  meshAdapter = AMP::Mesh::Mesh::buildMesh(params);

  boost::shared_ptr<AMP::Operator::DiffusionLinearFEOperator> diffOp;
  boost::shared_ptr<AMP::InputDatabase> diffLinFEOp_db =
          boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("LinearDiffusionOp"));

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel;
  diffOp = boost::dynamic_pointer_cast<AMP::Operator::DiffusionLinearFEOperator>(
										 AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
																"LinearDiffusionOp",
																input_db,
																elementModel));

  AMP::LinearAlgebra::Variable::shared_ptr diffSolVar = diffOp->getInputVariable();
  AMP::LinearAlgebra::Variable::shared_ptr diffRhsVar = diffOp->getOutputVariable();
  AMP::LinearAlgebra::Variable::shared_ptr diffResVar = diffOp->getOutputVariable();

  AMP::Discretization::DOFManager::shared_ptr NodalScalarDOF = AMP::Discretization::simpleDOFManager::create(meshAdapter,AMP::Mesh::Vertex,1,1,true);

  AMP::LinearAlgebra::Vector::shared_ptr diffSolVec = AMP::LinearAlgebra::createVector( NodalScalarDOF, diffSolVar, true );
  AMP::LinearAlgebra::Vector::shared_ptr diffRhsVec = AMP::LinearAlgebra::createVector( NodalScalarDOF, diffRhsVar, true );
  AMP::LinearAlgebra::Vector::shared_ptr diffResVec = AMP::LinearAlgebra::createVector( NodalScalarDOF, diffResVar, true );


  ut->passes(exeName);

  // Test apply
  for(int i = 0; i < 10; i++) {
    diffSolVec->setRandomValues();
    diffRhsVec->setRandomValues();
    diffResVec->setRandomValues();
    diffOp->apply(diffRhsVec, diffSolVec, diffResVec, 1.0, -1.0);
  }//end for i

  ut->passes(exeName);

  // Test reset
  boost::shared_ptr<AMP::Operator::DiffusionLinearFEOperatorParameters> diffOpParams(new
      AMP::Operator::DiffusionLinearFEOperatorParameters( diffLinFEOp_db ));
  diffOpParams->d_transportModel = boost::dynamic_pointer_cast<AMP::Operator::DiffusionTransportModel>(elementModel);
  diffOp->reset(diffOpParams);

  ut->passes(exeName);

  // Test eigenvalues (run output through mathematica)
  boost::shared_ptr<AMP::LinearAlgebra::Matrix> diffMat = diffOp->getMatrix();

  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
  int nranks = globalComm.getSize();

  size_t matdim=24;
  if(nranks==1) {
      std::cout << "cols={"<<std::endl;
      for(size_t i = 0; i < matdim; i++) {
        std::vector<unsigned int> matCols;
        std::vector<double> matVals;
        diffMat->getRowByGlobalID(i, matCols, matVals);
        std::cout <<"{";
        for (size_t j = 0; j < matCols.size(); j++) {
            std::cout<<matCols[j];
            if (j<matCols.size()-1) std::cout << ",";
        }
        std::cout <<"}";
        if (i<matdim-1) std::cout<<",";
        std::cout <<std::endl;
      }
      std::cout << "};"<<std::endl;

      std::cout<<"matrix = {"<<std::endl;

      for(size_t i = 0; i < matdim; i++) {
          std::vector<unsigned int> matCols;
          std::vector<double> matVals;
          diffMat->getRowByGlobalID(i, matCols, matVals);
          std::cout<<"{";
          size_t col=0;
          for(size_t j = 0; j < matCols.size(); j++) {
            while (col<matCols[j]) {
                std::cout << "0.";
                std::cout<<",";
                col++;
            }
            std::cout<<matVals[j];
            if (matCols[j]<matdim-1) std::cout<<",";
            col++;
          }//end for j
          while (col<matdim) {
            std::cout << "0";
            if (col<matdim-1)std::cout<<",";
            col++;
          }
          std::cout<<"}";
          if (i<matdim-1) std::cout<<","<<std::endl;
      }//end for i

      std::cout<<"};"<<std::endl;
  }

  ut->passes(exeName);
}

// void linearTestReset(AMP::UnitTest *ut, std::string exeName)
// {
//   // this tests creation from database and usage
// 
//   // Test create
//   std::string input_file = "input_" + exeName;
//   std::string log_file = "output_" + exeName;
// 
//   AMP::PIO::logOnlyNodeZero(log_file);
// 
//   boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
//   AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
//   input_db->printClassData(AMP::plog);
// 
//   AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
//   std::string mesh_file = input_db->getString("Mesh");
// 
//   //create a mesh adaptor
//   AMP::MeshManager::Adapter::shared_ptr meshAdapter =
//               AMP::MeshManager::Adapter::shared_ptr ( new AMP::MeshManager::Adapter () );
//   meshAdapter->readExodusIIFile ( mesh_file.c_str() );
// 
//   // create the linear diffusion operator
//   boost::shared_ptr<AMP::Operator::DiffusionLinearFEOperator> diffOp;
//   boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel;
//   boost::shared_ptr<AMP::InputDatabase> diffLinFEOp_db =
//           boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("LinearDiffusionOp"));
//   diffOp = boost::dynamic_pointer_cast<AMP::Operator::DiffusionLinearFEOperator>(
//                                        AMP::Operator::OperatorBuilder::createOperator(meshAdapter, diffLinFEOp_db, elementModel));
// 
//   // first test: reset with a NULL parameter object
//   boost::shared_ptr<OperatorParameters> resetParameters;
//   bool passed=false;
//   try
//     {
//       diffOp->reset(resetParameters);
//     }
//   catch(std::exceptions)
//     {
//       passed=true;
//     }
//   if(passed)
//     {
//       ut.passes(exeName+": DiffusionLinearFEOperator::reset with NULL parameter object");
//     }
//   else
//     {
//       ut.numFails++;
//     }
// 
//   // second test: create a non null parameter object but dont initialize it fully  passed=false;
//   try
//     {
//       diffOp->reset(resetParameters);
//     }
//   catch(std::exceptions)
//     {
//       passed=true;
//     }
//   if(passed)
//     {
//       ut.passes(exeName+": DiffusionLinearFEOperator::reset with parameter object having no valid physics and operation objects");
//     }
//   else
//     {
//       ut.numFails++;
//     }
// 
// 
//   boost::shared_ptr<AMP::Operator::DiffusionLinearFEOperatorParameters> diffusionOpParams(new AMP::Operator::DiffusionLinearFEOperatorParameters( diffLinFEOp_db ));
//   resetParameters = boost::dynamic_pointer_cast<OperatorParameters>(diffusionOpParams);
//   AMP_INSIST(resetParameters.get()!=NULL, "unable to create parameters");
//   
//   passed=false;
//   try
//     {
//       diffOp->reset(resetParameters);
//     }
//   catch(std::exceptions)
//     {
//       passed=true;
//     }
//   if(passed)
//     {
//       ut.passes(exeName+": DiffusionLinearFEOperator::reset with parameter object having no valid physics and operation objects");
//     }
//   else
//     {
//       ut.numFails++;
//     }
// 
//   // third test: use a parameter object without a element operation object
//   boost::shared_ptr<AMP::Database> transportModel_db = input_db->getDatabase("DiffusionTransportModel");
//   boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel = ElementPhysicsModelFactory::createElementPhysicsModel(transportModel_db);
//   boost::shared_ptr<AMP::Operator::DiffusionTransportModel> transportModel = boost::dynamic_pointer_cast<DiffusionTransportModel>(elementPhysicsModel);
//   AMP_INSIST(transportModel.get()!=NULL, "unable to create transport model");
//   diffusionOpParams->d_transportModel = transportModel;
// 
//   passed=false;
//   try
//     {
//       diffOp->reset(resetParameters);
//     }
//   catch(std::exceptions)
//     {
//       passed=true;
//     }
//   if(passed)
//     {
//       ut.passes(exeName+": DiffusionLinearFEOperator::reset with parameter object having no valid element operation objects");
//     }
//   else
//     {
//       ut.numFails++;
//     }
// 
//   // next create a ElementOperation object
//   AMP_INSIST(input_db->keyExists("DiffusionElement"), "Key ''DiffusionElement'' is missing!");
//   boost::shared_ptr<AMP::Operator::ElementOperation> diffusionLinElem = ElementOperationFactory::createElementOperation(input_db->getDatabase("DiffusionElement"));
// 
// }

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

  const int NUMFILES=8;
  std::string files[NUMFILES] = {
          "Diffusion-TUI-Thermal-1", "Diffusion-TUI-Fick-1", "Diffusion-TUI-Soret-1",
          "Diffusion-UO2MSRZC09-Thermal-1", "Diffusion-UO2MSRZC09-Fick-1", "Diffusion-UO2MSRZC09-Soret-1",
          "Diffusion-TUI-TensorFick-1", "Diffusion-CylindricalFick-1"
  };

    try {
    for (int i=0; i<NUMFILES; i++) {
        linearTest1(&ut, files[i]);    }
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



