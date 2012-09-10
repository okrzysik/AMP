
#include "utils/AMPManager.h"
#include "utils/InputManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>
#include "materials/Material.h"

#include "utils/AMPManager.h"
#include "utils/ReadTestMesh.h"
#include "utils/WriteSolutionToFile.h"

#include "ampmesh/libmesh/libMesh.h"
#include "ampmesh/SiloIO.h"
#include "ampmesh/Mesh.h"
#include "ampmesh/MultiMesh.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"

#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"

#include "operators/ColumnOperator.h"
#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/mechanics/MechanicsLinearFEOperator.h"

#include "solvers/TrilinosMLSolver.h"

#include "mesh_communication.h"

#include "MPCtestUtils.h"

void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

   AMP::PIO::logOnlyNodeZero(log_file);
   AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

   boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
   AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
   input_db->printClassData(AMP::plog);
   std::string mesh_file1 = input_db->getString("mesh_file1");

   boost::shared_ptr<AMP::Mesh::initializeLibMesh> libmeshInit( new AMP::Mesh::initializeLibMesh(globalComm) );
   {

      const unsigned int mesh_dim = 3;
      boost::shared_ptr< ::Mesh > meshMaster(new ::Mesh(mesh_dim));
      boost::shared_ptr< ::Mesh > meshSlave(new ::Mesh(mesh_dim));

      bool binaryMeshes = input_db->getBool("BinaryMeshes");

      if(binaryMeshes) {
        AMP::readBinaryTestMesh(mesh_file1, meshMaster);
        AMP::readBinaryTestMesh(mesh_file1, meshSlave);
      } else {
        AMP::readTestMesh(mesh_file1, meshMaster);
        AMP::readTestMesh(mesh_file1, meshSlave);
      }

      MeshCommunication().broadcast(*(meshMaster.get()));
      MeshCommunication().broadcast(*(meshSlave.get()));

      meshMaster->prepare_for_use(false);
      meshSlave->prepare_for_use(false);

      AMP::Mesh::Mesh::shared_ptr masterMeshAdapter( new AMP::Mesh::libMesh(meshMaster,"master") );
      AMP::Mesh::Mesh::shared_ptr slaveMeshAdapter( new AMP::Mesh::libMesh(meshMaster,"slave") );

      std::vector<AMP::Mesh::Mesh::shared_ptr> meshes;
      meshes.push_back( masterMeshAdapter );
      meshes.push_back( slaveMeshAdapter );
      AMP::Mesh::Mesh::shared_ptr manager( new AMP::Mesh::MultiMesh( globalComm, meshes ) );

      AMP::Discretization::DOFManager::shared_ptr NodalVectorDOF = 
        AMP::Discretization::simpleDOFManager::create(manager,AMP::Mesh::Vertex,1,3,true);

      std::vector<double> offset(3,0.0);
      offset[0] = input_db->getDouble("xOffset");
      offset[1] = input_db->getDouble("yOffset");
      offset[2] = input_db->getDouble("zOffset");
      slaveMeshAdapter->displaceMesh( offset );

      std::vector<AMP::Mesh::MeshElementID> masterNodes;
      std::vector<AMP::Mesh::MeshElementID> slaveNodes;

      getMasterAndSlaveNodes(masterNodes, slaveNodes, masterMeshAdapter, slaveMeshAdapter);

      boost::shared_ptr<AMP::Operator::ElementPhysicsModel> masterElementPhysicsModel;
      boost::shared_ptr<AMP::Operator::LinearBVPOperator> masterOperator = boost::dynamic_pointer_cast<
      AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(masterMeshAdapter,
										       "BVPOperator",
										       input_db,
										       masterElementPhysicsModel));

      boost::shared_ptr<AMP::Operator::ElementPhysicsModel> slaveElementPhysicsModel;
      boost::shared_ptr<AMP::Operator::MechanicsLinearFEOperator> slaveOperator = boost::dynamic_pointer_cast<
      AMP::Operator::MechanicsLinearFEOperator>(AMP::Operator::OperatorBuilder::createOperator(slaveMeshAdapter,
											       "MechanicsLinearFEOperator",
											       input_db,
											       slaveElementPhysicsModel));

      boost::shared_ptr<AMP::Operator::ElementPhysicsModel> slaveElementPhysicsModel2;
      boost::shared_ptr<AMP::Operator::LinearBVPOperator> slaveOperator2 = boost::dynamic_pointer_cast<
      AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(slaveMeshAdapter,
										       "BVPOperator",
										       input_db,
										       slaveElementPhysicsModel2));

      boost::shared_ptr<AMP::Operator::OperatorParameters> dummyParams;
      boost::shared_ptr<AMP::Operator::ColumnOperator> columnOperator(new AMP::Operator::ColumnOperator(dummyParams));
      columnOperator->append(masterOperator);
      columnOperator->append(slaveOperator);

      boost::shared_ptr<AMP::Database> mlPreconditioner_db = input_db->getDatabase("ML_Solver"); 

      boost::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> masterPreconditionerParams(new 
          AMP::Solver::TrilinosMLSolverParameters(mlPreconditioner_db));
      masterPreconditionerParams->d_pOperator = masterOperator;
      boost::shared_ptr<AMP::Solver::TrilinosMLSolver> masterPreconditioner(new AMP::Solver::TrilinosMLSolver(masterPreconditionerParams));

      boost::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> slavePreconditionerParams(new 
          AMP::Solver::TrilinosMLSolverParameters(mlPreconditioner_db));
      slavePreconditionerParams->d_pOperator = slaveOperator2;
      boost::shared_ptr<AMP::Solver::TrilinosMLSolver> slavePreconditioner(new AMP::Solver::TrilinosMLSolver(slavePreconditionerParams));

      AMP::LinearAlgebra::Variable::shared_ptr masterVar = masterOperator->getOutputVariable();
      AMP::LinearAlgebra::Variable::shared_ptr slaveVar = slaveOperator->getOutputVariable();
      //AMP::LinearAlgebra::Variable::shared_ptr columnVar = columnOperator->getOutputVariable();

      boost::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
      boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> loadOperator1 =
        boost::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
									      AMP::Operator::OperatorBuilder::createOperator(slaveMeshAdapter, "LoadOperator", input_db, dummyModel));
      loadOperator1->setVariable(slaveVar);

      AMP::LinearAlgebra::Variable::shared_ptr displacement( new AMP::LinearAlgebra::Variable("displacement") );

      AMP::LinearAlgebra::Vector::shared_ptr nullVec;
      AMP::LinearAlgebra::Vector::shared_ptr columnSolVec = AMP::LinearAlgebra::createVector(NodalVectorDOF,displacement);
      AMP::LinearAlgebra::Vector::shared_ptr columnRhsVec = columnSolVec->cloneVector();
      AMP::LinearAlgebra::Vector::shared_ptr columnResVec = columnSolVec->cloneVector();

      AMP::LinearAlgebra::VS_Mesh masterMeshSelector( masterMeshAdapter );
      AMP::LinearAlgebra::VS_Mesh slaveMeshSelector(  slaveMeshAdapter  );

      AMP::LinearAlgebra::Vector::shared_ptr masterSolVec = columnSolVec->select( masterMeshSelector, displacement->getName() );
      AMP::LinearAlgebra::Vector::shared_ptr slaveSolVec  = columnSolVec->select( slaveMeshSelector, displacement->getName()  );

      columnSolVec->zero();
      columnRhsVec->zero();
      columnResVec->zero();

      loadOperator1->apply(nullVec, nullVec, columnRhsVec, 1.0, 0.0);

      int maxIters = input_db->getInteger("maxIters");

      {
        std::cout<<"MPC-MLCG: "<<std::endl;
        AMP::LinearAlgebra::Vector::shared_ptr MatOutVec = columnSolVec->cloneVector();
        AMP::LinearAlgebra::Vector::shared_ptr solOldVec = columnSolVec->cloneVector();
        AMP::LinearAlgebra::Vector::shared_ptr resOldVec = columnSolVec->cloneVector();
        AMP::LinearAlgebra::Vector::shared_ptr pOldVec = columnSolVec->cloneVector();
        AMP::LinearAlgebra::Vector::shared_ptr pNewVec = columnSolVec->cloneVector();
        AMP::LinearAlgebra::Vector::shared_ptr zVec = columnSolVec->cloneVector();

        columnOperator->apply(nullVec, columnSolVec, MatOutVec, 1.0, 0.0);
        addSlaveToMaster(slaveVar, masterVar, MatOutVec, slaveMeshAdapter, masterMeshAdapter, slaveNodes, masterNodes);
        setSlaveToZero(slaveVar, MatOutVec, slaveMeshAdapter, slaveNodes);

        columnResVec->subtract(columnRhsVec, MatOutVec);

        AMP::LinearAlgebra::Vector::shared_ptr zMasterVec = masterOperator->subsetOutputVector( zVec );
        AMP::LinearAlgebra::Vector::shared_ptr zSlaveVec = slaveOperator->subsetOutputVector( zVec );
        AMP::LinearAlgebra::Vector::shared_ptr resMasterVec = masterOperator->subsetOutputVector( columnResVec );
        AMP::LinearAlgebra::Vector::shared_ptr resSlaveVec = slaveOperator->subsetOutputVector( columnResVec );
        masterPreconditioner->solve(resMasterVec, zMasterVec);
        slavePreconditioner->solve(resSlaveVec, zSlaveVec);

        pOldVec->copyVector(zVec);

        for(int iter = 0; iter < maxIters; iter++) {
          double resNorm = columnResVec->L2Norm();
          double solMaxNorm = columnSolVec->maxNorm();

          std::cout<<"Iter = "<<iter<<" ResNorm2 = "<<std::setprecision(15)<<resNorm<<
            " SolMaxNorm = "<<std::setprecision(15)<<solMaxNorm<<std::endl;

          copyMasterToSlave(slaveVar, masterVar, pOldVec, slaveMeshAdapter, masterMeshAdapter, slaveNodes, masterNodes);

          resOldVec->copyVector(columnResVec);
          solOldVec->copyVector(columnSolVec);

          double resOldDotZ = resOldVec->dot(zVec);

          columnOperator->apply(nullVec, pOldVec, MatOutVec, 1.0, 0.0);
          addSlaveToMaster(slaveVar, masterVar, MatOutVec, slaveMeshAdapter, masterMeshAdapter, slaveNodes, masterNodes);
          setSlaveToZero(slaveVar, MatOutVec, slaveMeshAdapter, slaveNodes);

          double alphaDenom = MatOutVec->dot(pOldVec);

          double alpha = resOldDotZ/alphaDenom;

          columnSolVec->axpy(alpha, pOldVec, solOldVec);

          columnResVec->axpy(-alpha, MatOutVec, resOldVec);

          zMasterVec = masterOperator->subsetOutputVector( zVec );
          zSlaveVec = slaveOperator->subsetOutputVector( zVec );
          resMasterVec = masterOperator->subsetOutputVector( columnResVec );
          resSlaveVec = slaveOperator->subsetOutputVector( columnResVec );
          masterPreconditioner->solve(resMasterVec, zMasterVec);
          slavePreconditioner->solve(resSlaveVec, zSlaveVec);

          double resNewDotZ = columnResVec->dot(zVec);

          double beta = resNewDotZ/resOldDotZ;

          std::cout<<"Iter = "<<iter<<" alpha = "<<std::setprecision(15)<<alpha<<
            " beta = "<<std::setprecision(15)<<beta<<std::endl<<std::endl;

          pNewVec->axpy(beta, pOldVec, zVec);

          pOldVec->copyVector(pNewVec);
        }
      }

      printSolution(masterMeshAdapter, masterSolVec, (exeName + "-master"));
      printSolution(slaveMeshAdapter, slaveSolVec, (exeName + "-slave"));

   }

   ut->passes(exeName);
}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::vector<std::string> exeNames;

  if(argc == 1) {
    exeNames.push_back("testMLCG_MPC-1");
    exeNames.push_back("testMLCG_MPC-2");
    exeNames.push_back("testMLCG_MPC-3");
    exeNames.push_back("testMLCG_MPC-4");
    exeNames.push_back("testMLCG_MPC-5");
  } else {
    for(int i = 1; i < argc; i++) {
      char inpName[100];
      sprintf(inpName, "testMLCG_MPC-%s", argv[i]);
      exeNames.push_back(inpName);
    }//end for i
  }

  for(size_t i = 0; i < exeNames.size(); i++) {
    try {
      myTest(&ut, exeNames[i]);
    } catch (std::exception &err) {
      std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
      ut.failure("ERROR: While testing");
    } catch( ... ) {
      std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
      ut.failure("ERROR: While testing");
    }
  }

  ut.report();

  int num_failed = ut.NumFailGlobal();
  AMP::AMPManager::shutdown();
  return num_failed;
}  

