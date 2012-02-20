
#include <iostream>

#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "mesh_communication.h"

#include "mpi.h"

#include "boost/shared_ptr.hpp"

#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>

#include "ampmesh/SiloIO.h"

#include "operators/OperatorBuilder.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/mechanics/MechanicsNonlinearFEOperator.h"
#include "operators/mechanics/VonMisesElastoPlasticModel.h"
#include "operators/boundary/DirichletVectorCorrection.h"

#include "solvers/PetscSNESSolverParameters.h"
#include "solvers/PetscSNESSolver.h"
#include "solvers/PetscKrylovSolverParameters.h"
#include "solvers/PetscKrylovSolver.h"
#include "solvers/TrilinosMLSolver.h"

#include "utils/ReadTestMesh.h"

void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::AMP_MPI globalComm = AMP::AMP_MPI(AMP_COMM_WORLD);
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  int numMeshes = input_db->getInteger("NumberOfMeshFiles");

  boost::shared_ptr<AMP::Database> ml_db = input_db->getDatabase("ML_Solver"); 

  int caseToRun = -1;
  if(input_db->keyExists("CaseToRun")) {
    caseToRun = input_db->getInteger("CaseToRun");
  }

  int caseCnt = 0;
  for(int meshId = 1; meshId <= numMeshes; meshId++) {
    AMP::Mesh::MeshManagerParameters::shared_ptr  meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
    AMP::Mesh::MeshManager::shared_ptr  manager ( new AMP::Mesh::MeshManager ( meshmgrParams ) );

    char meshFileKey[200];
    sprintf(meshFileKey, "mesh%d", meshId);

    std::string meshFile = input_db->getString(meshFileKey);

    const unsigned int mesh_dim = 3;
    boost::shared_ptr< ::Mesh > mesh(new ::Mesh(mesh_dim));

    if(globalComm.getRank() == 0) {
      AMP::readTestMesh(meshFile, mesh);
    }

    MeshCommunication().broadcast(*(mesh.get()));
    mesh->prepare_for_use(false);

    AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter ( new AMP::Mesh::MeshManager::Adapter (mesh) );

    manager->addMesh(meshAdapter, "mesh");

    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
    boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> loadOperator = boost::dynamic_pointer_cast<
      AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
            "LoadOperator", input_db, dummyModel));

    for(int useUL = 0; useUL < 2; useUL++) {
      std::string linOpDbName = "LinearBVP";
      if(useUL) {
        linOpDbName = linOpDbName + "_UL";
      } else {
        linOpDbName = linOpDbName + "_SS";
      }

      for(int useConsistent = 0; useConsistent < 2; useConsistent++) {
        std::string nlOpDbName = "NonlinearBVP";
        if(useUL) {
          nlOpDbName = nlOpDbName + "_UL";
        } else {
          nlOpDbName = nlOpDbName + "_SS";
        }
        if(useConsistent) {
          nlOpDbName = nlOpDbName + "_Consistent";
        } else {
          nlOpDbName = nlOpDbName + "_Continuum";
        }

        boost::shared_ptr<AMP::Operator::ElementPhysicsModel> materialModel;
        boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearBVPoperator = boost::dynamic_pointer_cast<
          AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
                nlOpDbName, input_db, materialModel));

        boost::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> nlVolOp = boost::dynamic_pointer_cast<
          AMP::Operator::MechanicsNonlinearFEOperator>(nonlinearBVPoperator->getVolumeOperator());

        boost::shared_ptr<AMP::Operator::VonMisesElastoPlasticModel> vonMisesModel = boost::dynamic_pointer_cast<
          AMP::Operator::VonMisesElastoPlasticModel>(nlVolOp->getMaterialModel());

        boost::shared_ptr<AMP::Operator::LinearBVPOperator> linearBVPoperator = boost::dynamic_pointer_cast<
          AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
                linOpDbName, input_db, materialModel));

        AMP::LinearAlgebra::Variable::shared_ptr dispVar = nonlinearBVPoperator->getOutputVariable();

        AMP::LinearAlgebra::Vector::shared_ptr nullVec;
        AMP::LinearAlgebra::Vector::shared_ptr solVec = manager->createVector(dispVar);
        AMP::LinearAlgebra::Vector::shared_ptr rhsVec = manager->createVector(dispVar);

        loadOperator->setVariable(dispVar);
        rhsVec->zero();
        loadOperator->apply(nullVec, nullVec, rhsVec, 1.0, 0.0);
        nonlinearBVPoperator->modifyRHSvector(rhsVec);

        for(int useJFNK = 0; useJFNK < 2; useJFNK++) {
          if(useConsistent) {
            if(useJFNK) {
              continue;
            }
          }
          for(int useEW = 0; useEW < 2; useEW++) {
            if( (caseToRun == -1) || (caseToRun == caseCnt) ) {
              std::string snesDbName = "SNES";
              if(useEW) {
                snesDbName = snesDbName + "_EW";
              }
              if(useJFNK) {
                snesDbName = snesDbName + "_JFNK";
              }
              boost::shared_ptr<AMP::Database> snes_db = input_db->getDatabase(snesDbName); 
              boost::shared_ptr<AMP::Database> ksp_db = snes_db->getDatabase("LinearSolver"); 

              solVec->zero();

              nonlinearBVPoperator->modifyInitialSolutionVector(solVec);

              if(globalComm.getRank() == 0) {
                std::cout<<"Solving ";
                if(useUL) {
                  std::cout<<"UL ";
                } else {
                  std::cout<<"SS ";
                }
                if(useConsistent) {
                  std::cout<<"Consistent ";
                } else {
                  std::cout<<"Continuum ";
                }
                if(useJFNK) {
                  std::cout<<"with JFNK ";
                }
                if(useEW) {
                  std::cout<<"with EW ";
                }
                std::cout<<" on mesh: "<<meshFile<<std::endl;
              }

              boost::shared_ptr<AMP::Solver::PetscSNESSolver> snesSolver;
              boost::shared_ptr<AMP::Solver::PetscKrylovSolver> kspSolver;
              boost::shared_ptr<AMP::Solver::TrilinosMLSolver> mlSolver;

              boost::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> mlParams(new AMP::Solver::TrilinosMLSolverParameters(ml_db));
              mlParams->d_pOperator = linearBVPoperator;
              mlSolver.reset(new AMP::Solver::TrilinosMLSolver(mlParams));

              boost::shared_ptr<AMP::Solver::PetscSNESSolverParameters> snesParams(new
                  AMP::Solver::PetscSNESSolverParameters(snes_db));
              snesParams->d_comm = globalComm;
              snesParams->d_pOperator = nonlinearBVPoperator;
              snesParams->d_pInitialGuess = solVec;
              if(!useJFNK) {
                boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> kspParams(new AMP::Solver::PetscKrylovSolverParameters(ksp_db));
                kspParams->d_pOperator = linearBVPoperator;
                kspParams->d_comm = globalComm;
                kspSolver.reset(new AMP::Solver::PetscKrylovSolver(kspParams));
                snesParams->d_pKrylovSolver = kspSolver;
              }
              snesSolver.reset(new AMP::Solver::PetscSNESSolver(snesParams));

              if(useJFNK) {
                kspSolver = snesSolver->getKrylovSolver();
              }

              kspSolver->setPreconditioner(mlSolver);

#ifdef USE_SILO
              manager->registerVectorAsData(solVec, "Displacement");
              char outFileName[200];
              sprintf(outFileName, "%s_case_%d", exeName.c_str(), caseCnt);
              manager->writeFile<AMP::Mesh::SiloIO>(outFileName, 1);
#endif

              snesSolver->solve(rhsVec, solVec);

              size_t numDofs = solVec->getGlobalSize();

              SNES snes = snesSolver->getSNESSolver();

              PetscInt nlFnCnt = 0;
              SNESGetNumberFunctionEvals(snes, &nlFnCnt);

              PetscInt snesItsCnt = 0;
              SNESGetIterationNumber(snes, &snesItsCnt);

              PetscInt kspItsCnt = 0;
              SNESGetLinearSolveIterations(snes, &kspItsCnt);

              SNESConvergedReason reason;
              SNESGetConvergedReason(snes, &reason);
              AMP_ASSERT( (reason == SNES_CONVERGED_FNORM_ABS) ||
                  (reason == SNES_CONVERGED_FNORM_RELATIVE) );

              unsigned int locPlCount = vonMisesModel->getLocalPlasticGaussPointCount();
              unsigned int locGpCount = vonMisesModel->getLocalGaussPointCount();

              unsigned int totPlCount = globalComm.sumReduce<unsigned int>(locPlCount);
              unsigned int totGpCount = globalComm.sumReduce<unsigned int>(locGpCount);

              double plFrac = 100.0*((double)totPlCount) / ((double)totGpCount);

              if(globalComm.getRank() == 0) {
                std::cout<<"Result: "<<numDofs<<" & "<<plFrac<<" & "<<nlFnCnt
                  <<" & "<<snesItsCnt<<" & "<<kspItsCnt<<" \\\\ "<<std::endl;
              }
            }//end if caseToRun
            caseCnt++;
          }//end useEW
        }//end useJFNK
      }//end useConsistent
    }//end useUL
  }//end for meshId

  ut->passes(exeName);
}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::string exeName = "jfnkPaperDriver";

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




