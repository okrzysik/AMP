#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>

#include "boost/shared_ptr.hpp"

#include "operators/VolumeIntegralOperator.h"
#include "operators/NeutronicsRhs.h"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"
#include "materials/Material.h"

#include "ampmesh/Mesh.h"
#include "vectors/VectorBuilder.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"

#include "utils/Writer.h"


#include "operators/mechanics/MechanicsLinearFEOperator.h"
#include "operators/mechanics/MechanicsNonlinearFEOperator.h"

#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"

#include "operators/boundary/DirichletVectorCorrection.h"

#include "operators/diffusion/FickSoretNonlinearFEOperator.h"
#include "operators/BVPOperatorParameters.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/ColumnOperator.h"
#include "operators/OperatorBuilder.h"

#include "solvers/ColumnSolver.h"
#include "solvers/petsc/PetscKrylovSolverParameters.h"
#include "solvers/petsc/PetscKrylovSolver.h"
#include "solvers/petsc/PetscSNESSolverParameters.h"
#include "solvers/petsc/PetscSNESSolver.h"

#include "solvers/trilinos/TrilinosMLSolver.h"


struct null_deleter {
    void operator()(void const *) const {}
};


void fickSoretTest(AMP::UnitTest *ut, std::string exeName, std::vector<double> &results)
{
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

//--------------------------------------------------
//   Create the Mesh.
//--------------------------------------------------
    AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
    boost::shared_ptr<AMP::Database>  mesh_db = input_db->getDatabase("Mesh");
    boost::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(new AMP::Mesh::MeshParameters(mesh_db));
    mgrParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
    boost::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh(mgrParams);
//--------------------------------------------------

//--------------------------------------------------
// Create a DOF manager for a nodal vector 
//--------------------------------------------------
  int DOFsPerNode = 1;
  int nodalGhostWidth = 1;
  bool split = true;
  AMP::Discretization::DOFManager::shared_ptr nodalDofMap      = AMP::Discretization::simpleDOFManager::create(meshAdapter, AMP::Mesh::Vertex, nodalGhostWidth,      DOFsPerNode,    split);
//--------------------------------------------------

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create a nonlinear BVP operator for nonlinear Fick-Soret diffusion
  AMP_INSIST( input_db->keyExists("testNonlinearFickSoretBVPOperator"), "key missing!" );

  // Create nonlinear FickSoret BVP operator and access volume nonlinear FickSoret operator
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
  boost::shared_ptr<AMP::Operator::Operator> nlinBVPOperator =
    AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
                           "testNonlinearFickSoretBVPOperator",
                           input_db,
                           elementPhysicsModel);
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nlinBVPOp =
          boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(nlinBVPOperator);
  boost::shared_ptr<AMP::Operator::FickSoretNonlinearFEOperator> nlinOp =
         boost::dynamic_pointer_cast<AMP::Operator::FickSoretNonlinearFEOperator>(nlinBVPOp->getVolumeOperator());
  boost::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> fickOp =
         boost::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(nlinOp->getFickOperator());
  boost::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> soretOp =
         boost::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(nlinOp->getSoretOperator());

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // use the linear BVP operator to create a Fick linear operator with bc's
  AMP_INSIST( input_db->keyExists("testLinearFickBVPOperator"), "key missing!" );

  boost::shared_ptr<AMP::Operator::Operator> linBVPOperator =
    AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
                           "testLinearFickBVPOperator",
                           input_db,
                           elementPhysicsModel);
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linBVPOp =
          boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(linBVPOperator);
  //boost::shared_ptr<AMP::Operator::DiffusionLinearFEOperator> linOp =
 //         boost::dynamic_pointer_cast<AMP::Operator::DiffusionLinearFEOperator>(linBVPOp->getVolumeOperator());

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // Set up input and output variables
  //AMP::LinearAlgebra::Variable::shared_ptr tVar(fickOp->getInputVariable(AMP::Operator::Diffusion::TEMPERATURE));
  //AMP::LinearAlgebra::Variable::shared_ptr cVar(fickOp->getInputVariable(AMP::Operator::Diffusion::CONCENTRATION));
  AMP::LinearAlgebra::Variable::shared_ptr tVar(new AMP::LinearAlgebra::Variable("temp") );
  AMP::LinearAlgebra::Variable::shared_ptr cVar(fickOp->getOutputVariable() );
  boost::shared_ptr<AMP::LinearAlgebra::Variable> fsOutVar(nlinBVPOp->getOutputVariable());

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create solution, rhs, and residual vectors

  AMP::LinearAlgebra::Vector::shared_ptr solVec = AMP::LinearAlgebra::createVector( nodalDofMap, cVar );
  AMP::LinearAlgebra::Vector::shared_ptr rhsVec = AMP::LinearAlgebra::createVector( nodalDofMap, fsOutVar );
  AMP::LinearAlgebra::Vector::shared_ptr resVec = AMP::LinearAlgebra::createVector( nodalDofMap, fsOutVar );

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create parameters for reset test and reset fick and soret operators

  AMP::LinearAlgebra::Vector::shared_ptr tVec = AMP::LinearAlgebra::createVector( nodalDofMap, tVar );
  
  fickOp->setVector(0, tVec);
  soretOp->setVector(0, tVec);

  std::vector<AMP::LinearAlgebra::Vector::shared_ptr> fickFrozen  = fickOp-> getFrozen();
  std::vector<AMP::LinearAlgebra::Vector::shared_ptr> soretFrozen = soretOp->getFrozen();

  double lenscale = input_db->getDouble("LengthScale");
  soretFrozen[AMP::Operator::Diffusion::TEMPERATURE]->setToScalar(300.);  // Fill in manufactured solution
  int zeroGhostWidth = 0;
  AMP::Mesh::MeshIterator  iterator = meshAdapter->getIterator(AMP::Mesh::Vertex, zeroGhostWidth);
  for( ; iterator != iterator.end(); iterator++ ) {
    double x, y;
    std::valarray<double> poly(10);
  x = ( iterator->coord() )[0];
  y = ( iterator->coord() )[1];
  std::vector<size_t> gid;
  nodalDofMap->getDOFs ( iterator->globalID() , gid);
    double value = 300. + 450*(1.-(x*x/lenscale/lenscale+y*y/lenscale/lenscale));
    fickFrozen[AMP::Operator::Diffusion::TEMPERATURE]->setValueByGlobalID(gid[0], value);
    soretFrozen[AMP::Operator::Diffusion::TEMPERATURE]->setValueByGlobalID(gid[0], value);
  }

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  //Initial guess

  double initialValue = input_db->getDouble("InitialValue");
  solVec->setToScalar(initialValue);
  double initialGuessNorm  = solVec->L2Norm();
  std::cout << "initial guess norm = " << initialGuessNorm <<"\n";

  nlinBVPOp->modifyInitialSolutionVector(solVec);

  initialGuessNorm  = solVec->L2Norm();
  std::cout << "initial guess norm  after apply = " << initialGuessNorm <<"\n";

  rhsVec->setToScalar(0.0);
  nlinBVPOp->modifyRHSvector(rhsVec);

  //----------------------------------------------------------------------------------------------------------------------------------------------/

  boost::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase("NonlinearSolver"); 
  boost::shared_ptr<AMP::Database> linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver"); 

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // initialize the nonlinear solver
  boost::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(new
      AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db));

  // change the next line to get the correct communicator out
  nonlinearSolverParams->d_comm = globalComm;
  nonlinearSolverParams->d_pOperator = nlinBVPOp;
  nonlinearSolverParams->d_pInitialGuess = solVec;

  boost::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  boost::shared_ptr<AMP::Database> fickPreconditioner_db = linearSolver_db->getDatabase("Preconditioner");
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> fickPreconditionerParams(new AMP::Solver::SolverStrategyParameters(fickPreconditioner_db));
  fickPreconditionerParams->d_pOperator = linBVPOp;
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> linearFickPreconditioner(new AMP::Solver::TrilinosMLSolver(fickPreconditionerParams));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // register the preconditioner with the Jacobian free Krylov solver
  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver = nonlinearSolver->getKrylovSolver();

  linearSolver->setPreconditioner(linearFickPreconditioner);

  nlinBVPOp->apply(rhsVec, solVec, resVec, 1.0, -1.0);
  double initialResidualNorm  = resVec->L2Norm();

  AMP::pout<<"Initial Residual Norm: "<<initialResidualNorm<<std::endl;

  nonlinearSolver->setZeroInitialGuess(false);

  nonlinearSolver->solve(rhsVec, solVec);

  nlinBVPOp->apply(rhsVec, solVec, resVec, 1.0, -1.0);

  double finalResidualNorm  = resVec->L2Norm();

  std::cout<<"Final Residual Norm: "<<finalResidualNorm<<std::endl;

  solVec->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
  resVec->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // evaluate and register material coefficients for graphical output

  AMP::LinearAlgebra::Variable::shared_ptr fickCoeffVar(new AMP::LinearAlgebra::Variable("FickCoefficient"));
  AMP::LinearAlgebra::Variable::shared_ptr soretCoeffVar(new AMP::LinearAlgebra::Variable("SoretCoefficient"));
  AMP::LinearAlgebra::Vector::shared_ptr fickCoeffVec  = AMP::LinearAlgebra::createVector( nodalDofMap, fickCoeffVar );
  AMP::LinearAlgebra::Vector::shared_ptr soretCoeffVec = AMP::LinearAlgebra::createVector( nodalDofMap, soretCoeffVar );
  boost::shared_ptr<AMP::Operator::DiffusionTransportModel> fickModel = fickOp->getTransportModel();
  boost::shared_ptr<AMP::Operator::DiffusionTransportModel> soretModel = soretOp->getTransportModel();

  {
      int zeroGhostWidth = 0;
      AMP::Mesh::MeshIterator  iterator = meshAdapter->getIterator(AMP::Mesh::Vertex, zeroGhostWidth);
      size_t nnodes = fickCoeffVec->getLocalSize(), node;
      std::vector<size_t> gids(nnodes);
      std::vector<double> temp(nnodes), conc(nnodes), fickCoeff(nnodes), soretCoeff(nnodes), burn(nnodes);
      for(node=0 ; iterator != iterator.end(); iterator++ ) {
      std::vector<size_t> gid;
      nodalDofMap->getDOFs ( iterator->globalID() , gid);
          gids[node] = gid[0];
          node++;
      }
      AMP_INSIST(node==nnodes, "invalid count");
      fickFrozen[AMP::Operator::Diffusion::TEMPERATURE]->getValuesByGlobalID(nnodes, &gids[0], &temp[0]);
      solVec->getValuesByGlobalID(nnodes, &gids[0], &conc[0]);
      // this is  used to plot the fick and soret coefficnets used.  commenting it out till someone finds out.
      // This is kevin - i found out because the vector wasn't used when silo is not enabled.
      std::map<std::string,  boost::shared_ptr<std::vector<double> > > args;
      args.insert(std::make_pair("temperature",   boost::shared_ptr<std::vector<double> >(&temp,null_deleter())));
      args.insert(std::make_pair("concentration", boost::shared_ptr<std::vector<double> >(&conc,null_deleter())));
      args.insert(std::make_pair("burnup",        boost::shared_ptr<std::vector<double> >(&burn,null_deleter())));
      fickModel->getTransport(fickCoeff, args); 
      soretModel->getTransport(soretCoeff, args ); 
      fickCoeffVec->setValuesByGlobalID(nnodes, &gids[0], &fickCoeff[0]);
      soretCoeffVec->setValuesByGlobalID(nnodes, &gids[0], &soretCoeff[0]);
  }

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // write graphical output

#ifdef USE_EXT_SILO
     AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter("Silo");
     siloWriter->registerMesh( meshAdapter );

     siloWriter->registerVector( solVec, meshAdapter, AMP::Mesh::Vertex, "Solution" );
     siloWriter->registerVector( resVec, meshAdapter, AMP::Mesh::Vertex, "Residual" );
     siloWriter->registerVector( fickFrozen[AMP::Operator::Diffusion::TEMPERATURE], meshAdapter, AMP::Mesh::Vertex, "Temperature" );
     siloWriter->registerVector( fickCoeffVec, meshAdapter, AMP::Mesh::Vertex, "FickCoefficient" );
     siloWriter->registerVector( soretCoeffVec, meshAdapter, AMP::Mesh::Vertex, "ThermalDiffusionCoefficient" );
 
     siloWriter->writeFile( exeName, 0 );
#endif

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // store result
  {
      int zeroGhostWidth = 0;
      AMP::Mesh::MeshIterator  iterator = meshAdapter->getIterator(AMP::Mesh::Vertex, zeroGhostWidth);
      iterator = iterator.begin();
      size_t numNodes = 0;
      for(; iterator != iterator.end(); iterator++ ) numNodes++;
      results.resize(numNodes);

      iterator = iterator.begin();
      size_t iNode=0;
      for(; iterator != iterator.end(); iterator++ ) {
      std::vector<size_t> gid;
      nodalDofMap->getDOFs ( iterator->globalID() , gid);
        results[iNode] = solVec->getValueByGlobalID(gid[0] );
        iNode++;
      }
  }

  if(finalResidualNorm>1.0e-08) {
    ut->failure(exeName);
  } else {
    ut->passes("PetscSNES Solver successfully solves a nonlinear mechanics equation with Jacobian provided, FGMRES for Krylov");
  }
  ut->passes(exeName);
}


int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    try {
      std::vector<double>  results;
      fickSoretTest(&ut, "testPetscSNESSolver-NonlinearFickSoret-cylinder-OxMSRZC09-1", results);
      ut.passes("fick-soret quadratic external T");
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


