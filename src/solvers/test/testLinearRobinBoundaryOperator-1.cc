#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "ampmesh/SiloIO.h"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"
#include "materials/Material.h"

#include <string>
#include "boost/shared_ptr.hpp"
#include "utils/Utilities.h"
#include "vectors/Variable.h"

#include "vectors/Vector.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"

#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionLinearElement.h"
#include "operators/diffusion/DiffusionTransportModel.h"
#include "operators/ElementPhysicsModelFactory.h"
#include "operators/ElementOperationFactory.h"
#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/VolumeIntegralOperator.h"
#include "operators/boundary/DirichletMatrixCorrection.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/boundary/RobinMatrixCorrection.h"
#include "operators/boundary/NeumannVectorCorrection.h"
#include "operators/boundary/ColumnBoundaryOperator.h"

#include "../PetscKrylovSolverParameters.h"
#include "../PetscKrylovSolver.h"
#include "../PetscSNESSolverParameters.h"
#include "../PetscSNESSolver.h"

#include "../TrilinosMLSolver.h"


#define __PI__ 3.14159265
#define __INIT_T0__(x,y,z,sg)   (sg * cos(0.1 * __PI__ * x) * cos(0.1 * __PI__ * y) * cos(0.1 * __PI__ * z) )
#define __INIT_dTdx__(x,y,z,sg) (sg * -0.1 * __PI__ * sin(0.1 * __PI__ * x) * cos(0.1 * __PI__ * y) * cos(0.1 * __PI__ * z) )
#define __INIT_dTdy__(x,y,z,sg) (sg * -0.1 * __PI__ * cos(0.1 * __PI__ * x) * sin(0.1 * __PI__ * y) * cos(0.1 * __PI__ * z) )
#define __INIT_dTdz__(x,y,z,sg) (sg * -0.1 * __PI__ * cos(0.1 * __PI__ * x) * cos(0.1 * __PI__ * y) * sin(0.1 * __PI__ * z) )
#define __INIT_rhs__(x,y,z,sg)  (sg * -0.03 * __PI__ * __PI__ * cos(0.1 * __PI__ * x) * cos(0.1 * __PI__ * y) * cos(0.1 * __PI__ * z) )


void linearRobinTest(AMP::UnitTest *ut, std::string exeName )
{
  // Input and output file names
  //  #include <string>
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;
  ////////////////////////////////////
  //    INITIALIZE THE PROBLEM      //
  ////////////////////////////////////

  // Create the map to get an available material from a string.
  //  #include "materials/Material.h"

  // Construct a smart pointer to a new database.
  //  #include "boost/shared_ptr.hpp"
  //  #include "utils/InputDatabase.h"
  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));


  // Fill the database from the input file.
  //  #include "utils/InputManager.h"
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);


  // Print from all cores into the output files
  //   #include "utils/PIO.h"
  AMP::PIO::logAllNodes(log_file);

//--------------------------------------------------
//   Create the Mesh.
//--------------------------------------------------
  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  boost::shared_ptr<AMP::Database>  mesh_db = input_db->getDatabase("Mesh");
  boost::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(new AMP::Mesh::MeshParameters(mesh_db));
  mgrParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
  boost::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh(mgrParams);

//--------------------------------------------------
//   Create DOF Managers.
//--------------------------------------------------
  int DOFsPerNode = 1;
  int ghostWidth = 1;
  bool split = true;
  AMP::Discretization::DOFManager::shared_ptr nodalDofMap      = AMP::Discretization::simpleDOFManager::create(meshAdapter, AMP::Mesh::Vertex, ghostWidth, DOFsPerNode,    split);

    // Create a shared pointer to a Variable - Power - Output because it will be used in the "residual" location of apply. 

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  //------------------------------------------
  //   CREATE THE THERMAL BVP OPERATOR  //
  //------------------------------------------
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> diffusionOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
																						       "DiffusionBVPOperator",
																						       input_db,
																						       transportModel));


  AMP::LinearAlgebra::Vector::shared_ptr TemperatureInKelvinVec = AMP::LinearAlgebra::createVector( nodalDofMap, diffusionOperator->getOutputVariable(), split );
  AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec       = TemperatureInKelvinVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr ResidualVec            = TemperatureInKelvinVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr variableFluxVec        = TemperatureInKelvinVec->cloneVector();

  RightHandSideVec->zero();
  variableFluxVec->zero();
  double rhsNorm = RightHandSideVec->L2Norm();

  //------------------------------------------

  AMP::Operator::Operator::shared_ptr boundaryOp;
  boundaryOp = diffusionOperator->getBoundaryOperator(); 

  AMP::Operator::Operator::shared_ptr  robinBoundaryOp;
  robinBoundaryOp = (boost::dynamic_pointer_cast<AMP::Operator::ColumnBoundaryOperator>(boundaryOp) )->getBoundaryOperator(0);

  boost::shared_ptr<AMP::InputDatabase> robinboundaryDatabase = boost::dynamic_pointer_cast<AMP::InputDatabase>( input_db->getDatabase("RobinMatrixCorrection"));

  robinboundaryDatabase->putBool("constant_flux", false);
  robinboundaryDatabase->putBool("skip_matrix_correction", true);
  boost::shared_ptr<AMP::Operator::RobinMatrixCorrectionParameters> correctionParameters (new AMP::Operator::RobinMatrixCorrectionParameters( robinboundaryDatabase ) );
  //------------------------------------------


  //------------------------------------------
  // check the solution
  int zeroGhostWidth = 0;
  AMP::Mesh::MeshIterator  node = meshAdapter->getIterator(AMP::Mesh::Vertex, zeroGhostWidth);
  AMP::Mesh::MeshIterator  end_node = node.end();

  for( ; node != end_node ; ++node)
  {
    std::vector<size_t> gid;
    nodalDofMap->getDOFs ( node->globalID() , gid);

    double px = (node->coord())[0];
    double py = (node->coord())[1];
    double pz = (node->coord())[2];

    double val, rhs;

    rhs =  __INIT_rhs__(px, py, pz, -1.0); 
    RightHandSideVec->setValueByGlobalID(gid[0], rhs);

    if(fabs(pz - 1.0) <= 1.0e-12){
      val = __INIT_dTdz__(px, py, pz, 1.0);
      val = val +  __INIT_T0__(px, py, pz, 1.0);
      variableFluxVec->setValueByGlobalID(gid[0], val);
    }else if(fabs(pz + 1.0) <= 1.0e-12){
      val = __INIT_dTdz__(px, py, pz, -1.0);
      val = val +  __INIT_T0__(px, py, pz, 1.0);
      variableFluxVec->setValueByGlobalID(gid[0], val);
    }else if(fabs(px -1.0) <= 1.0e-12){
      val = __INIT_dTdx__(px, py, pz, 1.0);
      val = val +  __INIT_T0__(px, py, pz, 1.0);
      variableFluxVec->setValueByGlobalID(gid[0], val);
    }else if(fabs(px + 1.0) <= 1.0e-12){
      val = __INIT_dTdx__(px, py, pz, -1.0);
      val = val +  __INIT_T0__(px, py, pz, 1.0);
      variableFluxVec->setValueByGlobalID(gid[0], val);
    }else if(fabs(py - 1.0) <= 1.0e-12){
      val = __INIT_dTdy__(px, py, pz, 1.0);
      val = val +  __INIT_T0__(px, py, pz, 1.0);
      variableFluxVec->setValueByGlobalID(gid[0], val);
    }else if(fabs(py + 1.0) <= 1.0e-12){
      val = __INIT_dTdy__(px, py, pz, -1.0);
      val = val +  __INIT_T0__(px, py, pz, 1.0);
      variableFluxVec->setValueByGlobalID(gid[0], val);
    }
  }//end for node

  correctionParameters->d_variableFlux = variableFluxVec;
  robinBoundaryOp->reset(correctionParameters);

  //----------------------------------------------------------
  //  Integrate Nuclear Rhs over Desnity * Volume //
  //----------------------------------------------------------

  AMP_INSIST( input_db->keyExists("VolumeIntegralOperator"), "key missing!" );

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
  boost::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator = boost::dynamic_pointer_cast<
    AMP::Operator::VolumeIntegralOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
											  "VolumeIntegralOperator",
											  input_db,
											  stransportModel));

  // Create the power (heat source) vector.
  AMP::LinearAlgebra::Variable::shared_ptr SourceVar = sourceOperator->getOutputVariable();
  AMP::LinearAlgebra::Vector::shared_ptr   SourceVec = AMP::LinearAlgebra::createVector( nodalDofMap, SourceVar, split );
  SourceVec->zero();

  // convert the vector of specific power to power for a given basis.
  sourceOperator->apply(nullVec, RightHandSideVec, SourceVec, 1., 0.);

  //------------------------------------------
  //   Add the boundary conditions corrections //
  //------------------------------------------

  std::cout<<"RHS Norm before BC Correction "<<SourceVec->L2Norm()<<std::endl;

  diffusionOperator->modifyRHSvector(SourceVec) ;

  rhsNorm = SourceVec->L2Norm();
  std::cout<<"RHS Norm after BC Correction "<<rhsNorm<<std::endl;

  //------------------------------------------

  // make sure the database on theinput file exists for the linear solver
  AMP_INSIST(input_db->keyExists("LinearSolver"),   "Key ''LinearSolver'' is missing!");

  // Read the input file onto a database.
  boost::shared_ptr<AMP::Database>                 mlSolver_db   = input_db->getDatabase("LinearSolver"); 

  // Fill in the parameters fo the class with the info on the database.
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> mlSolverParams (new AMP::Solver::SolverStrategyParameters(mlSolver_db));

  // Define the operature to be used by the Solver.
  mlSolverParams->d_pOperator = diffusionOperator;

  // Set initial guess
  TemperatureInKelvinVec->setToScalar(1.0);

  // Check the initial L2 norm of the solution
  double initSolNorm = TemperatureInKelvinVec->L2Norm();
  std::cout<<"Initial Solution Norm: "<<initSolNorm<<std::endl;

  // Create the ML Solver
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver>         mlSolver(new AMP::Solver::TrilinosMLSolver(mlSolverParams));

  // Use a random initial guess?
  mlSolver->setZeroInitialGuess(false);

  // Solve the prblem.
  mlSolver->solve(SourceVec, TemperatureInKelvinVec);

  // Compute the residual
  diffusionOperator->apply(SourceVec, TemperatureInKelvinVec, ResidualVec);

  // Check the L2 norm of the final residual.
  double finalResidualNorm = ResidualVec->L2Norm();
  std::cout<<"Final Residual Norm: "<<finalResidualNorm<<std::endl;

  node = node.begin();
  AMP::LinearAlgebra::Vector::shared_ptr diffVec  = TemperatureInKelvinVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr exactVec = TemperatureInKelvinVec->cloneVector();

  diffVec->zero();
  exactVec->zero();

  for( ; node != end_node ; ++node)
  {
    std::vector<size_t> gid;
    nodalDofMap->getDOFs ( node->globalID() , gid);

    double px = (node->coord())[0];
    double py = (node->coord())[1];
    double pz = (node->coord())[2];

    double exact;
    exact =  __INIT_T0__(px, py, pz, 1.0); 
    exactVec->setValueByGlobalID(gid[0], exact);
  }

  diffVec->subtract(exactVec, TemperatureInKelvinVec);

  double exactNorm = exactVec->L1Norm();
  std::cout<<"L2norm of exactVec "<<exactNorm<<std::endl;

  double solutionNorm = TemperatureInKelvinVec->L1Norm();
  std::cout<<"L2norm of solutionVec "<<solutionNorm<<std::endl;

  double errorNorm = diffVec->L1Norm();
  std::cout<<"L1norm of DiffVec "<<errorNorm<<std::endl;

  if(errorNorm>1.0) {
    ut->failure("linear robin boundary operator verification test-1.");
  } else {
    ut->passes("linear robin boundary operator verification test-1.");
  }

  // Plot the results
  AMP::AMP_MPI globalComm = AMP::AMP_MPI(AMP_COMM_WORLD);

#ifdef USE_SILO
  AMP::Mesh::SiloIO::shared_ptr  siloWriter( new AMP::Mesh::SiloIO);

  siloWriter->registerVector( TemperatureInKelvinVec, meshAdapter, AMP::Mesh::Vertex, "TemperatureInKelvin" );
  siloWriter->registerVector( exactVec,               meshAdapter, AMP::Mesh::Vertex, "Exact"    );
 
  siloWriter->writeFile( input_file , 0 );
#endif

  ut->passes(exeName);

}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.push_back("testLinearRobinBoundaryOperator-1");
    // exeNames.push_back("testLinearRobinBoundaryOperator-2");

    for(unsigned int i = 0; i < exeNames.size(); i++) {
        try {
            linearRobinTest(&ut,exeNames[i]);
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


