#include <string>
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "ampmesh/Mesh.h"
#include "boost/shared_ptr.hpp"
#include "utils/InputDatabase.h"
#include "utils/Utilities.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/Database.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/SimpleVector.h"
#include "operators/subchannel/SubchannelTwoEqNonlinearOperator.h"
#include "operators/subchannel/SubchannelTwoEqLinearOperator.h"
#include "operators/OperatorBuilder.h"
#include "solvers/ColumnSolver.h"
#include "solvers/PetscKrylovSolverParameters.h"
#include "solvers/PetscKrylovSolver.h"
#include "solvers/PetscSNESSolverParameters.h"
#include "solvers/PetscSNESSolver.h"
#include "solvers/TrilinosMLSolver.h"

#include "ampmesh/SiloIO.h"
#include "vectors/VectorBuilder.h"
#include "discretization/simpleDOF_Manager.h"
#include "ampmesh/StructuredMeshHelper.h"

void flowTest(AMP::UnitTest *ut, std::string exeName )
{
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;
    std::string silo_name = exeName;
  AMP::PIO::logAllNodes(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  // Read the input file
  boost::shared_ptr<AMP::InputDatabase>  input_db ( new AMP::InputDatabase ( "input_db" ) );
  AMP::InputManager::getManager()->parseInputFile ( input_file , input_db );

  // Get the Mesh database and create the mesh parameters
  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  boost::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
  boost::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(mesh_db));
  meshParams->setComm(globalComm);

  // Create the meshes from the input database
  boost::shared_ptr<AMP::Mesh::Mesh> subchannelMesh = AMP::Mesh::Mesh::buildMesh(meshParams);
  AMP::Mesh::Mesh::shared_ptr xyFaceMesh;
  xyFaceMesh = subchannelMesh->Subset( AMP::Mesh::StructuredMeshHelper::getXYFaceIterator( subchannelMesh , 0 ) );

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;

  // get subchannel physics model
  boost::shared_ptr<AMP::Database> subchannelPhysics_db = input_db->getDatabase("SubchannelPhysicsModel");
  boost::shared_ptr<AMP::Operator::ElementPhysicsModelParameters> params( new AMP::Operator::ElementPhysicsModelParameters(subchannelPhysics_db));
  boost::shared_ptr<AMP::Operator::SubchannelPhysicsModel>  subchannelPhysicsModel (new AMP::Operator::SubchannelPhysicsModel(params));

  // get nonlinear operator database
  boost::shared_ptr<AMP::Database> nonlinearOperator_db = input_db->getDatabase("SubchannelTwoEqNonlinearOperator");
  // create parameters
  boost::shared_ptr<AMP::Operator::SubchannelOperatorParameters> subchannelOpParams(new AMP::Operator::SubchannelOperatorParameters( nonlinearOperator_db ));
  // put mesh into parameters
  subchannelOpParams->d_Mesh = xyFaceMesh ;
  // put subchannel physics model into parameters
  subchannelOpParams->d_subchannelPhysicsModel = subchannelPhysicsModel;

  // create nonlinear operator
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel;
  boost::shared_ptr<AMP::Operator::SubchannelTwoEqNonlinearOperator> nonlinearOperator =
      boost::dynamic_pointer_cast<AMP::Operator::SubchannelTwoEqNonlinearOperator>(AMP::Operator::OperatorBuilder::createOperator(
      subchannelMesh ,"SubchannelTwoEqNonlinearOperator",input_db,elementModel ));
  

  // create linear operator
  boost::shared_ptr<AMP::Operator::SubchannelTwoEqLinearOperator> linearOperator =
      boost::dynamic_pointer_cast<AMP::Operator::SubchannelTwoEqLinearOperator>(AMP::Operator::OperatorBuilder::createOperator(
      subchannelMesh ,"SubchannelTwoEqLinearOperator",input_db,elementModel ));

  
  // pass creation test
  ut->passes(exeName+": creation");
  std::cout.flush();

  // get input and output variables
  AMP::LinearAlgebra::Variable::shared_ptr inputVariable  = nonlinearOperator->getInputVariable();
  AMP::LinearAlgebra::Variable::shared_ptr outputVariable = nonlinearOperator->getOutputVariable();

  int DofsPerFace =  2;
  AMP::Discretization::DOFManager::shared_ptr faceDOFManager = AMP::Discretization::simpleDOFManager::create( subchannelMesh, 
      AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,1), AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,0), DofsPerFace );
//  AMP::Discretization::DOFManager::shared_ptr faceDOFManager = AMP::Discretization::simpleDOFManager::create( xyFaceMesh, AMP::Mesh::Face, 1, DofsPerFace, true);

  // create solution, rhs, and residual vectors
  AMP::LinearAlgebra::Vector::shared_ptr manufacturedVec = AMP::LinearAlgebra::createVector( faceDOFManager , inputVariable  , true );
  AMP::LinearAlgebra::Vector::shared_ptr solVec = AMP::LinearAlgebra::createVector( faceDOFManager , inputVariable  , true );
  AMP::LinearAlgebra::Vector::shared_ptr rhsVec = AMP::LinearAlgebra::createVector( faceDOFManager , outputVariable  , true );
  AMP::LinearAlgebra::Vector::shared_ptr resVec = AMP::LinearAlgebra::createVector( faceDOFManager , outputVariable  , true );

  // get exit pressure
  double Pout = nonlinearOperator_db->getDouble("Exit_Pressure");
  // calculate initial guess
  // get inlet temperature
  double Tin = nonlinearOperator_db->getDouble("Inlet_Temperature");
  // compute inlet enthalpy
  std::map<std::string, boost::shared_ptr<std::vector<double> > > enthalpyArgMap;
  enthalpyArgMap.insert(std::make_pair("temperature",new std::vector<double>(1,Tin)));
  enthalpyArgMap.insert(std::make_pair("pressure",   new std::vector<double>(1,Pout)));
  std::vector<double> enthalpyResult(1);
  subchannelPhysicsModel->getProperty("Enthalpy",enthalpyResult,enthalpyArgMap); 
  double hin = enthalpyResult[0];
  std::cout<< "Enthalpy Solution:"<< hin <<std::endl;

  // set manufactured solution
  // get the Iterators for the subchannel mesh
  AMP::Mesh::MeshIterator face     = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
  AMP::Mesh::MeshIterator end_face = face.end();

  std::vector<size_t> dofs;

  //initial guess
  for( ; face != end_face; ++face){
    faceDOFManager->getDOFs( face->globalID(), dofs );
    solVec->setValueByGlobalID(dofs[0], 1000);
    solVec->setValueByGlobalID(dofs[1], Pout);
  }

  // get nonlinear solver database
  boost::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase("NonlinearSolver"); 
  
  // get linear solver database
  boost::shared_ptr<AMP::Database> linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver"); 
 

  // put manufactured RHS into resVec
  nonlinearOperator->reset(subchannelOpParams);
//  nonlinearOperator->apply(rhsVec, manufacturedVec, resVec, 1.0, 0.0);
  
  linearOperator->reset(nonlinearOperator->getJacobianParameters(solVec));
  linearOperator->apply(rhsVec, solVec, resVec, 1.0, -1.0);

  
  // create nonlinear solver parameters
  boost::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(new AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db));


  // change the next line to get the correct communicator out
  nonlinearSolverParams->d_comm = globalComm;
  nonlinearSolverParams->d_pOperator = nonlinearOperator;
  nonlinearSolverParams->d_pInitialGuess = solVec;


  // create nonlinear solver
  boost::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));


  // create linear solver
  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver = nonlinearSolver->getKrylovSolver();

  // create preconditioner
  boost::shared_ptr<AMP::Database> Preconditioner_db =  linearSolver_db->getDatabase("Preconditioner");
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> PreconditionerParams(new AMP::Solver::SolverStrategyParameters(Preconditioner_db));
  PreconditionerParams->d_pOperator = linearOperator;
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> linearFlowPreconditioner(new AMP::Solver::TrilinosMLSolver(PreconditionerParams));
  // set preconditioner
  linearSolver->setPreconditioner(linearFlowPreconditioner);


  // don't use zero initial guess
  nonlinearSolver->setZeroInitialGuess(false);

  // solve
  nonlinearSolver->solve(rhsVec, solVec);

  face  = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
  faceDOFManager->getDOFs( face->globalID(), dofs );
  std::cout<< "Inlet Computed Enthalpy = "<<solVec->getValueByGlobalID(dofs[0]) << " Computed Pressure = "<<solVec->getValueByGlobalID(dofs[1]) ;
  // compute inlet temperature 
  std::map<std::string, boost::shared_ptr<std::vector<double> > > temperatureArgMap;
  temperatureArgMap.insert(std::make_pair("enthalpy",new std::vector<double>(1,solVec->getValueByGlobalID(dofs[0]))));
  temperatureArgMap.insert(std::make_pair("pressure",new std::vector<double>(1,solVec->getValueByGlobalID(dofs[1]))));
  std::vector<double> temperatureResult(1);
  subchannelPhysicsModel->getProperty("Temperature", temperatureResult, temperatureArgMap); 
  double compTin = temperatureResult[0];
  std::cout<< "Temperature Inlet :"<< compTin <<std::endl;

  face = --end_face;
  faceDOFManager->getDOFs( face->globalID(), dofs );
  std::cout<< "Outlet Computed Enthalpy = "<<solVec->getValueByGlobalID(dofs[0]) << " Computed Pressure = "<<solVec->getValueByGlobalID(dofs[1]) ;
  std::map<std::string, boost::shared_ptr<std::vector<double> > > outTemperatureArgMap;
  outTemperatureArgMap.insert(std::make_pair("enthalpy",new std::vector<double>(1,solVec->getValueByGlobalID(dofs[0]))));
  outTemperatureArgMap.insert(std::make_pair("pressure",new std::vector<double>(1,solVec->getValueByGlobalID(dofs[1]))));
  std::vector<double> outTemperatureResult(1);
  subchannelPhysicsModel->getProperty("Temperature", outTemperatureResult, outTemperatureArgMap); 
  compTin = outTemperatureResult[0];
  std::cout<< "Temperature Outlet:"<< compTin <<std::endl;


  // print final solution
  int j=1;
  face  = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
  ++end_face;

  for( ; face != end_face; ++face,++j){
    faceDOFManager->getDOFs( face->globalID(), dofs );
     std::cout<< "Computed Enthalpy["<<j<<"] = "<<solVec->getValueByGlobalID(dofs[0]) << "Computed Pressure["<<j<<"] = "<<solVec->getValueByGlobalID(dofs[1]) ;
     std::cout<<" Manufactured Pressure["<<j<<"] = "<<manufacturedVec->getValueByGlobalID(dofs[1]) ;
     std::cout<<std::endl;
  }
  // print manufactured solution
  for( int i=0; i<5; i++) {
     std::cout<<std::endl;
  }

  // print absolute error between manufactured solution and final solution
  AMP::LinearAlgebra::Vector::shared_ptr absErrorVec = AMP::LinearAlgebra::createVector( faceDOFManager , inputVariable  , true );
  face  = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
  j=1;
  for( ; face != end_face; ++face,++j){
    faceDOFManager->getDOFs( face->globalID(), dofs );
    absErrorVec->setValueByGlobalID(dofs[1], manufacturedVec->getValueByGlobalID(dofs[1]) - solVec->getValueByGlobalID(dofs[1]));
    absErrorVec->setValueByGlobalID(dofs[0], manufacturedVec->getValueByGlobalID(dofs[0]) - solVec->getValueByGlobalID(dofs[0]));
    std::cout<<"Absolute_Error["<<j<<"] = "<<absErrorVec->getValueByGlobalID(dofs[1]);
    std::cout<<std::endl;
  }
  double absErrorNorm = absErrorVec->L2Norm();
  std::cout<<"L2 Norm of Absolute Error: "<<absErrorNorm<<std::endl;
  
  // print relative error between manufactured solution and final solution
  AMP::LinearAlgebra::Vector::shared_ptr relErrorVec = AMP::LinearAlgebra::createVector( faceDOFManager , inputVariable  , true );
  face  = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
  j=1;
  for( ; face != end_face; ++face,++j){
    faceDOFManager->getDOFs( face->globalID(), dofs );
     double relError;
     if (std::abs(manufacturedVec->getValueByGlobalID(dofs[1])) < 1.0e-15){
        relError = 0.0;
     } else {
        relError = absErrorVec->getValueByGlobalID(dofs[1])/manufacturedVec->getValueByGlobalID(dofs[1]);
     }
     relErrorVec->setValueByGlobalID(dofs[1], relError);
     std::cout<<"Relative_Error["<<j<<"] = "<<relError;
     std::cout<<std::endl;
  }
  double relErrorNorm = relErrorVec->L2Norm();
  std::cout<<"L2 Norm of Relative Error: "<<relErrorNorm<<std::endl;

  // check that norm of relative error is less than tolerance
  if(relErrorNorm > 1.0e-7){
     ut->failure(exeName+": manufactured solution test");
  } else {
     ut->passes(exeName+": manufactured solution test");
  }

  input_db.reset();

#ifdef USES_SILO
    // Register the quantities to plot
    AMP::Mesh::SiloIO::shared_ptr  siloWriter( new AMP::Mesh::SiloIO );
    siloWriter->registerVector( manufacturedVec, xyFaceMesh, AMP::Mesh::Face, "ManufacturedSolution" );
    siloWriter->registerVector( solVec, xyFaceMesh, AMP::Mesh::Face, "ComputedSolution" );
    siloWriter->writeFile( silo_name , 0 );
#endif


}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    try {
        flowTest(&ut, "testSubchannelSolution");
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

