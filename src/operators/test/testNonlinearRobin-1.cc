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

#include "ampmesh/SiloIO.h"

#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "operators/BVPOperatorParameters.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/ColumnOperator.h"
#include "operators/OperatorBuilder.h"

#include "ampmesh/Mesh.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"


extern "C"{
#ifdef MPICH_SKIP_MPICXX
#define _FIX_FOR_PETSC_MPICH_CXX
#undef MPICH_SKIP_MPICXX
#endif

#ifdef OMPI_SKIP_MPICXX
#define _FIX_FOR_PETSC_OMPI_CXX
#undef OMPI_SKIP_MPICXX
#endif

#include "petsc.h"

#ifdef _FIX_FOR_PETSC_MPICH_CXX
#ifndef MPICH_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#endif
#endif

#ifdef _FIX_FOR_PETSC_OMPI_CXX
#ifndef OMPI_SKIP_MPICXX
#define OMPI_SKIP_MPICXX
#endif
#endif
}


void myTest(AMP::UnitTest *ut, std::string exeName)
{
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::AMP_MPI globalComm = AMP::AMP_MPI(AMP_COMM_WORLD);
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  // Get the Mesh database and create the mesh parameters
  boost::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
  boost::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(database));
  meshParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));

  // Create the meshes from the input database
  AMP::Mesh::Mesh::shared_ptr manager = AMP::Mesh::Mesh::buildMesh(meshParams);
  AMP::Mesh::Mesh::shared_ptr meshAdapter = manager->Subset( "cylinder" );

  AMP::pout<<"Constructing Nonlinear Thermal Operator..."<<std::endl;

  //-------------------------------------------------------------------------------------------//
  // create a nonlinear BVP operator for nonlinear thermal diffusion
  AMP_INSIST( input_db->keyExists("testNonlinearThermalOperator"), "key missing!" );

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearThermalOperator = boost::dynamic_pointer_cast<
    AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
											"testNonlinearThermalOperator",
											input_db,
											thermalTransportModel));

  //-------------------------------------------------------------------------------------------//
  // initialize the input variable
  boost::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> thermalVolumeOperator =
    boost::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(nonlinearThermalOperator->getVolumeOperator());

  boost::shared_ptr<AMP::LinearAlgebra::Variable> thermalVariable = thermalVolumeOperator->getOutputVariable();

  // create solution, rhs, and residual vectors
  AMP::Discretization::DOFManager::shared_ptr NodalScalarDOF = AMP::Discretization::simpleDOFManager::create(meshAdapter,AMP::Mesh::Vertex,1,1,true);
  AMP::LinearAlgebra::Vector::shared_ptr solVec = AMP::LinearAlgebra::createVector( NodalScalarDOF, thermalVariable, true );
  AMP::LinearAlgebra::Vector::shared_ptr rhsVec = AMP::LinearAlgebra::createVector( NodalScalarDOF, thermalVariable, true );
  AMP::LinearAlgebra::Vector::shared_ptr resVec = AMP::LinearAlgebra::createVector( NodalScalarDOF, thermalVariable, true );

  // create the following shared pointers for ease of use
  AMP::LinearAlgebra::Vector::shared_ptr nullVec;

#ifdef USES_SILO
  //-------------------------------------------------------------------------------------------//
  // Create the silo writer and register the data
  AMP::Mesh::SiloIO::shared_ptr  siloWriter( new AMP::Mesh::SiloIO);
  siloWriter->registerVector( solVec, meshAdapter, AMP::Mesh::Vertex, "Solution" );
  siloWriter->registerVector( resVec, meshAdapter, AMP::Mesh::Vertex, "Residual" );
#endif

  AMP::pout<<"Constructing Linear Thermal Operator..."<<std::endl;

  //-------------------------------------------------------------------------------------------//
  // now construct the linear BVP operator for thermal
  AMP_INSIST( input_db->keyExists("testLinearThermalOperator"), "key missing!" );
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
																	    AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
																							   "testLinearThermalOperator",
																							   input_db,
																							   thermalTransportModel));

  ////////////////////////////////////
  //  CREATE THE NEUTRONICS SOURCE  //
  ////////////////////////////////////
  AMP_INSIST(input_db->keyExists("NeutronicsOperator"), "Key ''NeutronicsOperator'' is missing!");
  boost::shared_ptr<AMP::Database>  neutronicsOp_db = input_db->getDatabase("NeutronicsOperator");
  boost::shared_ptr<AMP::Operator::NeutronicsRhsParameters> neutronicsParams(new AMP::Operator::NeutronicsRhsParameters( neutronicsOp_db ));
  neutronicsParams->d_Mesh = meshAdapter;
  boost::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOperator(new AMP::Operator::NeutronicsRhs( neutronicsParams ));

  // Create a DOF manager for a gauss point vector 
  int DOFsPerNode = 8;
  int ghostWidth = 1;
  bool split = true;
  AMP::Discretization::DOFManager::shared_ptr gauss_dof_map = AMP::Discretization::simpleDOFManager::create(meshAdapter, AMP::Mesh::Volume, ghostWidth, DOFsPerNode, split);

  AMP::LinearAlgebra::Variable::shared_ptr SpecificPowerVar = neutronicsOperator->getOutputVariable();
  AMP::LinearAlgebra::Vector::shared_ptr   SpecificPowerVec = AMP::LinearAlgebra::createVector( gauss_dof_map, SpecificPowerVar, split );

  neutronicsOperator->apply(nullVec, nullVec, SpecificPowerVec, 1., 0.);

  /////////////////////////////////////////////////////
  //  Integrate Nuclear Rhs over Desnity * Volume //
  /////////////////////////////////////////////////////

  AMP_INSIST( input_db->keyExists("VolumeIntegralOperator"), "key missing!" );

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
  boost::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator = boost::dynamic_pointer_cast<
    AMP::Operator::VolumeIntegralOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
											  "VolumeIntegralOperator",
											  input_db,
											  stransportModel));

  // Create the power (heat source) vector.
  AMP::LinearAlgebra::Variable::shared_ptr PowerInWattsVar = sourceOperator->getOutputVariable();
  AMP::LinearAlgebra::Vector::shared_ptr   PowerInWattsVec = AMP::LinearAlgebra::createVector( NodalScalarDOF, PowerInWattsVar, true );
  PowerInWattsVec->zero();

  // convert the vector of specific power to power for a given basis.
  sourceOperator->apply(nullVec, SpecificPowerVec, PowerInWattsVec, 1., 0.);

  rhsVec->copyVector(PowerInWattsVec);

  AMP::pout << "RHS L2 norm before corrections = " << (rhsVec->L2Norm()) <<"\n";
  AMP::pout << "RHS max before corrections = " << (rhsVec->max()) <<"\n";
  AMP::pout << "RHS min before corrections = " << (rhsVec->min()) <<"\n";

  nonlinearThermalOperator->modifyRHSvector(rhsVec);

  AMP::pout << "RHS L2 norm after corrections = " << (rhsVec->L2Norm()) <<"\n";
  AMP::pout << "RHS max after corrections = " << (rhsVec->max()) <<"\n";
  AMP::pout << "RHS min after corrections = " << (rhsVec->min()) <<"\n";

  //---------------------------------------------------------------------------------------------//
  //Initial guess

  double initGuess = input_db->getDoubleWithDefault("InitialGuess", 400.0);
  solVec->setToScalar(initGuess);

  AMP::pout << "initial guess L2 norm before corrections = " << (solVec->L2Norm()) <<"\n";
  AMP::pout << "initial guess max before corrections = " << (solVec->max()) <<"\n";
  AMP::pout << "initial guess min before corrections = " << (solVec->min()) <<"\n";

  nonlinearThermalOperator->modifyInitialSolutionVector(solVec);

  AMP::pout << "initial guess L2 norm after corrections = " << (solVec->L2Norm()) <<"\n";
  AMP::pout << "initial guess max after corrections = " << (solVec->max()) <<"\n";
  AMP::pout << "initial guess min after corrections = " << (solVec->min()) <<"\n";

  //---------------------------------------------------------------------------------------------/

  nonlinearThermalOperator->modifyInitialSolutionVector(solVec);
  linearThermalOperator->reset(nonlinearThermalOperator->getJacobianParameters(solVec));

  AMP::pout<<"Finished reseting the jacobian."<<std::endl;

  nonlinearThermalOperator->apply(rhsVec, solVec, resVec, 1.0, -1.0);

  double initialResidualNorm  = resVec->L2Norm();
  AMP::pout<<"Initial Residual Norm: "<< initialResidualNorm <<std::endl;

#ifdef USES_SILO
  siloWriter->writeFile( exeName , 0 );
#endif

  if(initialResidualNorm > 1.0e-08)
  {
    ut->failure("Nonlinear Diffusion Operator with stand alone Robin BC ");
  }
  else
  {
    ut->passes("Nonlinear Diffusion Operator with stand alone Robin BC ");
  }
  ut->passes(exeName);

  if( globalComm.getSize() == 1 ) {
#ifdef USES_SILO
    siloWriter->writeFile( exeName , 0 );
#endif
  }

}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.push_back("testNonlinearRobin-1");

    for(unsigned int i = 0; i < exeNames.size(); i++) {
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


