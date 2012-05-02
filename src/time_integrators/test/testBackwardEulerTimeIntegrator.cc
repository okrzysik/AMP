#include <string>
#include "boost/shared_ptr.hpp"

#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/InputDatabase.h"
#include "utils/Utilities.h"
#include "utils/PIO.h"
#include "utils/Database.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/SiloIO.h"
#include "vectors/VectorBuilder.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"

#include "materials/Material.h"

#include "vectors/Variable.h"
#include "vectors/Vector.h"

#include "operators/NeutronicsRhs.h"
#include "operators/diffusion/DiffusionLinearElement.h"
#include "operators/diffusion/DiffusionTransportModel.h"
#include "operators/VolumeIntegralOperator.h"
#include "operators/NeutronicsRhs.h"
#include "operators/MassLinearElement.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/MassLinearFEOperator.h"
#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"

#include "solvers/TrilinosMLSolver.h"

#include "time_integrators/ImplicitTimeIntegratorParameters.h"
#include "time_integrators/BackwardEulerTimeIntegrator.h"
#include "time_integrators/BackwardEulerTimeOperator.h"

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

#define __PI__ 3.14159265
#define __INIT_FN__(x,y,z,t) ( 750.0+ 10000.0*(0.5+ x) * (0.5 -x) *(0.5+ y) * (0.5 -y) *(0.5+ z) * (0.5 -z) )

void BackwardEulerTimeIntegrator(AMP::UnitTest *ut )
{
    std::string input_file = "input_testBackwardEulerTimeIntegrator";
    std::string log_file = "output_testBackwardEulerTimeIntegrator";
    AMP::PIO::logOnlyNodeZero(log_file);

    // Read the input file
    boost::shared_ptr<AMP::InputDatabase>  input_db ( new AMP::InputDatabase ( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile ( input_file , input_db );

    // Get the Mesh database and create the mesh parameters
    boost::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
    boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
    params->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));

    // Create the meshes from the input database
    AMP::Mesh::Mesh::shared_ptr manager = AMP::Mesh::Mesh::buildMesh(params);
    AMP::Mesh::Mesh::shared_ptr meshAdapter = manager->Subset( "Cube" );

//--------------------------------------------------
// Create a DOF manager for a nodal vector 
//--------------------------------------------------
    int DOFsPerNode = 1;
    int DOFsPerElement = 8;
    int nodalGhostWidth = 1;
    int gaussPointGhostWidth = 1;
    bool split = true;
    AMP::Discretization::DOFManager::shared_ptr nodalDofMap      = AMP::Discretization::simpleDOFManager::create(meshAdapter, AMP::Mesh::Vertex, nodalGhostWidth,      DOFsPerNode,    split);
    AMP::Discretization::DOFManager::shared_ptr gaussPointDofMap = AMP::Discretization::simpleDOFManager::create(meshAdapter, AMP::Mesh::Volume, gaussPointGhostWidth, DOFsPerElement, split);

    // create a linear BVP operator
    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel;
    boost::shared_ptr<AMP::Operator::LinearBVPOperator> diffusionOperator ;
    diffusionOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter, "LinearOperator", input_db, elementModel));
 
    // create a mass linear BVP operator
    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> massElementModel;
    boost::shared_ptr<AMP::Operator::LinearBVPOperator> massOperator ;
    massOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter, "MassLinearOperator", input_db, massElementModel));

    //  create neutronics source
    AMP_INSIST(input_db->keyExists("NeutronicsOperator"), "Key ''NeutronicsOperator'' is missing!");
    boost::shared_ptr<AMP::Database>  neutronicsOp_db = input_db->getDatabase("NeutronicsOperator");
    boost::shared_ptr<AMP::Operator::NeutronicsRhsParameters> neutronicsParams(new AMP::Operator::NeutronicsRhsParameters( neutronicsOp_db ));
    neutronicsParams->d_Mesh = meshAdapter;
    boost::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOperator(new AMP::Operator::NeutronicsRhs( neutronicsParams ));
  
    AMP::LinearAlgebra::Variable::shared_ptr SpecificPowerVar = neutronicsOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr   SpecificPowerVec = AMP::LinearAlgebra::createVector( gaussPointDofMap, SpecificPowerVar );

    // create the following shared pointers for ease of use
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    neutronicsOperator->apply(nullVec, nullVec, SpecificPowerVec, 1., 0.);
  
    //  Integrate Nuclear Rhs over Density * Volume //
    AMP_INSIST( input_db->keyExists("VolumeIntegralOperator"), "key missing!" );
  
    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> sourceTransportModel;
    boost::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator = boost::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
																							      "VolumeIntegralOperator",
																							      input_db,
																							      sourceTransportModel));
  
    // Create the power (heat source) vector.
    AMP::LinearAlgebra::Variable::shared_ptr powerInWattsVar = sourceOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr   powerInWattsVec = AMP::LinearAlgebra::createVector( nodalDofMap, powerInWattsVar );
    powerInWattsVec->zero();
  
    // convert the vector of specific power to power for a given basis.
    sourceOperator->apply(nullVec, SpecificPowerVec, powerInWattsVec, 1., 0.);
  
    AMP::LinearAlgebra::Variable::shared_ptr outputVar = diffusionOperator->getOutputVariable();
  
    AMP::LinearAlgebra::Vector::shared_ptr initialCondition      = AMP::LinearAlgebra::createVector( nodalDofMap, outputVar );
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec                = AMP::LinearAlgebra::createVector( nodalDofMap, outputVar );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // set initial conditions, initialize created vectors
    int zeroGhostWidth = 0;
    AMP::Mesh::MeshIterator  node = meshAdapter->getIterator(AMP::Mesh::Vertex, zeroGhostWidth);
    AMP::Mesh::MeshIterator  end_node = node.end();

    int counter=0;     
    for( ; node != end_node ; ++node)
    {
      counter+=1;

      std::vector<size_t> gid;
      nodalDofMap->getDOFs ( node->globalID() , gid);

      double px = ( node->coord() )[0];
      double py = ( node->coord() )[1];
      double pz = ( node->coord() )[2];

      double val = __INIT_FN__(px, py, pz, 0);
      for(unsigned int i = 0; i < gid.size(); i++)
      {
        initialCondition->setValueByGlobalID(gid[i], val);
      }//end for i
    }//end for node

    diffusionOperator->modifyRHSvector( rhsVec );

    boost::shared_ptr<AMP::Database> pcSolver_db = input_db->getDatabase("Solver");
    boost::shared_ptr<AMP::Solver::SolverStrategyParameters> pcSolverParams(new AMP::Solver::SolverStrategyParameters(pcSolver_db));
    boost::shared_ptr<AMP::Solver::TrilinosMLSolver> pcSolver(new AMP::Solver::TrilinosMLSolver(pcSolverParams));

    boost::shared_ptr<AMP::Database> timeIntegrator_db = input_db->getDatabase("BDFTimeIntegrator");
    boost::shared_ptr<AMP::TimeIntegrator::ImplicitTimeIntegratorParameters> time_Params( new AMP::TimeIntegrator::ImplicitTimeIntegratorParameters(timeIntegrator_db));
    time_Params->d_pMassOperator = massOperator;
    time_Params->d_operator = diffusionOperator;
    time_Params->d_solver   = pcSolver;

    time_Params->d_ic_vector = initialCondition;    

    time_Params->d_pSourceTerm = rhsVec ;
    time_Params->d_object_name = "ImplicitTimeIntegratorParameters";

    boost::shared_ptr<AMP::TimeIntegrator::BackwardEulerTimeIntegrator> BDFTimeIntegrator(new AMP::TimeIntegrator::BackwardEulerTimeIntegrator(time_Params));

    if(BDFTimeIntegrator.get() == NULL) {
      ut->failure("Testing BDFTimeIntegrator's constructor");
    } else {
      ut->passes("Tested BDFTimeIntegrator's constructor");
    }

    double current_time=0, max, min;
    int j=0;
    while(BDFTimeIntegrator->getCurrentTime() < BDFTimeIntegrator->getFinalTime())
    {
      BDFTimeIntegrator->advanceSolution(BDFTimeIntegrator->getCurrentDt(), 0);
      current_time = BDFTimeIntegrator->getCurrentTime();

      cout << j++ << "-th timestep" << endl;

      max = BDFTimeIntegrator->getCurrentSolution()->max();
      min = BDFTimeIntegrator->getCurrentSolution()->min();

      cout << "current_time = " << current_time << endl;
      cout << "max val of the current solution = " << max << endl;
      cout << "min val of the current solution = " << min << endl;
    }

    if (ut->NumFailLocal() == 0)
    {
      ut->passes("test Backward Euler Time Intgrator successful");
    }
}


//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  try {        
    BackwardEulerTimeIntegrator(&ut);        
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
