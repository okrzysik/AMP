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
#include "operators/boundary/DirichletMatrixCorrection.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/LinearBVPOperator.h"

#include "solvers/TrilinosMLSolver.h"

#include "time_integrators/IDATimeIntegrator.h"
#include "time_integrators/IDATimeOperator.h"

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

#define __PI__ 3.14159265

//#define __INIT_FN__(x, y, z, t) ( exp(- __PI__ * __PI__ * t) * sin(__PI__ * x) * sin(__PI__ * y) * sin(__PI__ * z) )
//#define __INIT_FN__(x,y,z,t) ( exp(-0.015 *  __PI__ * __PI__ * t) * cos(0.1 * __PI__ * x) * cos(0.1 * __PI__ * y) * cos(0.05 * __PI__ * z) )
#define __INIT_FN__(x,y,z,t) ( 750.0+ 10000.0*(0.5+ x) * (0.5 -x) *(0.5+ y) * (0.5 -y) *(0.5+ z) * (0.5 -z) )

void IDATimeIntegratorTest(AMP::UnitTest *ut )
{
    std::string input_file = "input_testIDA-LinearBVPOperator-1";
    std::string log_file = "output_testIDA-LinearBVPOperator-1";
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
    AMP::Mesh::Mesh::shared_ptr meshAdapter = manager->Subset( "ida" );

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
//--------------------------------------------------
  
    // create a linear BVP operator
    boost::shared_ptr<AMP::Operator::LinearBVPOperator> IDARhsOperator;
    boost::shared_ptr<AMP::Operator::LinearBVPOperator> linearPCOperator;
    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel;
    IDARhsOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
																"LinearOperator",
																input_db,
																elementModel));
 
    // create a mass linear BVP operator
    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> massElementModel;
    boost::shared_ptr<AMP::Operator::LinearBVPOperator> massOperator;
    massOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
															      "MassLinearOperator",
															      input_db,
															      massElementModel));

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
  
    // ---------------------------------------------------------------------------------------
    // create vectors for initial conditions (IC) and time derivative at IC
    AMP::LinearAlgebra::Variable::shared_ptr inputVar = IDARhsOperator->getInputVariable();
    AMP::LinearAlgebra::Variable::shared_ptr outputVar = IDARhsOperator->getOutputVariable();
  
    AMP::LinearAlgebra::Vector::shared_ptr initialCondition      = AMP::LinearAlgebra::createVector( nodalDofMap, outputVar );
    AMP::LinearAlgebra::Vector::shared_ptr initialConditionPrime = AMP::LinearAlgebra::createVector( nodalDofMap, outputVar );
    AMP::LinearAlgebra::Vector::shared_ptr   f                   = AMP::LinearAlgebra::createVector( nodalDofMap, outputVar );

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // set initial conditions, initialize created vectors
  int zeroGhostWidth = 1;
  AMP::Mesh::MeshIterator  node = meshAdapter->getSurfaceIterator(AMP::Mesh::Vertex, zeroGhostWidth);
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
      cout << "val = " << val << endl;
      
      cout << "counter = " << counter << "gid.size() = " << gid.size() << endl;
      for(unsigned int i = 0; i < gid.size(); i++)
    {
      initialCondition->setValueByGlobalID(gid[i], val);
    }//end for i
    }//end for node
  initialConditionPrime->zero();
  
  // create a copy of the rhs which can be modified at each time step (maybe)
  f->copyVector(powerInWattsVec);
  // modify the rhs to take into account boundary conditions
  IDARhsOperator->modifyRHSvector(f);

  // ---------------------------------------------------------------------------------------
  // create a preconditioner

  // get the ida database
  AMP_INSIST(input_db->keyExists("IDATimeIntegrator"), "Key ''IDATimeIntegrator'' is missing!");
  boost::shared_ptr<AMP::Database> ida_db = input_db->getDatabase("IDATimeIntegrator");
  boost::shared_ptr<AMP::Database> pcSolver_db = ida_db->getDatabase("Preconditioner");
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> pcSolverParams(new AMP::Solver::SolverStrategyParameters(pcSolver_db));
  
  if(pcSolverParams.get() == NULL) {
    ut->failure("Testing SolverStrategyParameters's constructor: FAIL");
  } else {
    ut->passes("Testing SolverStrategyParameters's constructor: PASS");
  }
  
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> pcSolver(new AMP::Solver::TrilinosMLSolver(pcSolverParams));
  
  if(pcSolver.get() == NULL) {
    ut->failure("Testing TrilinosMLSolver's constructor: FAIL");
  } else {
    ut->passes("Testing TrilinosMLSolver's constructor: PASS");
  }

  // ---------------------------------------------------------------------------------------
  // create the IDA time integrator
  boost::shared_ptr<AMP::TimeIntegrator::IDATimeIntegratorParameters> time_Params( new AMP::TimeIntegrator::IDATimeIntegratorParameters(ida_db));
  
  if( (time_Params.get()) == NULL ) {
    ut->failure("Testing IDATimeIntegratorParameters' Constructor");
  } else {
    ut->passes("Testing IDATimeIntegratorParameters' Constructor");
  }
    
  time_Params->d_pMassOperator = massOperator;
  time_Params->d_operator = IDARhsOperator;
  time_Params->d_pPreconditioner = pcSolver;
  
  time_Params->d_ic_vector = initialCondition;    
  time_Params->d_ic_vector_prime = initialConditionPrime;
  
  time_Params->d_pSourceTerm = f;
  time_Params->d_object_name = "IDATimeIntegratorParameters";
    
  cout << "Before IDATimeIntegrator" << endl;    
#ifdef USE_SUNDIALS
  boost::shared_ptr<AMP::TimeIntegrator::IDATimeIntegrator> pIDATimeIntegrator(new AMP::TimeIntegrator::IDATimeIntegrator(time_Params));
  
  if(pIDATimeIntegrator.get() == NULL) {
    ut->failure("Testing IDATimeIntegrator's constructor");
  } else {
    ut->passes("Tested IDATimeIntegrator's constructor");
  }
  // ---------------------------------------------------------------------------------------
  // step in time
  int retval=0;
  double current_time=0;
  double max=0;
  //double abs_error=0.0;
  double min=0;
  //double rel_error=0.0;
  //double exact_sol=0.0;

  int j=1;
  while(pIDATimeIntegrator->getCurrentTime() < pIDATimeIntegrator->getFinalTime())
    {
      retval = pIDATimeIntegrator->advanceSolution(pIDATimeIntegrator->getCurrentDt(), 0);
      //pIDATimeIntegrator->updateSolution();
      current_time = pIDATimeIntegrator->getCurrentTime();
      
      cout << j++ << "-th timestep" << endl;
      if(retval == 0) {
    ut->passes("Testing IDATimeIntegrator's advanceSolution. PASS!!");
      } else {
    ut->failure("Tested IDATimeIntegrator's advanceSolution. FAIL!!");
      }
      
      max = pIDATimeIntegrator->getCurrentSolution()->max();
      min = pIDATimeIntegrator->getCurrentSolution()->min();
      
      //      exact_sol = exp(-0.015 *__PI__ * __PI__ * current_time);
      //exact_sol = exp(
      //      abs_error = exact_sol-max;
      //      rel_error = abs_error/exact_sol;
      
      cout << "current_time = " << current_time << endl;
      cout << "max val of the current solution = " << max << endl;
      cout << "min val of the current solution = " << min << endl;
      //      cout << "exact solution = " << exact_sol << endl;
      //      cout << "absolute error = " << abs_error << endl;
      //      cout << "relative error = " << rel_error << endl;
    }
  
#else
    ut->passes("IDA will not fail a test if there is no IDA.");
#endif  
  
    if (ut->NumFailLocal() == 0)
    {
      ut->passes("testIDATimeIntegrator successful");
    }
}


//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;
    
    try {        
        IDATimeIntegratorTest(&ut);        
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

//---------------------------------------------------------------------------//
//                        end of SundialsVectorTest.cc
//---------------------------------------------------------------------------//








