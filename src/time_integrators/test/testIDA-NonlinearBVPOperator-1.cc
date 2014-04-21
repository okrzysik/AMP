#include <string>
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/AMP_MPI.h"
#include "materials/Material.h"
#include "boost/shared_ptr.hpp"
#include "utils/InputDatabase.h"
#include "utils/Utilities.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/Database.h"
#include "operators/NeutronicsRhs.h"
#include "vectors/Variable.h"

#include "utils/Writer.h"
#include "ampmesh/Mesh.h"

#include "vectors/Vector.h"
#include "operators/libmesh/VolumeIntegralOperator.h"
#include "operators/NeutronicsRhs.h"
#include "operators/NonlinearBVPOperator.h"
/* libMesh files */
#include "libmesh_config.h"
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "equation_systems.h"
#include "fe.h"
#include "quadrature_gauss.h"
#include "dof_map.h"
#include "sparse_matrix.h"
#include "petsc_matrix.h"
#include "petsc_vector.h"
#include "dense_matrix.h"
#include "linear_implicit_system.h"
#include "elem.h"

#include "vectors/VectorBuilder.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"

#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/libmesh/MassLinearFEOperator.h"

#include "operators/LinearBVPOperator.h"

#include "solvers/trilinos/TrilinosMLSolver.h"
#include "time_integrators/sundials/IDATimeIntegrator.h"
#include "time_integrators/sundials/IDATimeOperator.h"

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
    std::string input_file = "input_testIDA-NonlinearBVPOperator-1";
    std::string log_file = "output_testIDA-NonlinearBVPOperator-1";
    
    AMP::PIO::logOnlyNodeZero(log_file);
    
    boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
    AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
    input_db->printClassData(AMP::plog);
    
    AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
    boost::shared_ptr<AMP::Database>  mesh_db = input_db->getDatabase("Mesh");
    boost::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(new AMP::Mesh::MeshParameters(mesh_db));
    mgrParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
    boost::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh(mgrParams);
    
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

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a nonlinear BVP operator for nonlinear BVP operator
    AMP_INSIST( input_db->keyExists("NonlinearOperator"), "key missing!" );

    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel;
    boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearOperator = boost::dynamic_pointer_cast<
      AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
											"NonlinearOperator",
											input_db,
											elementModel));

    AMP::LinearAlgebra::Variable::shared_ptr outputVar = nonlinearOperator->getOutputVariable();
    // ---------------------------------------------------------------------------------------
    // create a linear BVP operator
    boost::shared_ptr<AMP::Operator::LinearBVPOperator> linearOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
                                                                                               AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
																	      "LinearOperator",
																	      input_db,
																	      elementModel));

    // ---------------------------------------------------------------------------------------
    // create a mass linear BVP operator
    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> massElementModel;
    boost::shared_ptr<AMP::Operator::LinearBVPOperator> massOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
                                                                                               AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
																	      "MassLinearOperator",
																	      input_db,
																	      massElementModel));

    // ---------------------------------------------------------------------------------------
    //  create neutronics source
    AMP_INSIST(input_db->keyExists("NeutronicsOperator"), "Key ''NeutronicsOperator'' is missing!");
    boost::shared_ptr<AMP::Database>  neutronicsOp_db = input_db->getDatabase("NeutronicsOperator");
    boost::shared_ptr<AMP::Operator::NeutronicsRhsParameters> neutronicsParams(new AMP::Operator::NeutronicsRhsParameters( neutronicsOp_db ));
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
    
    AMP::LinearAlgebra::Vector::shared_ptr  initialCondition = AMP::LinearAlgebra::createVector( nodalDofMap, outputVar );
    AMP::LinearAlgebra::Vector::shared_ptr  initialConditionPrime = AMP::LinearAlgebra::createVector( nodalDofMap, outputVar );
    AMP::LinearAlgebra::Vector::shared_ptr   f = AMP::LinearAlgebra::createVector( nodalDofMap, outputVar );
    
    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // set initial conditions, initialize created vectors
    
    AMP::Mesh::MeshIterator node = meshAdapter->getIterator( AMP::Mesh::Vertex, 0 );
    AMP::Mesh::MeshIterator end_node = node.end();
    
    int counter=0;     
    for( ; node != end_node ; ++node)
    {
        counter++;
        
        std::vector<size_t> bndGlobalIds;
        nodalDofMap->getDOFs( node->globalID() , bndGlobalIds );
        
        std::vector<double> pt = node->coord();
        double px = pt[0];
        double py = pt[1];
        double pz = pt[2];
        
        double val = __INIT_FN__(px, py, pz, 0);

        for(unsigned int i = 0; i < bndGlobalIds.size(); i++)
        {
          initialCondition->setValueByGlobalID(bndGlobalIds[i], val);
          // ** please do not set the time derivative to be non-zero!!
          // ** as this causes trouble with the boundary - BP, 07/16/2010
          initialConditionPrime->setValueByGlobalID(bndGlobalIds[i], 0.0);

        }//end for i
    }//end for node
    initialCondition->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
    initialConditionPrime->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
    
    std::cout << "With Counter "<< counter << " Max initial temp "<< initialCondition->max()<< " Min initial temp "<<initialCondition->min() << std::endl;

    // create a copy of the rhs which can be modified at each time step (maybe)
    f->copyVector(powerInWattsVec);
    // modify the rhs to take into account boundary conditions
    nonlinearOperator->modifyRHSvector(f);
    nonlinearOperator->modifyInitialSolutionVector(initialCondition);
    
    // ---------------------------------------------------------------------------------------
    // create a linear time operator
    boost::shared_ptr<AMP::InputDatabase> timeOperator_db(new AMP::InputDatabase("TimeOperatorDatabase"));
    timeOperator_db->putDouble("CurrentDt", 0.01);
    timeOperator_db->putString("name", "TimeOperator");
    timeOperator_db->putBool("bLinearMassOperator", true);
    timeOperator_db->putBool("bLinearRhsOperator", false);
    timeOperator_db->putDouble("ScalingFactor", 1.0/0.01);
    timeOperator_db->putDouble("CurrentTime", .0);
    
    boost::shared_ptr<AMP::TimeIntegrator::TimeOperatorParameters> timeOperatorParameters(new AMP::TimeIntegrator::TimeOperatorParameters(timeOperator_db));
    timeOperatorParameters->d_pRhsOperator = linearOperator;
    timeOperatorParameters->d_pMassOperator = massOperator;
    //timeOperatorParameters->d_pMassOperator = massLinearOperator;
    boost::shared_ptr<AMP::Operator::Operator> linearTimeOperator( new AMP::TimeIntegrator::LinearTimeOperator(timeOperatorParameters));

    AMP::LinearAlgebra::Vector::shared_ptr  residualVec = AMP::LinearAlgebra::createVector( nodalDofMap, outputVar );
    
    linearOperator->apply(nullVec, initialCondition, residualVec, 1., 0 );
    std::cout << "Residual Norm of linearTimeOp apply : "<< residualVec->L2Norm()  << std::endl;
    
    massOperator->apply(nullVec, initialCondition, residualVec, 1., 0 );
    std::cout << "Residual Norm of linearTimeOp apply : "<< residualVec->L2Norm()  << std::endl;
    
    linearTimeOperator->apply(nullVec, initialCondition, residualVec, 1., 0 );
    std::cout << "Residual Norm of linearTimeOp apply : "<< residualVec->L2Norm()  << std::endl;

    // ---------------------------------------------------------------------------------------
    // create a preconditioner
    
    // get the ida database
    AMP_INSIST(input_db->keyExists("IDATimeIntegrator"), "Key ''IDATimeIntegrator'' is missing!");
    boost::shared_ptr<AMP::Database> ida_db = input_db->getDatabase("IDATimeIntegrator");
    boost::shared_ptr<AMP::Database> pcSolver_db = ida_db->getDatabase("Preconditioner");
    boost::shared_ptr<AMP::Solver::SolverStrategyParameters> pcSolverParams(new AMP::Solver::SolverStrategyParameters(pcSolver_db));
    pcSolverParams->d_pOperator = linearTimeOperator;
    
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
    
#ifdef USE_EXT_SILO
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter("Silo");
    siloWriter->registerMesh( meshAdapter );

    siloWriter->registerVector( initialCondition,                 meshAdapter, AMP::Mesh::Vertex, "InitialSolution" );

    siloWriter->writeFile( input_file , 0 );
#endif

    // ---------------------------------------------------------------------------------------
    // create the IDA time integrator
    boost::shared_ptr<AMP::TimeIntegrator::IDATimeIntegratorParameters> time_Params( new AMP::TimeIntegrator::IDATimeIntegratorParameters(ida_db));
    
    if( (time_Params.get()) == NULL ) {
        ut->failure("Testing IDATimeIntegratorParameters' Constructor");
    } else {
        ut->passes("Testing IDATimeIntegratorParameters' Constructor");
    }
    
    time_Params->d_pMassOperator = massOperator;
    //time_Params->d_pMassOperator = massLinearOperator;
    time_Params->d_operator = nonlinearOperator;
    time_Params->d_pPreconditioner = pcSolver;
    
    time_Params->d_ic_vector = initialCondition;    
    time_Params->d_ic_vector_prime = initialConditionPrime;
    
    time_Params->d_pSourceTerm = f;
    time_Params->d_object_name = "IDATimeIntegratorParameters";
    
    std::cout << "Before IDATimeIntegrator" << std::endl;    
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
    double min=0;
    int j=1;
    while(pIDATimeIntegrator->getCurrentTime() < pIDATimeIntegrator->getFinalTime())
    {
        retval = pIDATimeIntegrator->advanceSolution(pIDATimeIntegrator->getCurrentDt(), 0);
        //pIDATimeIntegrator->updateSolution();
        current_time = pIDATimeIntegrator->getCurrentTime();
        
        std::cout << j++ << "-th timestep" << std::endl;
        if(retval == 0) {
            ut->passes("Testing IDATimeIntegrator's advanceSolution. PASS!!");
        } else {
            ut->failure("Tested IDATimeIntegrator's advanceSolution. FAIL!!");
        }
        
        max = pIDATimeIntegrator->getCurrentSolution()->max();
        min = pIDATimeIntegrator->getCurrentSolution()->min();
        
        std::cout << "current_time = " << current_time << std::endl;
        std::cout << "max val of the current solution = " << max << std::endl;
        std::cout << "min val of the current solution = " << min << std::endl;
    }
    
    
#ifdef USE_EXT_SILO

    AMP::LinearAlgebra::Vector::shared_ptr pSolution=pIDATimeIntegrator->getCurrentSolution();
    siloWriter->registerVector( pSolution,                 meshAdapter, AMP::Mesh::Vertex, "Solution" );

    siloWriter->writeFile( input_file , 1 );
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









