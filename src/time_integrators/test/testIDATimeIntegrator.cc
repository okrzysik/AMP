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
#include "utils/Writer.h"
#include "vectors/VectorBuilder.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"

#include "materials/Material.h"

#include "vectors/Variable.h"
#include "vectors/Vector.h"

#include "operators/diffusion/DiffusionLinearElement.h"
#include "operators/diffusion/DiffusionTransportModel.h"
#include "operators/VolumeIntegralOperator.h"
#include "operators/MassLinearElement.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/MassLinearFEOperator.h"
#include "operators/boundary/DirichletMatrixCorrection.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/LinearBVPOperator.h"

#include "solvers/TrilinosMLSolver.h"

#include "time_integrators/sundials/IDATimeIntegrator.h"
#include "time_integrators/sundials/IDATimeOperator.h"


void IDATimeIntegratorTest(AMP::UnitTest *ut )
{
    boost::weak_ptr<AMP::Mesh::Mesh>  mesh_ref;
    std::string input_file = "input_testIDA-LinearBVPOperator-1";
    std::string log_file = "output_testIDA-LinearBVPOperator-1";
    AMP::PIO::logOnlyNodeZero(log_file);

    {
        // Read the input file
        boost::shared_ptr<AMP::InputDatabase>  input_db ( new AMP::InputDatabase ( "input_db" ) );
        AMP::InputManager::getManager()->parseInputFile ( input_file , input_db );

        // Get the Mesh database and create the mesh parameters
        boost::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
        boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
        params->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));

        // Create the meshes from the input database
        AMP::Mesh::Mesh::shared_ptr mesh = AMP::Mesh::Mesh::buildMesh(params);
        mesh_ref = boost::weak_ptr<AMP::Mesh::Mesh>( mesh );

        // Create a DOF manager for a nodal vector 
        int DOFsPerNode = 1;
        //int DOFsPerElement = 8;
        int nodalGhostWidth = 1;
        //int gaussPointGhostWidth = 1;
        bool split = true;
        AMP::Discretization::DOFManager::shared_ptr nodalDofMap      = AMP::Discretization::simpleDOFManager::create(mesh, AMP::Mesh::Vertex, nodalGhostWidth,      DOFsPerNode,    split);
        //AMP::Discretization::DOFManager::shared_ptr gaussPointDofMap = AMP::Discretization::simpleDOFManager::create(mesh, AMP::Mesh::Volume, gaussPointGhostWidth, DOFsPerElement, split);
        
        // create a linear BVP operator
        boost::shared_ptr<AMP::Operator::LinearBVPOperator> IDARhsOperator;
        boost::shared_ptr<AMP::Operator::LinearBVPOperator> linearPCOperator;
        boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel;
        IDARhsOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(mesh,"LinearOperator",input_db,elementModel) );
     
        // create a mass linear BVP operator
        boost::shared_ptr<AMP::Operator::ElementPhysicsModel> massElementModel;
        boost::shared_ptr<AMP::Operator::LinearBVPOperator> massOperator;
        massOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(mesh,"MassLinearOperator",input_db,massElementModel) );
      
        // create vectors for initial conditions (IC) and time derivative at IC
        AMP::LinearAlgebra::Variable::shared_ptr outputVar = IDARhsOperator->getOutputVariable();
        AMP::LinearAlgebra::Vector::shared_ptr initialCondition      = AMP::LinearAlgebra::createVector( nodalDofMap, outputVar );
        AMP::LinearAlgebra::Vector::shared_ptr initialConditionPrime = AMP::LinearAlgebra::createVector( nodalDofMap, outputVar );
        initialCondition->zero();
        initialConditionPrime->zero();

        // get the ida database
        AMP_INSIST(input_db->keyExists("IDATimeIntegrator"), "Key ''IDATimeIntegrator'' is missing!");
        boost::shared_ptr<AMP::Database> ida_db = input_db->getDatabase("IDATimeIntegrator");
        boost::shared_ptr<AMP::Database> pcSolver_db = ida_db->getDatabase("Preconditioner");
        boost::shared_ptr<AMP::Solver::SolverStrategyParameters> pcSolverParams(new AMP::Solver::SolverStrategyParameters(pcSolver_db));
        boost::shared_ptr<AMP::Solver::TrilinosMLSolver> pcSolver(new AMP::Solver::TrilinosMLSolver(pcSolverParams));

        // create the IDA time integrator
        boost::shared_ptr<AMP::TimeIntegrator::IDATimeIntegratorParameters> time_Params( 
            new AMP::TimeIntegrator::IDATimeIntegratorParameters(ida_db));
        time_Params->d_pMassOperator = massOperator;
        time_Params->d_operator = IDARhsOperator;
        time_Params->d_pPreconditioner = pcSolver;
        time_Params->d_ic_vector = initialCondition;    
        time_Params->d_ic_vector_prime = initialConditionPrime;
        time_Params->d_object_name = "IDATimeIntegratorParameters";
        
        boost::shared_ptr<AMP::TimeIntegrator::IDATimeIntegrator> pIDATimeIntegrator(new AMP::TimeIntegrator::IDATimeIntegrator(time_Params));
        if ( pIDATimeIntegrator.get() == NULL ) {
            ut->failure("Testing IDATimeIntegrator's constructor");
        } else {
            ut->passes("Tested IDATimeIntegrator's constructor");
        }

    }

    // Check that everything was destroyed
    if ( mesh_ref.expired() )
        ut->passes("IDATimeIntegrator delete");
    else
        ut->failure("IDATimeIntegrator delete");

}



// Main
int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;
    
    IDATimeIntegratorTest(&ut);        

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;    
    
}   








