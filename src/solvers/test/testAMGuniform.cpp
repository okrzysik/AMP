
#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/libmesh/initializeLibMesh.h"
#include "AMP/ampmesh/libmesh/libmeshMesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/ReadTestMesh.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include "libmesh/mesh_communication.h"
#include <string>

#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"

void myTest( AMP::UnitTest *ut )
{
    std::string exeName( "testAMGuniform" );
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;
    AMP::logAllNodes( log_file );

    AMP::AMP_MPI globalComm = AMP::AMP_MPI( AMP_COMM_WORLD );


    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    const unsigned int mesh_dim = 3;
    libMesh::Parallel::Communicator comm( globalComm.getCommunicator() );
    auto mesh             = std::make_shared<libMesh::Mesh>( comm, mesh_dim );
    std::string mesh_file = input_db->getString( "mesh_file" );
    if ( globalComm.getRank() == 0 ) {
        AMP::readTestMesh( mesh_file, mesh );
    } // end if root processor
    libMesh::MeshCommunication().broadcast( *( mesh.get() ) );
    // mesh->prepare_for_use(false);
    mesh->prepare_for_use( true );
    auto meshAdapter = std::make_shared<AMP::Mesh::libmeshMesh>( mesh, "uniform" );

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
    auto bvpOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "LinearBVPOperator", input_db, elementPhysicsModel ) );

    auto mat     = bvpOperator->getMatrix();
    size_t matSz = mat->numGlobalRows();
    for ( size_t i = 0; i < matSz; ++i ) {
        std::vector<size_t> cols;
        std::vector<double> vals;
        mat->getRowByGlobalID( i, cols, vals );
        for ( size_t j = 0; j < cols.size(); ++j ) {
            std::cout << "A[" << i << "][" << ( cols[j] ) << "] = " << std::setprecision( 15 )
                      << ( vals[j] ) << std::endl;
        } // end j
        std::cout << std::endl;
    } // end i

    int DOFsPerNode = 1;
    // int DOFsPerElement = 8;
    int nodalGhostWidth = 1;
    bool split          = true;
    auto nodalDofMap    = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    auto solVec = AMP::LinearAlgebra::createVector( nodalDofMap, bvpOperator->getOutputVariable() );
    auto rhsVec = solVec->cloneVector();

    solVec->setRandomValues();
    bvpOperator->apply( solVec, rhsVec );
    solVec->zero();

    auto mlSolver_db    = input_db->getDatabase( "LinearSolver" );
    auto mlSolverParams = std::make_shared<AMP::Solver::SolverStrategyParameters>( mlSolver_db );
    mlSolverParams->d_pOperator = bvpOperator;
    auto mlSolver               = std::make_shared<AMP::Solver::TrilinosMLSolver>( mlSolverParams );

    mlSolver->setZeroInitialGuess( true );

    mlSolver->solve( rhsVec, solVec );

    ut->passes( exeName );
}


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    auto libmeshInit =
        std::make_shared<AMP::Mesh::initializeLibMesh>( AMP::AMP_MPI( AMP_COMM_WORLD ) );

    myTest( &ut );

    ut.report();
    int num_failed = ut.NumFailGlobal();

    libmeshInit.reset();

    AMP::AMPManager::shutdown();
    return num_failed;
}
