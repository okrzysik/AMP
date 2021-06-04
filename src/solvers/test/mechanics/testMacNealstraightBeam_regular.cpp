#include "AMP/ampmesh/libmesh/initializeLibMesh.h"
#include "AMP/ampmesh/libmesh/libmeshMesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/solvers/petsc/PetscKrylovSolver.h"
#include "AMP/solvers/petsc/PetscKrylovSolverParameters.h"
#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/ReadTestMesh.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/VectorSelector.h"
#include <memory>

#include "libmesh/mesh_communication.h"

#include <fstream>
#include <iostream>
#include <string>
#include <sys/stat.h>


static void linearElasticTest( AMP::UnitTest *ut, std::string exeName, int exampleNum )
{
    NULL_USE( exampleNum );
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName + ".txt";

    AMP::logOnlyNodeZero( log_file );

    AMP::AMP_MPI globalComm = AMP::AMP_MPI( AMP_COMM_WORLD );
    auto input_db           = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    const unsigned int mesh_dim = 3;
    libMesh::Parallel::Communicator comm( globalComm.getCommunicator() );
    auto mesh = std::make_shared<libMesh::Mesh>( comm, mesh_dim );

    std::string mesh_file = input_db->getString( "mesh_file" );
    if ( globalComm.getRank() == 0 )
        AMP::readTestMesh( mesh_file, mesh );

    libMesh::MeshCommunication().broadcast( *( mesh.get() ) );
    mesh->prepare_for_use( false );
    auto meshAdapter = std::make_shared<AMP::Mesh::libmeshMesh>( mesh, "beam" );

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
    auto bvpOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "MechanicsBVPOperator", input_db, elementPhysicsModel ) );

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
    auto dirichletVecOp = std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "Load_Boundary", input_db, dummyModel ) );
    // This has an in-place apply. So, it has an empty input variable and
    // the output variable is the same as what it is operating on.
    dirichletVecOp->setVariable( bvpOperator->getOutputVariable() );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    auto DOF_vector = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );
    auto mechSolVec =
        AMP::LinearAlgebra::createVector( DOF_vector, bvpOperator->getOutputVariable(), true );
    auto mechRhsVec = mechSolVec->cloneVector();
    auto mechResVec = mechSolVec->cloneVector();

    mechSolVec->setToScalar( 0.5 );
    mechRhsVec->setToScalar( 0.0 );
    mechResVec->setToScalar( 0.0 );

    dirichletVecOp->apply( nullVec, mechRhsVec );


    AMP::pout << "RHS Norm: " << mechRhsVec->L2Norm() << std::endl;
    AMP::pout << "Initial Solution Norm: " << mechSolVec->L2Norm() << std::endl;

    bvpOperator->residual( mechRhsVec, mechSolVec, mechResVec );

    double initResidualNorm = static_cast<double>( mechResVec->L2Norm() );
    AMP::pout << "Initial Residual Norm: " << initResidualNorm << std::endl;

    std::shared_ptr<AMP::Database> linearSolver_db = input_db->getDatabase( "LinearSolver" );

    // ---- first initialize the preconditioner
    auto pcSolver_db    = linearSolver_db->getDatabase( "Preconditioner" );
    auto pcSolverParams = std::make_shared<AMP::Solver::TrilinosMLSolverParameters>( pcSolver_db );
    pcSolverParams->d_pOperator = bvpOperator;
    auto pcSolver               = std::make_shared<AMP::Solver::TrilinosMLSolver>( pcSolverParams );

    // initialize the linear solver
    auto linearSolverParams =
        std::make_shared<AMP::Solver::PetscKrylovSolverParameters>( linearSolver_db );
    linearSolverParams->d_pOperator       = bvpOperator;
    linearSolverParams->d_comm            = globalComm;
    linearSolverParams->d_pPreconditioner = pcSolver;
    auto linearSolver = std::make_shared<AMP::Solver::PetscKrylovSolver>( linearSolverParams );

    linearSolver->setZeroInitialGuess( false );

    linearSolver->apply( mechRhsVec, mechSolVec );

    AMP::pout << "Final Solution Norm: " << mechSolVec->L2Norm() << std::endl;

    auto mechUvec = mechSolVec->select( AMP::LinearAlgebra::VS_Stride( 0, 3 ), "U" );
    auto mechVvec = mechSolVec->select( AMP::LinearAlgebra::VS_Stride( 1, 3 ), "V" );
    auto mechWvec = mechSolVec->select( AMP::LinearAlgebra::VS_Stride( 2, 3 ), "W" );

    AMP::pout << "Maximum U displacement: " << mechUvec->maxNorm() << std::endl;
    AMP::pout << "Maximum V displacement: " << mechVvec->maxNorm() << std::endl;
    AMP::pout << "Maximum W displacement: " << mechWvec->maxNorm() << std::endl;

    bvpOperator->residual( mechRhsVec, mechSolVec, mechResVec );

    double finalResidualNorm = static_cast<double>( mechResVec->L2Norm() );

    AMP::pout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    if ( finalResidualNorm > ( 1e-10 * initResidualNorm ) ) {
        ut->failure( exeName );
    } else {
        ut->passes( exeName );
    }

#ifdef USE_EXT_SILO
    auto siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    siloWriter->registerVector( mechSolVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Solution" );
    char outFileName1[256];
    sprintf( outFileName1, "undeformedBeam_%d", exampleNum );
    siloWriter->writeFile( outFileName1, 0 );
    meshAdapter->displaceMesh( mechSolVec );
    char outFileName2[256];
    sprintf( outFileName2, "deformedBeam_%d", exampleNum );
    siloWriter->writeFile( outFileName2, 0 );
#endif
}

int testMacNealstraightBeam_regular( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    auto libmeshInit =
        std::make_shared<AMP::Mesh::initializeLibMesh>( AMP::AMP_MPI( AMP_COMM_WORLD ) );

    std::vector<std::string> exeNames;

    if ( argc == 1 ) {
        exeNames.emplace_back( "testMacNealstraightBeam-regular-X-normal-mesh0" );
        exeNames.emplace_back( "testMacNealstraightBeam-regular-X-reduced-mesh0" );
        exeNames.emplace_back( "testMacNealstraightBeam-regular-Y-normal-mesh0" );
        exeNames.emplace_back( "testMacNealstraightBeam-regular-Y-reduced-mesh0" );
        exeNames.emplace_back( "testMacNealstraightBeam-regular-Z-normal-mesh0" );
        exeNames.emplace_back( "testMacNealstraightBeam-regular-Z-reduced-mesh0" );

        exeNames.emplace_back( "testMacNealstraightBeam-regular-X-normal-mesh1" );
        exeNames.emplace_back( "testMacNealstraightBeam-regular-X-reduced-mesh1" );
        exeNames.emplace_back( "testMacNealstraightBeam-regular-Y-normal-mesh1" );
        exeNames.emplace_back( "testMacNealstraightBeam-regular-Y-reduced-mesh1" );
        exeNames.emplace_back( "testMacNealstraightBeam-regular-Z-normal-mesh1" );
        exeNames.emplace_back( "testMacNealstraightBeam-regular-Z-reduced-mesh1" );

        exeNames.emplace_back( "testMacNealstraightBeam-regular-X-normal-mesh2" );
        exeNames.emplace_back( "testMacNealstraightBeam-regular-X-reduced-mesh2" );
        exeNames.emplace_back( "testMacNealstraightBeam-regular-Y-normal-mesh2" );
        exeNames.emplace_back( "testMacNealstraightBeam-regular-Y-reduced-mesh2" );
        exeNames.emplace_back( "testMacNealstraightBeam-regular-Z-normal-mesh2" );
        exeNames.emplace_back( "testMacNealstraightBeam-regular-Z-reduced-mesh2" );
    } else {
        for ( int i = 1; i < argc; i += 3 ) {
            char inpName[100];
            sprintf( inpName,
                     "testMacNealstraightBeam-regular-%s-%s-mesh%d",
                     argv[i],
                     argv[i + 1],
                     atoi( argv[i + 2] ) );
            exeNames.emplace_back( inpName );
        } // end for i
    }

    for ( size_t i = 0; i < exeNames.size(); i++ ) {
        try {
            linearElasticTest( &ut, exeNames[i], i );
            AMP::pout << exeNames[i] << " had " << ut.NumFailGlobal() << " failures." << std::endl;
        } catch ( std::exception &err ) {
            AMP::pout << "ERROR: " << err.what() << std::endl;
        } catch ( ... ) {
            AMP::pout << "ERROR: "
                      << "An unknown exception was thrown." << std::endl;
        }
    } // end for i

    ut.report();
    int num_failed = ut.NumFailGlobal();

    libmeshInit.reset();

    AMP::AMPManager::shutdown();
    return num_failed;
}
