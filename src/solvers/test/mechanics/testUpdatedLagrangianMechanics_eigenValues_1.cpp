#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/libmesh/libmeshMesh.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/operators/BVPOperatorParameters.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletMatrixCorrection.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/mechanics/MechanicsLinearFEOperator.h"
#include "AMP/operators/mechanics/MechanicsNonlinearFEOperator.h"
#include "AMP/operators/mechanics/ThermalStrainMaterialModel.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/ReadTestMesh.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/VectorBuilder.h"
#include "libmesh/mesh_communication.h"

#include <iostream>
#include <string>

static void myTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "log_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );

    // Read the input file
    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    std::string mesh_file       = input_db->getString( "mesh_file" );
    const unsigned int mesh_dim = 3;
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    libMesh::Parallel::Communicator comm( globalComm.getCommunicator() );
    auto mesh = std::make_shared<libMesh::Mesh>( comm, mesh_dim );
    AMP::readTestMesh( mesh_file, mesh );
    libMesh::MeshCommunication().broadcast( *( mesh.get() ) );
    mesh->prepare_for_use( false );

    auto meshAdapter = std::make_shared<AMP::Mesh::libmeshMesh>( mesh, "TestMesh" );

    AMP_INSIST( input_db->keyExists( "NumberOfLoadingSteps" ),
                "Key ''NumberOfLoadingSteps'' is missing!" );

    AMP_INSIST( input_db->keyExists( "OutputFileName" ), "Key ''OutputFileName'' is missing!" );
    std::string outFileName = input_db->getString( "OutputFileName" );

    FILE *fp;
    fp = fopen( outFileName.c_str(), "w" );
    fprintf( fp, "clc; \n clear; \n A = zeros(24, 24); \n \n" );

    // Create a nonlinear BVP operator for mechanics
    AMP_INSIST( input_db->keyExists( "NonlinearMechanicsOperator" ), "key missing!" );
    auto nonlinearMechanicsBVPoperator =
        std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "NonlinearMechanicsOperator", input_db ) );
    auto nonlinearMechanicsVolumeOperator =
        std::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
            nonlinearMechanicsBVPoperator->getVolumeOperator() );
    auto mechanicsMaterialModel = nonlinearMechanicsVolumeOperator->getMaterialModel();

    // Create a Linear BVP operator for mechanics
    AMP_INSIST( input_db->keyExists( "LinearMechanicsOperator" ), "key missing!" );
    auto linearMechanicsBVPoperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "LinearMechanicsOperator", input_db, mechanicsMaterialModel ) );

    // Create the variables
    auto mechanicsNonlinearVolumeOperator =
        std::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
            nonlinearMechanicsBVPoperator->getVolumeOperator() );
    auto dispVar = mechanicsNonlinearVolumeOperator->getOutputVariable();

    // For RHS (Point Forces)
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
    auto dirichletLoadVecOp = std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "Load_Boundary", input_db, dummyModel ) );
    dirichletLoadVecOp->setVariable( dispVar );

    auto dofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );

    // Create the vectors
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    auto solVec = AMP::LinearAlgebra::createVector( dofMap, dispVar, true );
    auto rhsVec = solVec->cloneVector();
    auto resVec = solVec->cloneVector();
    // AMP::LinearAlgebra::Vector::shared_ptr scaledRhsVec = meshAdapter->createVector( dispVar );

    // Initial guess
    solVec->zero(*solVec);
    nonlinearMechanicsBVPoperator->modifyInitialSolutionVector( solVec );

    // RHS
    rhsVec->zero(*rhsVec);
    dirichletLoadVecOp->apply( nullVec, rhsVec );
    nonlinearMechanicsBVPoperator->modifyRHSvector( rhsVec );

    // We need to reset the linear operator before the solve since TrilinosML does
    // the factorization of the matrix during construction and so the matrix must
    // be correct before constructing the TrilinosML object.
    nonlinearMechanicsBVPoperator->apply( solVec, resVec );
    linearMechanicsBVPoperator->reset(
        nonlinearMechanicsBVPoperator->getParameters( "Jacobian", solVec ) );

    std::shared_ptr<AMP::LinearAlgebra::Matrix> mechMat = linearMechanicsBVPoperator->getMatrix();

    for ( int i = 0; i < 24; i++ ) {
        std::vector<size_t> matCols;
        std::vector<double> matVals;
        mechMat->getRowByGlobalID( i, matCols, matVals );
        for ( unsigned int j = 0; j < matCols.size(); j++ ) {
            fprintf(
                fp, "A(%d, %d) = %.15f ; \n", ( i + 1 ), (int) ( matCols[j] + 1 ), matVals[j] );
        } // end for j
        fprintf( fp, "\n" );
    } // end for i

    ut->passes( exeName );
}

int testUpdatedLagrangianMechanics_eigenValues_1( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    auto libmeshInit = std::make_shared<AMP::Mesh::initializeLibMesh>( AMP_COMM_WORLD );

    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.emplace_back( "testUpdatedLagrangianMechanics-eigenValues-1" );

    for ( auto &exeName : exeNames )
        myTest( &ut, exeName );

    ut.report();
    int num_failed = ut.NumFailGlobal();

    libmeshInit.reset();
    AMP::AMPManager::shutdown();
    return num_failed;
}
