
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"

#include <cstdlib>
#include <iostream>
#include <string>

#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/libmesh/libMesh.h"

#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/vectors/VectorBuilder.h"

#include "libmesh/libmesh.h"
#include "libmesh/mesh_communication.h"

#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/mechanics/MechanicsNonlinearFEOperator.h"

#include "AMP/utils/ReadTestMesh.h"

static void myTest( AMP::UnitTest *ut, const std::string &exeName, int callLinReset )
{
    auto input_file = "input_" + exeName;
    auto log_file   = "output_" + exeName;
    auto msgName    = exeName;
    if ( callLinReset ) {
        log_file = log_file + "-1";
        msgName  = msgName + "-1";
    } else {
        log_file = log_file + "-0";
        msgName  = msgName + "-0";
    }

    AMP::PIO::logOnlyNodeZero( log_file );


    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    std::string mesh_file = input_db->getString( "mesh_file" );

    const unsigned int mesh_dim = 3;
    std::shared_ptr<::Mesh> mesh( new ::Mesh( mesh_dim ) );

    if ( ut->rank() == 0 ) {
        AMP::readTestMesh( mesh_file, mesh );
    } // end if root processor

    MeshCommunication().broadcast( *( mesh.get() ) );

    mesh->prepare_for_use( false );

    auto meshAdapter = std::make_shared<AMP::Mesh::libMesh>( mesh, "TestMesh" );

    auto nonlinOperator = std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "NonlinearMechanicsOperator", input_db ) );

    auto mechNonlinOperator =
        std::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
            nonlinOperator->getVolumeOperator() );
    auto elementPhysicsModel = mechNonlinOperator->getMaterialModel();

    auto linOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "LinearMechanicsOperator", input_db, elementPhysicsModel ) );

    auto dofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );

    auto var = nonlinOperator->getOutputVariable();

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    auto solVec       = AMP::LinearAlgebra::createVector( dofMap, var, true );
    auto resVecNonlin = solVec->cloneVector();
    auto resVecLin    = solVec->cloneVector();
    auto resDiffVec   = solVec->cloneVector();
    auto tmpNonlinVec = solVec->cloneVector();
    auto tmpLinVec    = solVec->cloneVector();

    solVec->setToScalar( 0.0 );
    double solNorm = solVec->L2Norm();
    AMP::pout << "sol-Norm-2 = " << solNorm << std::endl;

    nonlinOperator->apply( solVec, resVecNonlin );
    linOperator->reset( nonlinOperator->getParameters( "Jacobian", solVec ) );
    linOperator->apply( solVec, resVecLin );
    resDiffVec->subtract( resVecNonlin, resVecLin );

    double epsilon = 1.0e-13 * ( ( ( linOperator->getMatrix() )->extractDiagonal() )->L1Norm() );
    AMP::pout << "epsilon = " << epsilon << std::endl;

    double nonLinNorm = resVecNonlin->L1Norm();
    double linNorm    = resVecLin->L1Norm();
    double diffNorm   = resDiffVec->L1Norm();

    AMP::pout << "nonLin-Norm-1 = " << nonLinNorm << std::endl;
    AMP::pout << "lin-Norm-1 = " << linNorm << std::endl;
    AMP::pout << "diff-Norm-1 = " << diffNorm << std::endl;

    if ( nonLinNorm > epsilon ) {
        ut->failure( msgName );
    }

    if ( linNorm > epsilon ) {
        ut->failure( msgName );
    }

    if ( diffNorm > epsilon ) {
        ut->failure( msgName );
    }

    solVec->setRandomValues();
    solVec->scale( 100.0 );
    nonlinOperator->modifyInitialSolutionVector( solVec );
    solNorm = solVec->L2Norm();
    AMP::pout << "sol-Norm-2 = " << solNorm << std::endl;

    nonlinOperator->apply( solVec, resVecNonlin );
    if ( callLinReset ) {
        linOperator->reset( nonlinOperator->getParameters( "Jacobian", solVec ) );
    }
    linOperator->apply( solVec, resVecLin );
    resDiffVec->subtract( resVecNonlin, resVecLin );

    nonLinNorm = resVecNonlin->L2Norm();
    linNorm    = resVecLin->L2Norm();
    diffNorm   = resDiffVec->L1Norm();

    AMP::pout << "nonLin-Norm-2 = " << nonLinNorm << std::endl;
    AMP::pout << "lin-Norm-2 = " << linNorm << std::endl;
    AMP::pout << "diff-Norm-1 = " << diffNorm << std::endl;

    if ( diffNorm > epsilon ) {
        ut->failure( msgName );
    }

    tmpNonlinVec->copyVector( resVecNonlin );
    tmpNonlinVec->scale( 10.0 );

    tmpLinVec->copyVector( resVecLin );
    tmpLinVec->scale( 10.0 );

    solVec->scale( 10.0 );
    nonlinOperator->modifyInitialSolutionVector( solVec );
    solNorm = solVec->L2Norm();
    AMP::pout << "sol-Norm-2 = " << solNorm << std::endl;

    nonlinOperator->apply( solVec, resVecNonlin );
    if ( callLinReset ) {
        linOperator->reset( nonlinOperator->getParameters( "Jacobian", solVec ) );
    }
    linOperator->apply( solVec, resVecLin );

    nonLinNorm = resVecNonlin->L2Norm();
    linNorm    = resVecLin->L2Norm();
    AMP::pout << "nonLin-Norm-2 = " << nonLinNorm << std::endl;
    AMP::pout << "lin-Norm-2 = " << linNorm << std::endl;

    resDiffVec->subtract( resVecNonlin, tmpNonlinVec );
    diffNorm = resDiffVec->L1Norm();

    if ( diffNorm > ( 10.0 * epsilon ) ) {
        ut->failure( msgName );
    }

    resDiffVec->subtract( resVecLin, tmpLinVec );
    diffNorm = resDiffVec->L1Norm();

    if ( diffNorm > ( 10.0 * epsilon ) ) {
        ut->failure( msgName );
    }

    ut->passes( msgName );
}


int testConsistentTangentBVP( int argc, char *argv[] )
{

    AMP::AMPManager::startup( argc, argv );
    auto libmeshInit = std::make_shared<AMP::Mesh::initializeLibMesh>( AMP_COMM_WORLD );

    AMP::UnitTest ut;

    std::vector<std::string> exeNames;

    exeNames.emplace_back( "testConsistentTangentBVP-1-mesh1-normal" );
    exeNames.emplace_back( "testConsistentTangentBVP-2-mesh1-normal" );
    exeNames.emplace_back( "testConsistentTangentBVP-3-mesh1-normal" );
    exeNames.emplace_back( "testConsistentTangentBVP-4-mesh1-normal" );

    exeNames.emplace_back( "testConsistentTangentBVP-4-mesh1-reduced" );

    for ( int j = 0; j < 2; j++ ) {
        for ( auto &exeName : exeNames )
            myTest( &ut, exeName, j );
    }


    ut.report();
    int num_failed = ut.NumFailGlobal();

    libmeshInit.reset();
    AMP::AMPManager::shutdown();
    return num_failed;
}
