#include "AMP/IO/PIO.h"
#include "AMP/mesh/MeshFactory.h"
#include "AMP/mesh/MeshParameters.h"
#include "AMP/solvers/libmesh/PelletStackHelpers.h"
#include "AMP/solvers/petsc/PetscSNESSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>


static void myTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    auto global_input_db = AMP::Database::parseInputFile( input_file );
    global_input_db->print( AMP::plog );

    unsigned int NumberOfLoadingSteps = global_input_db->getScalar<int>( "NumberOfLoadingSteps" );
    bool usePointLoad                 = global_input_db->getScalar<bool>( "USE_POINT_LOAD" );
    bool useThermalLoad               = global_input_db->getScalar<bool>( "USE_THERMAL_LOAD" );

    AMP_INSIST( global_input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto mesh_db    = global_input_db->getDatabase( "Mesh" );
    auto meshParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    meshParams->setComm( globalComm );
    auto manager = AMP::Mesh::MeshFactory::create( meshParams );

    std::shared_ptr<AMP::Operator::CoupledOperator> coupledOp;
    std::shared_ptr<AMP::Operator::ColumnOperator> linearColumnOperator;
    std::shared_ptr<AMP::Operator::PelletStackOperator> pelletStackOp;
    helperCreateAllOperatorsForPelletMechanics(
        manager, globalComm, global_input_db, coupledOp, linearColumnOperator, pelletStackOp );

    AMP::LinearAlgebra::Vector::shared_ptr solVec, rhsVec, scaledRhsVec;
    helperCreateVectorsForPelletMechanics( manager, coupledOp, solVec, rhsVec, scaledRhsVec );

    if ( usePointLoad ) {
        helperBuildPointLoadRHSForPelletMechanics( global_input_db, coupledOp, rhsVec );
    } else {
        rhsVec->zero();
    }

    AMP::LinearAlgebra::Vector::shared_ptr initialTemperatureVec, finalTemperatureVec;
    if ( useThermalLoad ) {
        helperCreateTemperatureVectorsForPelletMechanics(
            manager, initialTemperatureVec, finalTemperatureVec );
    }

    if ( useThermalLoad ) {
        auto initialTemp = global_input_db->getScalar<double>( "InitialTemperature" );
        initialTemperatureVec->setToScalar( initialTemp );
        helperSetReferenceTemperatureForPelletMechanics( coupledOp, initialTemperatureVec );
    }

    solVec->zero();
    helperApplyBoundaryCorrectionsForPelletMechanics( coupledOp, solVec, rhsVec );

    auto nonlinearSolver_db   = global_input_db->getDatabase( "NonlinearSolver" );
    auto linearSolver_db      = nonlinearSolver_db->getDatabase( "LinearSolver" );
    auto pelletStackSolver_db = linearSolver_db->getDatabase( "PelletStackSolver" );

    std::shared_ptr<AMP::Solver::SolverStrategy> pelletStackSolver;
    helperBuildStackSolverForPelletMechanics(
        pelletStackSolver_db, pelletStackOp, linearColumnOperator, pelletStackSolver );

    auto nonlinearSolverParams =
        std::make_shared<AMP::Solver::SolverStrategyParameters>( nonlinearSolver_db );
    nonlinearSolverParams->d_comm          = globalComm;
    nonlinearSolverParams->d_pOperator     = coupledOp;
    nonlinearSolverParams->d_pInitialGuess = solVec;
    auto nonlinearSolver = std::make_shared<AMP::Solver::PetscSNESSolver>( nonlinearSolverParams );

    auto linearSolver = nonlinearSolver->getKrylovSolver();
    linearSolver->setNestedSolver( pelletStackSolver );

    for ( unsigned int step = 0; step < NumberOfLoadingSteps; step++ ) {
        AMP::pout << "########################################" << std::endl;
        AMP::pout << "The current loading step is " << ( step + 1 ) << std::endl;

        double scaleValue =
            ( static_cast<double>( step + 1 ) ) / ( static_cast<double>( NumberOfLoadingSteps ) );
        scaledRhsVec->scale( scaleValue, *rhsVec );

        if ( useThermalLoad ) {
            auto initialTemp = global_input_db->getScalar<double>( "InitialTemperature" );
            auto finalTemp   = global_input_db->getScalar<double>( "FinalTemperature" );
            double deltaTemp =
                initialTemp + ( ( static_cast<double>( step + 1 ) ) * ( finalTemp - initialTemp ) /
                                ( static_cast<double>( NumberOfLoadingSteps ) ) );
            finalTemperatureVec->setToScalar( deltaTemp );
            helperSetFinalTemperatureForPelletMechanics( coupledOp, finalTemperatureVec );
        }

        auto resVec = solVec->clone();
        resVec->zero();
        coupledOp->residual( scaledRhsVec, solVec, resVec );
        AMP::pout << "initial, rhsVec: " << scaledRhsVec->L2Norm() << std::endl;
        AMP::pout << "initial, solVec: " << solVec->L2Norm() << std::endl;
        AMP::pout << "initial, resVec: " << resVec->L2Norm() << std::endl;
        nonlinearSolver->apply( scaledRhsVec, solVec );
        AMP::pout << "solved,  rhsVec: " << scaledRhsVec->L2Norm() << std::endl;
        AMP::pout << "solved,  solVec: " << solVec->L2Norm() << std::endl;
        coupledOp->residual( scaledRhsVec, solVec, resVec );
        AMP::pout << "final,   rhsVec: " << scaledRhsVec->L2Norm() << std::endl;
        AMP::pout << "final,   solVec: " << solVec->L2Norm() << std::endl;
        AMP::pout << "final,   resVec: " << resVec->L2Norm() << std::endl;

        helperResetNonlinearOperatorForPelletMechanics( coupledOp );
    } // end for step

    ut->passes( exeName );
}

int testPelletStackMechanicsSolver( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    const char *exeName[] = { "testPelletStackMechanicsSolver-1",
                              "testPelletStackMechanicsSolver-2" };
    for ( auto exe : exeName )
        myTest( &ut, exe );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
