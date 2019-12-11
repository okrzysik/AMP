#include "AMP/ampmesh/libmesh/initializeLibMesh.h"
#include "AMP/ampmesh/libmesh/libMesh.h"
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
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include <memory>

DISABLE_WARNINGS
#include "libmesh/mesh_communication.h"
#include "petsc.h"
#include "petscksp.h"
ENABLE_WARNINGS

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

    AMP::PIO::logOnlyNodeZero( log_file );


    AMP::AMP_MPI globalComm = AMP::AMP_MPI( AMP_COMM_WORLD );
    auto input_db           = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    int numMeshes = input_db->getScalar<int>( "NumberOfMeshFiles" );
    std::shared_ptr<AMP::Mesh::initializeLibMesh> libmeshInit(
        new AMP::Mesh::initializeLibMesh( globalComm ) );

    std::shared_ptr<AMP::Database> ml_db   = input_db->getDatabase( "ML_Solver" );
    std::shared_ptr<AMP::Database> lu_db   = input_db->getDatabase( "LU_Solver" );
    std::shared_ptr<AMP::Database> cg_db   = input_db->getDatabase( "CG_Solver" );
    std::shared_ptr<AMP::Database> rich_db = input_db->getDatabase( "Richardson_Solver" );

    for ( int meshId = 1; meshId <= numMeshes; meshId++ ) {
        std::cout << "Working on mesh " << meshId << std::endl;

        char meshFileKey[200];
        sprintf( meshFileKey, "mesh%d", meshId );

        std::string meshFile = input_db->getString( meshFileKey );

        const unsigned int mesh_dim = 3;
        std::shared_ptr<::Mesh> mesh( new ::Mesh( mesh_dim ) );

        if ( globalComm.getRank() == 0 ) {
            AMP::readBinaryTestMesh( meshFile, mesh );
        }

        MeshCommunication().broadcast( *( mesh.get() ) );
        mesh->prepare_for_use( false );

        AMP::Mesh::Mesh::shared_ptr meshAdapter( new AMP::Mesh::libMesh( mesh, "mesh" ) );

        std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
        std::shared_ptr<AMP::Operator::LinearBVPOperator> bvpOperator =
            std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
                AMP::Operator::OperatorBuilder::createOperator(
                    meshAdapter, "BVPOperator", input_db, elementPhysicsModel ) );

        AMP::LinearAlgebra::Variable::shared_ptr dispVar = bvpOperator->getOutputVariable();

        std::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
        std::shared_ptr<AMP::Operator::DirichletVectorCorrection> loadOperator =
            std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
                AMP::Operator::OperatorBuilder::createOperator(
                    meshAdapter, "LoadOperator", input_db, dummyModel ) );
        loadOperator->setVariable( dispVar );

        AMP::Discretization::DOFManager::shared_ptr NodalVectorDOF =
            AMP::Discretization::simpleDOFManager::create(
                meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3 );

        AMP::LinearAlgebra::Vector::shared_ptr nullVec;
        AMP::LinearAlgebra::Vector::shared_ptr solVec =
            AMP::LinearAlgebra::createVector( NodalVectorDOF, dispVar );
        AMP::LinearAlgebra::Vector::shared_ptr rhsVec =
            AMP::LinearAlgebra::createVector( NodalVectorDOF, dispVar );
        AMP::LinearAlgebra::Vector::shared_ptr resVec =
            AMP::LinearAlgebra::createVector( NodalVectorDOF, dispVar );

        rhsVec->zero();
        loadOperator->apply( nullVec, rhsVec );

        solVec->zero();
        resVec->zero();

        size_t numDofs = solVec->getGlobalSize();

        if ( globalComm.getRank() == 0 ) {
            std::cout << "Solving using LU" << std::endl;
        }
        globalComm.barrier();
        double luStartTime = AMP::AMP_MPI::time();

        std::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> luParams(
            new AMP::Solver::TrilinosMLSolverParameters( lu_db ) );
        luParams->d_pOperator = bvpOperator;
        std::shared_ptr<AMP::Solver::TrilinosMLSolver> luPC(
            new AMP::Solver::TrilinosMLSolver( luParams ) );

        std::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> richParams(
            new AMP::Solver::PetscKrylovSolverParameters( rich_db ) );
        richParams->d_pOperator       = bvpOperator;
        richParams->d_comm            = globalComm;
        richParams->d_pPreconditioner = luPC;
        std::shared_ptr<AMP::Solver::PetscKrylovSolver> richSolver(
            new AMP::Solver::PetscKrylovSolver( richParams ) );
        richSolver->setZeroInitialGuess( true );

        richSolver->solve( rhsVec, solVec );

        globalComm.barrier();
        double luEndTime = AMP::AMP_MPI::time();

        PetscInt richIters = 0;
        KSP richKsp        = richSolver->getKrylovSolver();
        KSPGetIterationNumber( richKsp, &richIters );
        AMP_INSIST( richIters <= 1, "Should not need more than 1 LU-Richardson iteration." );

        KSPConvergedReason richReason;
        KSPGetConvergedReason( richKsp, &richReason );
        AMP_INSIST(
            ( ( richReason == KSP_CONVERGED_RTOL ) || ( richReason == KSP_CONVERGED_ATOL ) ),
            "KSP did not converge properly." );

        richSolver.reset();
        luPC.reset();

        solVec->zero();
        resVec->zero();

        if ( globalComm.getRank() == 0 ) {
            std::cout << "Solving using ML" << std::endl;
        }
        globalComm.barrier();
        double mlStartTime = AMP::AMP_MPI::time();

        std::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> mlParams(
            new AMP::Solver::TrilinosMLSolverParameters( ml_db ) );
        mlParams->d_pOperator = bvpOperator;
        std::shared_ptr<AMP::Solver::TrilinosMLSolver> mlPC(
            new AMP::Solver::TrilinosMLSolver( mlParams ) );

        std::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> cgParams(
            new AMP::Solver::PetscKrylovSolverParameters( cg_db ) );
        cgParams->d_pOperator       = bvpOperator;
        cgParams->d_comm            = globalComm;
        cgParams->d_pPreconditioner = mlPC;
        std::shared_ptr<AMP::Solver::PetscKrylovSolver> cgSolver(
            new AMP::Solver::PetscKrylovSolver( cgParams ) );
        cgSolver->setZeroInitialGuess( true );

        cgSolver->solve( rhsVec, solVec );

        globalComm.barrier();
        double mlEndTime = AMP::AMP_MPI::time();

        PetscInt cgIters = 0;
        KSP cgKsp        = cgSolver->getKrylovSolver();
        KSPGetIterationNumber( cgKsp, &cgIters );

        KSPConvergedReason cgReason;
        KSPGetConvergedReason( cgKsp, &cgReason );
        AMP_INSIST( ( ( cgReason == KSP_CONVERGED_RTOL ) || ( cgReason == KSP_CONVERGED_ATOL ) ),
                    "KSP did not converge properly." );

        cgSolver.reset();
        mlPC.reset();

        if ( globalComm.getRank() == 0 ) {
            std::cout << "Result: " << numDofs << " & " << globalComm.getSize() << " & " << cgIters
                      << " & " << ( luEndTime - luStartTime ) << " & "
                      << ( mlEndTime - mlStartTime ) << " \\\\ " << std::endl;
        }

    } // end for meshId

    ut->passes( exeName );
}

int testMechElasticLUvsML( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::string exeName = "testMechElasticLUvsML";

    myTest( &ut, exeName );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
