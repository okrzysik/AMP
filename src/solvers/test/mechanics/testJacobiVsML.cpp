#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/solvers/petsc/PetscKrylovSolver.h"
#include "AMP/solvers/petsc/PetscKrylovSolverParameters.h"
#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/WriteSolutionToFile.h"
#include "AMP/vectors/VectorBuilder.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>


static void myTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    // Read the input file
    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Read the mesh
    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    std::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    std::shared_ptr<AMP::Mesh::MeshParameters> meshParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh( meshParams );

    std::cout << "Mesh has " << ( meshAdapter->numLocalElements( AMP::Mesh::GeomType::Vertex ) )
              << " nodes." << std::endl;

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
    std::shared_ptr<AMP::Operator::LinearBVPOperator> bvpOperator =
        std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "MechanicsBVPOperator", input_db, elementPhysicsModel ) );

    AMP::LinearAlgebra::Variable::shared_ptr displacementVariable =
        bvpOperator->getOutputVariable();

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
    std::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletVecOp =
        std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "Load_Boundary", input_db, dummyModel ) );
    // This has an in-place apply. So, it has an empty input variable and
    // the output variable is the same as what it is operating on.
    dirichletVecOp->setVariable( displacementVariable );

    AMP::Discretization::DOFManager::shared_ptr dofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    AMP::LinearAlgebra::Vector::shared_ptr mechSolVec =
        AMP::LinearAlgebra::createVector( dofMap, displacementVariable, true );
    AMP::LinearAlgebra::Vector::shared_ptr mechRhsVec = mechSolVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr mechResVec = mechSolVec->cloneVector();

    mechRhsVec->zero(mechRhsVec);
    mechResVec->zero(mechResVec);

    dirichletVecOp->apply( nullVec, mechRhsVec );

    for ( int type = 1; type < 4; type++ ) {
        if ( type == 0 ) {
            std::cout << "Solving using CG algorithm (Own Implementation)..." << std::endl;

            std::shared_ptr<AMP::Database> linearSolver_db = input_db->getDatabase( "CGsolver" );

            int maxIters = linearSolver_db->getScalar<int>( "max_iterations" );

            AMP::LinearAlgebra::Vector::shared_ptr matOutVec = mechSolVec->cloneVector();
            AMP::LinearAlgebra::Vector::shared_ptr pVec      = mechSolVec->cloneVector();

            mechSolVec->zero(mechSolVec);

            bvpOperator->apply( mechSolVec, matOutVec );

            mechResVec->subtract( mechRhsVec, matOutVec, mechResVec );

            pVec->copyVector( mechResVec );

            for ( int iter = 0; iter <= maxIters; iter++ ) {
                double resNorm = mechResVec->L2Norm(mechResVec);
                std::cout << "Iter = " << iter << " ResNorm2 = " << std::setprecision( 15 )
                          << resNorm << std::endl;

                bvpOperator->apply( pVec, matOutVec );

                double matOutNorm = matOutVec->L2Norm(matOutVec);
                std::cout << "CG-Iter = " << iter << " MatOutNorm2 = " << std::setprecision( 15 )
                          << matOutNorm << std::endl;

                double resOldDot = mechResVec->dot( mechResVec, mechResVec );

                double alphaDenom = matOutVec->dot( pVec, matOutVec );

                double alpha = resOldDot / alphaDenom;

                mechSolVec->axpy( alpha, pVec, mechSolVec, mechSolVec );

                mechResVec->axpy( -alpha, matOutVec, mechResVec, mechResVec );

                double resNewDot = mechResVec->dot( mechResVec, mechResVec );

                double beta = resNewDot / resOldDot;

                std::cout << "Iter = " << iter << " resOldDot = " << std::setprecision( 15 )
                          << resOldDot << " alphaDenom = " << std::setprecision( 15 ) << alphaDenom
                          << " alpha = " << std::setprecision( 15 ) << alpha
                          << " resNewDot = " << std::setprecision( 15 ) << resNewDot
                          << " beta = " << std::setprecision( 15 ) << beta << std::endl
                          << std::endl;

                pVec->axpy( beta, pVec, mechResVec, pVec );
            }

            std::cout << std::endl << std::endl;
        } else if ( type == 1 ) {
            std::cout << "Solving using CG algorithm (Petsc Implementation)..." << std::endl;

            std::shared_ptr<AMP::Database> linearSolver_db = input_db->getDatabase( "CGsolver" );

            // initialize the linear solver
            std::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(
                new AMP::Solver::PetscKrylovSolverParameters( linearSolver_db ) );
            linearSolverParams->d_pOperator = bvpOperator;
            linearSolverParams->d_comm      = globalComm;
            std::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(
                new AMP::Solver::PetscKrylovSolver( linearSolverParams ) );

            linearSolver->solve( mechRhsVec, mechSolVec );

            std::cout << std::endl << std::endl;
        } else if ( type == 2 ) {
            std::cout << "Solving using Jacobi preconditioned CG algorithm..." << std::endl;

            std::shared_ptr<AMP::Database> linearSolver_db =
                input_db->getDatabase( "JacobiCGsolver" );

            // initialize the linear solver
            std::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(
                new AMP::Solver::PetscKrylovSolverParameters( linearSolver_db ) );
            linearSolverParams->d_pOperator = bvpOperator;
            linearSolverParams->d_comm      = globalComm;
            std::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(
                new AMP::Solver::PetscKrylovSolver( linearSolverParams ) );

            linearSolver->solve( mechRhsVec, mechSolVec );

            std::cout << std::endl << std::endl;
        } else {
            std::cout << "Solving using ML preconditioned CG algorithm..." << std::endl;

            std::shared_ptr<AMP::Database> linearSolver_db = input_db->getDatabase( "MLCGsolver" );

            // ---- first initialize the preconditioner
            std::shared_ptr<AMP::Database> pcSolver_db = linearSolver_db->getDatabase( "MLsolver" );
            std::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> pcSolverParams(
                new AMP::Solver::TrilinosMLSolverParameters( pcSolver_db ) );
            pcSolverParams->d_pOperator = bvpOperator;
            std::shared_ptr<AMP::Solver::TrilinosMLSolver> pcSolver(
                new AMP::Solver::TrilinosMLSolver( pcSolverParams ) );

            // initialize the linear solver
            std::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(
                new AMP::Solver::PetscKrylovSolverParameters( linearSolver_db ) );
            linearSolverParams->d_pOperator       = bvpOperator;
            linearSolverParams->d_comm            = globalComm;
            linearSolverParams->d_pPreconditioner = pcSolver;
            std::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(
                new AMP::Solver::PetscKrylovSolver( linearSolverParams ) );

            linearSolver->solve( mechRhsVec, mechSolVec );

            std::cout << std::endl << std::endl;
        }
    }

    ut->passes( exeName );
}


int testJacobiVsML( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::string exeName = "testJacobiVsML";

    myTest( &ut, exeName );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
