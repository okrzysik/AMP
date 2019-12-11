#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/CustomConstraintsEliminationOperator.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletMatrixCorrection.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/mechanics/IsotropicElasticModel.h"
#include "AMP/operators/mechanics/MechanicsLinearElement.h"
#include "AMP/operators/mechanics/MechanicsLinearFEOperator.h"
#include "AMP/operators/petsc/PetscMatrixShellOperator.h"
#include "AMP/solvers/ColumnSolver.h"
#include "AMP/solvers/ConstraintsEliminationSolver.h"
#include "AMP/solvers/petsc/PetscKrylovSolver.h"
#include "AMP/solvers/petsc/PetscKrylovSolverParameters.h"
#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/VectorBuilder.h"
#include <memory>

#include <iostream>
#include <string>


static void myTest( AMP::UnitTest *ut )
{
    std::string exeName( "testConstraintsEliminationOperator" );
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );


    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    std::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    std::shared_ptr<AMP::Mesh::MeshParameters> meshParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    meshParams->setComm( globalComm );
    AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh( meshParams );

    int const dofsPerNode = 3;
    int const gostWidth   = 1;
    bool const split      = true;
    AMP::Discretization::DOFManager::shared_ptr dofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, gostWidth, dofsPerNode, split );

    AMP::LinearAlgebra::Vector::shared_ptr vec1;
    AMP::LinearAlgebra::Vector::shared_ptr vec2;

    for ( int dummy = 0; dummy < 2; ++dummy ) {
        std::shared_ptr<AMP::Operator::OperatorParameters> emptyParams;
        std::shared_ptr<AMP::Operator::ColumnOperator> colOp(
            new AMP::Operator::ColumnOperator( emptyParams ) );
        std::shared_ptr<AMP::Database> linearSolver_db = input_db->getDatabase( "LinearSolver" );
        std::shared_ptr<AMP::Database> preconditioner_db =
            linearSolver_db->getDatabase( "Preconditioner" );
        std::shared_ptr<AMP::Solver::ColumnSolverParameters> preconditionerParams(
            new AMP::Solver::ColumnSolverParameters( preconditioner_db ) );
        preconditionerParams->d_pOperator = colOp;
        std::shared_ptr<AMP::Solver::ColumnSolver> colPre(
            new AMP::Solver::ColumnSolver( preconditionerParams ) );

        std::shared_ptr<AMP::Operator::ElementPhysicsModel> physicsModel;
        std::shared_ptr<AMP::Operator::LinearBVPOperator> bvpOp =
            std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
                AMP::Operator::OperatorBuilder::createOperator(
                    meshAdapter,
                    ( dummy ? "DummyBVPOperator" : "BVPOperator" ),
                    input_db,
                    physicsModel ) );
        colOp->append( bvpOp );

        std::shared_ptr<AMP::Database> solver_db = preconditioner_db->getDatabase( "Solver" );
        std::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> solverParams(
            new AMP::Solver::PetscKrylovSolverParameters( solver_db ) );
        solverParams->d_pOperator = bvpOp;
        solverParams->d_comm      = globalComm;
        std::shared_ptr<AMP::Solver::PetscKrylovSolver> solver(
            new AMP::Solver::PetscKrylovSolver( solverParams ) );
        colPre->append( solver );

        std::shared_ptr<AMP::Operator::CustomConstraintsEliminationOperator> dirOp;
        if ( dummy ) {
            std::shared_ptr<AMP::Database> dummyOperator_db( new AMP::Database( "DummyOperator" ) );
            dummyOperator_db->putScalar( "InputVariable", "displacement" );
            dummyOperator_db->putScalar( "OutputVariable", "displacement" );
            std::shared_ptr<AMP::Operator::OperatorParameters> dummyOperatorParams(
                new AMP::Operator::OperatorParameters( dummyOperator_db ) );
            dirOp = std::make_shared<AMP::Operator::CustomConstraintsEliminationOperator>(
                dummyOperatorParams );
            std::vector<size_t> slaveIndices;
            std::vector<double> slaveValues;
            int const boundaryID = 2;
            AMP::Mesh::MeshIterator it =
                meshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, boundaryID );
            size_t const numDOFs = it.size() * dofsPerNode;
            slaveIndices.resize( numDOFs );
            slaveValues.resize( numDOFs );
            AMP::Mesh::MeshIterator it_begin = it.begin();
            AMP::Mesh::MeshIterator it_end   = it.end();
            std::vector<size_t> dofIndices;
            double dirichletValues[3] = { 2.0, 1.0, 0.0 };
            size_t p                  = 0;
            for ( it = it_begin; it != it_end; ++it ) {
                dofMap->getDOFs( it->globalID(), dofIndices );
                AMP_ASSERT( static_cast<int>( dofIndices.size() ) == dofsPerNode );
                std::copy( &( dofIndices[0] ),
                           &( dofIndices[0] ) + dofsPerNode,
                           &( slaveIndices[0] ) + p );
                std::copy( &( dirichletValues[0] ),
                           &( dirichletValues[0] ) + dofsPerNode,
                           &( slaveValues[0] ) + p );
                p += dofsPerNode;
            } // end for
            AMP_ASSERT( p == numDOFs );
            dirOp->initialize( slaveIndices, slaveValues );
            //    colOp->append(dirOp);

            //    std::shared_ptr<AMP::Database> dummySolver_db =
            //    preconditioner_db->getDatabase("ContactPreconditioner");
            std::shared_ptr<AMP::Database> dummySolver_db( new AMP::Database( "DummySolver" ) );
            dummySolver_db->putScalar( "print_info_level", 1 );
            dummySolver_db->putScalar( "max_iterations", 1 );
            dummySolver_db->putScalar( "max_error", 1.0e-16 );
            std::shared_ptr<AMP::Solver::ConstraintsEliminationSolverParameters> dummySolverParams(
                new AMP::Solver::ConstraintsEliminationSolverParameters( dummySolver_db ) );
            dummySolverParams->d_pOperator = dirOp;
            std::shared_ptr<AMP::Solver::ConstraintsEliminationSolver> dirSolver(
                new AMP::Solver::ConstraintsEliminationSolver( dummySolverParams ) );
            colPre->append( dirSolver );
        } // end if

        std::shared_ptr<AMP::Operator::DirichletVectorCorrection> loadOp =
            std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
                AMP::Operator::OperatorBuilder::createOperator(
                    meshAdapter, "LoadOperator", input_db, physicsModel ) );
        AMP::LinearAlgebra::Variable::shared_ptr var = bvpOp->getOutputVariable();
        loadOp->setVariable( var );

        std::shared_ptr<AMP::Database> shell_db( new AMP::Database( "MatrixShellOperator" ) );
        shell_db->putScalar( "name", "MatShellOperator" );
        shell_db->putScalar( "print_info_level", 1 );
        std::shared_ptr<AMP::Operator::OperatorParameters> shellParams(
            new AMP::Operator::OperatorParameters( shell_db ) );
        std::shared_ptr<AMP::Operator::PetscMatrixShellOperator> shellOp(
            new AMP::Operator::PetscMatrixShellOperator( shellParams ) );
        int const numLocalNodes = meshAdapter->numLocalElements( AMP::Mesh::GeomType::Vertex );
        int const matLocalSize  = dofsPerNode * numLocalNodes;
        AMP_ASSERT( matLocalSize == static_cast<int>( dofMap->numLocalDOF() ) );
        shellOp->setComm( globalComm );
        shellOp->setMatLocalRowSize( matLocalSize );
        shellOp->setMatLocalColumnSize( matLocalSize );
        shellOp->setOperator( colOp );

        AMP::LinearAlgebra::Vector::shared_ptr nullVec;
        AMP::LinearAlgebra::Vector::shared_ptr solVec =
            AMP::LinearAlgebra::createVector( dofMap, var, split );
        AMP::LinearAlgebra::Vector::shared_ptr rhsVec = solVec->cloneVector();
        AMP::LinearAlgebra::Vector::shared_ptr resVec = solVec->cloneVector();
        AMP::LinearAlgebra::Vector::shared_ptr dirVec = solVec->cloneVector();
        AMP::LinearAlgebra::Vector::shared_ptr corVec = solVec->cloneVector();

        solVec->zero();
        rhsVec->zero();
        loadOp->apply( nullVec, rhsVec );
        bvpOp->modifyRHSvector( rhsVec );
        if ( dummy ) {
            dirVec->zero();
            dirOp->addShiftToSlave( dirVec );
            colOp->apply( dirVec, corVec );
            corVec->scale( -1.0 );
            colOp->append( dirOp );
            rhsVec->add( rhsVec, corVec );
            dirOp->addSlaveToMaster( rhsVec );
            dirOp->setSlaveToZero( rhsVec );
            dirOp->copyMasterToSlave( solVec );
        } // end if

        //  shellOp->apply(rhsVec, solVec, resVec, -1.0, 1.0);
        //  double const initialResidualNorm = resVec->L2Norm();
        //  std::cout<<"initial residual norm = "<<std::setprecision(15)<<initialResidualNorm<<"\n";

        std::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(
            new AMP::Solver::PetscKrylovSolverParameters( linearSolver_db ) );
        linearSolverParams->d_pOperator       = shellOp;
        linearSolverParams->d_comm            = globalComm;
        linearSolverParams->d_pPreconditioner = colPre;
        std::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(
            new AMP::Solver::PetscKrylovSolver( linearSolverParams ) );
        linearSolver->setInitialGuess( solVec );

        linearSolver->solve( rhsVec, solVec );

        if ( dummy ) {
            dirOp->copyMasterToSlave( solVec );
            dirOp->addShiftToSlave( solVec );
        } // end if

        //  shellOp->apply(rhsVec, solVec, resVec, -1.0, 1.0);
        //  double const finalResidualNorm = resVec->L2Norm();
        //  std::cout<<"final residual norm = "<<std::setprecision(15)<<finalResidualNorm<<"\n";
        if ( dummy ) {
            vec1 = solVec;
        } else {
            vec2 = solVec;
        } // end i

    } // end for

    vec1->subtract( vec1, vec2 );
    double const solutionL2Norm = vec2->L2Norm();
    double const errorL2Norm    = vec1->L2Norm();
    std::cout << "solution L2 norm = " << std::setprecision( 15 ) << solutionL2Norm << "\n";
    std::cout << "error L2 norm = " << std::setprecision( 15 ) << errorL2Norm << "\n";
    std::cout << "relative error = " << errorL2Norm / solutionL2Norm << "\n";
    AMP_ASSERT( errorL2Norm / solutionL2Norm < 1.0e-12 );

    ut->passes( exeName );
}

int testConstraintsEliminationOperator( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    myTest( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
