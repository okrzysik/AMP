
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>

#include "utils/shared_ptr.h"

#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"

#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"

#include "operators/ColumnOperator.h"
#include "operators/CustomConstraintsEliminationOperator.h"
#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/boundary/DirichletMatrixCorrection.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/mechanics/IsotropicElasticModel.h"
#include "operators/mechanics/MechanicsLinearElement.h"
#include "operators/mechanics/MechanicsLinearFEOperator.h"
#include "operators/petsc/PetscMatrixShellOperator.h"


#include "solvers/ColumnSolver.h"
#include "solvers/ConstraintsEliminationSolver.h"
#include "solvers/petsc/PetscKrylovSolver.h"
#include "solvers/petsc/PetscKrylovSolverParameters.h"
#include "solvers/trilinos/TrilinosMLSolver.h"

void myTest( AMP::UnitTest *ut )
{
    std::string exeName( "testConstraintsEliminationOperator" );
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> meshParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    meshParams->setComm( globalComm );
    AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh( meshParams );

    int const dofsPerNode = 3;
    int const gostWidth   = 1;
    bool const split      = true;
    AMP::Discretization::DOFManager::shared_ptr dofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::Vertex, gostWidth, dofsPerNode, split );

    AMP::LinearAlgebra::Vector::shared_ptr vec1;
    AMP::LinearAlgebra::Vector::shared_ptr vec2;

    for ( int dummy = 0; dummy < 2; ++dummy ) {
        AMP::shared_ptr<AMP::Operator::OperatorParameters> emptyParams;
        AMP::shared_ptr<AMP::Operator::ColumnOperator> colOp(
            new AMP::Operator::ColumnOperator( emptyParams ) );
        AMP::shared_ptr<AMP::Database> linearSolver_db = input_db->getDatabase( "LinearSolver" );
        AMP::shared_ptr<AMP::Database> preconditioner_db =
            linearSolver_db->getDatabase( "Preconditioner" );
        AMP::shared_ptr<AMP::Solver::ColumnSolverParameters> preconditionerParams(
            new AMP::Solver::ColumnSolverParameters( preconditioner_db ) );
        preconditionerParams->d_pOperator = colOp;
        AMP::shared_ptr<AMP::Solver::ColumnSolver> colPre(
            new AMP::Solver::ColumnSolver( preconditionerParams ) );

        AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> physicsModel;
        AMP::shared_ptr<AMP::Operator::LinearBVPOperator> bvpOp =
            AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
                AMP::Operator::OperatorBuilder::createOperator(
                    meshAdapter,
                    ( dummy ? "DummyBVPOperator" : "BVPOperator" ),
                    input_db,
                    physicsModel ) );
        colOp->append( bvpOp );

        AMP::shared_ptr<AMP::Database> solver_db = preconditioner_db->getDatabase( "Solver" );
        AMP::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> solverParams(
            new AMP::Solver::PetscKrylovSolverParameters( solver_db ) );
        solverParams->d_pOperator = bvpOp;
        solverParams->d_comm      = globalComm;
        AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> solver(
            new AMP::Solver::PetscKrylovSolver( solverParams ) );
        colPre->append( solver );

        AMP::shared_ptr<AMP::Operator::CustomConstraintsEliminationOperator> dirOp;
        if ( dummy ) {
            AMP::shared_ptr<AMP::InputDatabase> dummyOperator_db(
                new AMP::InputDatabase( "DummyOperator" ) );
            dummyOperator_db->putString( "InputVariable", "displacement" );
            dummyOperator_db->putString( "OutputVariable", "displacement" );
            AMP::shared_ptr<AMP::Operator::OperatorParameters> dummyOperatorParams(
                new AMP::Operator::OperatorParameters( dummyOperator_db ) );
            dirOp = AMP::shared_ptr<AMP::Operator::CustomConstraintsEliminationOperator>(
                new AMP::Operator::CustomConstraintsEliminationOperator( dummyOperatorParams ) );
            std::vector<size_t> slaveIndices;
            std::vector<double> slaveValues;
            int const boundaryID = 2;
            AMP::Mesh::MeshIterator it =
                meshAdapter->getBoundaryIDIterator( AMP::Mesh::Vertex, boundaryID );
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

            //    AMP::shared_ptr<AMP::Database> dummySolver_db =
            //    preconditioner_db->getDatabase("ContactPreconditioner");
            AMP::shared_ptr<AMP::InputDatabase> dummySolver_db(
                new AMP::InputDatabase( "DummySolver" ) );
            dummySolver_db->putInteger( "print_info_level", 1 );
            dummySolver_db->putInteger( "max_iterations", 1 );
            dummySolver_db->putDouble( "max_error", 1.0e-16 );
            AMP::shared_ptr<AMP::Solver::ConstraintsEliminationSolverParameters> dummySolverParams(
                new AMP::Solver::ConstraintsEliminationSolverParameters( dummySolver_db ) );
            dummySolverParams->d_pOperator = dirOp;
            AMP::shared_ptr<AMP::Solver::ConstraintsEliminationSolver> dirSolver(
                new AMP::Solver::ConstraintsEliminationSolver( dummySolverParams ) );
            colPre->append( dirSolver );
        } // end if

        AMP::shared_ptr<AMP::Operator::DirichletVectorCorrection> loadOp =
            AMP::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
                AMP::Operator::OperatorBuilder::createOperator(
                    meshAdapter, "LoadOperator", input_db, physicsModel ) );
        AMP::LinearAlgebra::Variable::shared_ptr var = bvpOp->getOutputVariable();
        loadOp->setVariable( var );

        AMP::shared_ptr<AMP::Database> shell_db( new AMP::MemoryDatabase( "MatrixShellOperator" ) );
        shell_db->putString( "name", "MatShellOperator" );
        shell_db->putInteger( "print_info_level", 1 );
        AMP::shared_ptr<AMP::Operator::OperatorParameters> shellParams(
            new AMP::Operator::OperatorParameters( shell_db ) );
        AMP::shared_ptr<AMP::Operator::PetscMatrixShellOperator> shellOp(
            new AMP::Operator::PetscMatrixShellOperator( shellParams ) );
        int const numLocalNodes = meshAdapter->numLocalElements( AMP::Mesh::Vertex );
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

        AMP::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(
            new AMP::Solver::PetscKrylovSolverParameters( linearSolver_db ) );
        linearSolverParams->d_pOperator       = shellOp;
        linearSolverParams->d_comm            = globalComm;
        linearSolverParams->d_pPreconditioner = colPre;
        AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(
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

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    myTest( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
