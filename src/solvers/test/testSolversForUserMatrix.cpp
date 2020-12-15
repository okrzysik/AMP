#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MeshParameters.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/matrices/MatrixBuilder.h"
#include "AMP/operators/LinearOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/OperatorParameters.h"
#include "AMP/solvers/KrylovSolverParameters.h"
#include "AMP/solvers/SolverFactory.h"
#include "AMP/solvers/hypre/BoomerAMGSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include <iostream>
#include <memory>
#include <string>

std::shared_ptr<AMP::Solver::SolverStrategy>
buildSolver( const std::shared_ptr<AMP::Database> &input_db,
             const std::string &solver_name,
             const AMP::AMP_MPI &comm,
             std::shared_ptr<AMP::Operator::LinearOperator> &op )
{

    std::shared_ptr<AMP::Solver::SolverStrategy> solver;
    std::shared_ptr<AMP::Solver::SolverStrategyParameters> parameters;

    AMP_INSIST( input_db->keyExists( solver_name ), "Key " + solver_name + " is missing!" );

    const auto &db = input_db->getDatabase( solver_name );

    if ( db->keyExists( "name" ) ) {

        auto name = db->getString( "name" );

        if ( ( name == "GMRESSolver" ) || ( name == "CGSolver" ) || ( name == "BiCGSTABSolver" ) ) {

            // check if we need to construct a preconditioner
            auto use_preconditioner = db->getWithDefault<bool>( "use_preconditioner", false );
            std::shared_ptr<AMP::Solver::SolverStrategy> pcSolver;

            if ( use_preconditioner ) {

                auto pc_name = db->getWithDefault<std::string>( "pc_name", "Preconditioner" );

                auto pcSolver = buildSolver( input_db, pc_name, comm, op );

                AMP_INSIST( pcSolver.get() != nullptr, "null preconditioner" );
            }

            auto params               = std::make_shared<AMP::Solver::KrylovSolverParameters>( db );
            params->d_comm            = comm;
            params->d_pPreconditioner = pcSolver;
            parameters                = params;

        } else {
            parameters = std::make_shared<AMP::Solver::SolverStrategyParameters>( db );
        }

        AMP_INSIST( parameters != nullptr, "null parameter object" );
        parameters->d_pOperator = op;

    } else {
        AMP_ERROR( "Key name does not exist in solver database" );
    }

    solver = AMP::Solver::SolverFactory::create( parameters );

    return solver;
}

void userLinearOperatorTest( AMP::UnitTest *const ut, const std::string &inputFileName )
{
    // Test create
    const std::string input_file = inputFileName;
    const std::string log_file   = "output_" + inputFileName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    // read the input file into a database
    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // extract the Mesh database and create the mesh parameters
    const auto meshDB = input_db->getDatabase( "Mesh" );
    auto params       = std::make_shared<AMP::Mesh::MeshParameters>( meshDB );
    params->setComm( globalComm );

    // create the mesh
    const auto meshAdapter = AMP::Mesh::Mesh::buildMesh( params );

    // create a linear diffusion operator
    auto linearOperator = AMP::Operator::OperatorBuilder::createOperator(
        meshAdapter, "DiffusionBVPOperator", input_db );
    auto diffOp = std::dynamic_pointer_cast<AMP::Operator::LinearOperator>( linearOperator );

    // extract the internal matrix
    const auto &userMat = diffOp->getMatrix();

    AMP_INSIST( userMat->numGlobalColumns() == userMat->numGlobalRows(), "matrix is not square" );

    // extract the right vector
    const auto userVector = userMat->getRightVector();

    // concludes creation of a native linear operator
    // ************************************************************************************************

    // extract information about the local size and mpi comm
    const auto localSize = userVector->getLocalSize();
    const auto ampComm   = userVector->getComm();

    // construct a dof manager
    const auto dofManager = std::make_shared<AMP::Discretization::DOFManager>( localSize, ampComm );
    const auto copyVariable = std::make_shared<AMP::LinearAlgebra::Variable>( "copyVariable" );

    // create a vector based on the dofs and variable
    auto ampVector = AMP::LinearAlgebra::createVector( dofManager, copyVariable );
    AMP_INSIST( ampVector != nullptr, "ampVector is null" );

    // copy values from one vector to another
    std::copy( userVector->begin(), userVector->end(), ampVector->begin() );
    ampVector->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
    // concludes demonstrating how to initialize an AMP vector from a user vector
    // ************************************************************************************************

    // create a lambda that returns non zero column ids given a global row id
    auto getColumnIDS = [userMat]( size_t row ) { return userMat->getColumnIDs( row ); };

    // create a matrix based on the dimensions of the copied vector
    auto ampMat = AMP::LinearAlgebra::createMatrix( ampVector, ampVector, "auto", getColumnIDS );

    // construct a LinearOperator and set its matrix
    const auto linearOpDB = std::make_shared<AMP::Database>( "linearOperatorDB" );
    linearOpDB->putScalar<int>( "print_info_level", 0 );
    auto linearOpParameters = std::make_shared<AMP::Operator::OperatorParameters>( linearOpDB );
    auto linearOp           = std::make_shared<AMP::Operator::LinearOperator>( linearOpParameters );
    linearOp->setMatrix( ampMat );
    linearOp->setVariables( copyVariable, copyVariable );

    // copy the user matrix into the amp matrix
    std::vector<double> coefficients;
    std::vector<size_t> cols;
    const size_t numRows = 1;
    for ( auto row = userMat->beginRow(); row < userMat->endRow(); ++row ) {
        userMat->getRowByGlobalID( row, cols, coefficients );
        ampMat->setValuesByGlobalID( numRows, cols.size(), &row, cols.data(), coefficients.data() );
    }
    // concludes demonstrating how to initialize an AMP linear operator from a user matrix
    // ************************************************************************************************

    auto f  = ampVector->cloneVector();
    auto u  = ampVector->cloneVector();
    auto ri = ampVector->cloneVector();
    auto rf = ampVector->cloneVector();

    f->zero();
    u->setRandomValues();

    // ************************************************************************************************
    // make sure the database on theinput file exists for the linear solver
    AMP_INSIST( input_db->keyExists( "LinearSolver" ), "Key ''LinearSolver'' is missing!" );

    auto linearSolver = buildSolver( input_db, "LinearSolver", globalComm, linearOp );

    // Use a random initial guess
    linearSolver->setZeroInitialGuess( false );

    // Compute the initial residual
    linearOp->residual( f, u, ri );

    // Solve the problem.
    linearSolver->solve( f, u );

    // Compute the final residual
    linearOp->residual( f, u, rf );

    // Check the L2 norm of the residuals.
    const double finalResidualNorm   = static_cast<double>( rf->L2Norm() );
    const double initialResidualNorm = static_cast<double>( ri->L2Norm() );
    AMP::pout << "Initial Residual Norm: " << initialResidualNorm << std::endl;
    AMP::pout << "Final Residual Norm  : " << finalResidualNorm << std::endl;

    if ( finalResidualNorm / initialResidualNorm > 1.0e-7 ) {
        ut->failure( "solver could NOT solve a linear thermal problem" );
    }
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    AMP::Solver::registerSolverFactories();

    std::vector<std::string> files;

    if ( argc > 1 ) {

        files.push_back( argv[1] );

    } else {

        files.push_back( "input_testSolversForUserMatrix-ML" );
#ifdef USE_EXT_HYPRE
        files.push_back( "input_testSolversForUserMatrix-BoomerAMG" );
#endif
#ifdef USE_TRILINOS_MUELU
        files.push_back( "input_testSolversForUserMatrix-MueLu" );
#endif
    }

    for ( const auto &file : files )
        userLinearOperatorTest( &ut, file );

    ut.report();

    const int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
