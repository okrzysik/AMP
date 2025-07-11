#include "AMP/IO/PIO.h"
#include "AMP/utils/AMPManager.h"

#include "AMP/vectors/CommunicationList.h"
#include "AMP/matrices/petsc/NativePetscMatrix.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/data/VectorData.h"
#include "AMP/vectors/data/VectorDataNull.h"
#include "AMP/vectors/operations/default/VectorOperationsDefault.h"

#include "AMP/discretization/boxMeshDOFManager.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshParameters.h"
#include "AMP/mesh/MeshElement.h"
#include "AMP/mesh/structured/BoxMesh.h"

#include "AMP/matrices/CSRMatrix.h"
#include "AMP/matrices/CSRPolicy.h"
#include "AMP/matrices/MatrixBuilder.h"
#include "AMP/matrices/RawCSRMatrixParameters.h"

#include "AMP/operators/Operator.h"
#include "AMP/operators/LinearOperator.h"
#include "AMP/operators/OperatorFactory.h"

#include "AMP/solvers/SolverFactory.h"
#include "AMP/solvers/SolverStrategy.h"
#include "AMP/solvers/testHelpers/SolverTestParameters.h"
#include "AMP/solvers/SolverStrategyParameters.h"
#include "AMP/solvers/SolverFactory.h"

#include "AMP/time_integrators/TimeIntegratorParameters.h"
#include "AMP/time_integrators/TimeIntegratorFactory.h"
#include "AMP/time_integrators/TimeIntegrator.h"
#include "AMP/time_integrators/TimeOperator.h"

#include <iostream>


/* Linear system we'll have to solve is [I - gamma*dt]*u = f. For gamma*dt < 1 the matrix is SPD. */
void updateDatabaseIfImplicit( std::shared_ptr<AMP::Database> db, std::string implicitSolverName )
{
    AMP_ASSERT( db );
    auto imp_ti      = { "Backward Euler", "BDF1", "BDF2", "BDF3", "BDF4", "BDF5", "BDF6" };
    auto name        = db->getScalar<std::string>( "name" );
    auto is_implicit = ( std::find( imp_ti.begin(), imp_ti.end(), name ) != imp_ti.end() );
    db->putScalar<bool>( "is_implicit", is_implicit );
    if ( is_implicit ) {
        db->putScalar<std::string>( "implicit_integrator", name );
        db->putScalar<std::string>( "solver_name", "Solver" );
        db->putScalar<std::string>( "timestep_selection_strategy", "constant" );
        db->putScalar<bool>( "use_predictor", false );
        int print_info_level = 0;
        if (implicitSolverName == "CGSolver") 
            print_info_level = 2;
        auto solver_db = AMP::Database::create( "name",
                                                implicitSolverName,
                                                "print_info_level",
                                                print_info_level,
                                                "max_iterations",
                                                100,
                                                "absolute_tolerance",
                                                1.0e-12,
                                                "relative_tolerance",
                                                1.0e-12 );
        db->putDatabase( "Solver", std::move( solver_db ) );
    }
}

/* 
    Looking in ExplicitEuler, we solve ODEs of the form
        u'(t) = f(t, u)
    So we must provide an apply for the function f. Note that 
        f == d_operator + d_pSourceTerm
    with d_operator and d_pSourceTerm defined in TimeIntegratorParameters
*/

/* We want to solve an ODE of the form u'(t) = f(t, u) == -u(t) + g_source(t) */
double g_source( double a, double t )      { return a * ( 2*M_PI*cos(2*M_PI*t) + sin(2*M_PI*t) );  }
double exactSolution( double a, double t ) { return std::exp( -t ) + a*sin(2*M_PI*t); }

/* This operator represents u(t) -> -u(t) */
class MyOperator : public AMP::Operator::LinearOperator {
public:
    std::shared_ptr<AMP::Discretization::DOFManager> d_DOFManager;

    // Constructor
    MyOperator( std::shared_ptr<const AMP::Operator::OperatorParameters> params_) : 
            AMP::Operator::LinearOperator( params_ ) { }

    // Create and set DOFManager for single variable
    void setDOFManager( AMP::AMP_MPI comm ) {
        size_t num_dofs = 1;
        d_DOFManager = std::make_shared<AMP::Discretization::DOFManager>( num_dofs, comm );
    }  

    // out <- -1*in
    void apply( std::shared_ptr<const AMP::LinearAlgebra::Vector> in, std::shared_ptr<AMP::LinearAlgebra::Vector> out ) override {
        out->scale (-1.0, *in);
    }

    std::string type() const override { return "MyOperator"; }

    std::shared_ptr<AMP::LinearAlgebra::Vector> getRightVector() const override {
        AMP_INSIST( d_DOFManager, "Must call setDOFManager first!" );
        auto tempVar = std::make_shared<AMP::LinearAlgebra::Variable> (" ");
        return AMP::LinearAlgebra::createVector( this->d_DOFManager, tempVar );
    }
};



void driver( AMP::AMP_MPI comm, double dt, std::string integratorName, std::string implicitSolverName ) {

    double a = 0.5;

    // Print what we're doing
    AMP::pout << "-------------------------------------------" << std::endl;
    AMP::pout << "Time integrator=" << integratorName << ", with implicit solver=" << implicitSolverName << std::endl;
    AMP::pout << "-------------------------------------------" << std::endl;

    // Create a MyOperator
    const auto opDB = std::make_shared<AMP::Database>( "OperatorDB" );
    opDB->putScalar<std::string>( "name", "MyOperator" );  
    auto opParams = std::make_shared<AMP::Operator::OperatorParameters>( opDB );
    auto myOp = std::make_shared<MyOperator>( opParams );
    myOp->setDOFManager( comm );

    // Parameters for time integraor
    double finalTime       = 1;
    // dt              = 0.25 / 2;
    int maxIntegratorSteps = 300;

    // Create initial condition vector
    auto ic = myOp->getRightVector();
    ic->setToScalar( 1.0 );
    ic->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );

    // Create vectors to hold current and new solutions
    auto sol_old = ic->clone( );
    sol_old->copy( *ic ); 
    auto sol_new = ic->clone( );
    
    // Create Database to construct TimeIntegratorParameters from
    std::shared_ptr<AMP::Database> db = AMP::Database::create( "name",
                                                               integratorName,
                                                               "initial_time",
                                                               0.0,
                                                               "final_time",
                                                               finalTime,
                                                               "max_integrator_steps",
                                                               maxIntegratorSteps,
                                                               "print_info_level",
                                                               0 );
    auto tiParams           = std::make_shared<AMP::TimeIntegrator::TimeIntegratorParameters>( db );
    tiParams->d_ic_vector   = ic;
    tiParams->d_operator    = myOp;
    // Create a source vector
    tiParams->d_pSourceTerm = myOp->getRightVector();
    
    updateDatabaseIfImplicit( db, implicitSolverName );

    // Create timeIntegrator from factory
    auto timeIntegrator = AMP::TimeIntegrator::TimeIntegratorFactory::create( tiParams );
    
    // Integrate!
    double T = 0.0;
    timeIntegrator->setInitialDt( dt );
    while ( T < finalTime ) {

        // Set the solution-independent source term; note that this approach can only work for multistep methods, since then the new solution 
        double g;
        if ( db->getScalar<bool>( "is_implicit" ) ) {
            g = g_source( a, T + dt ); // Set source to new time...
        } else {
            g = g_source( a, T ); // Set source to current time...
        }
        tiParams->d_pSourceTerm->setValueByGlobalID<double>( 0, g );

        // Advance the solution
        timeIntegrator->advanceSolution( dt, T == 0, sol_old, sol_new );
        if ( timeIntegrator->checkNewSolution() ) {
            timeIntegrator->updateSolution();
            sol_old->copyVector( sol_new );
        } else {
            AMP_ERROR( "Solution didn't converge" );
        }

        // Update time
        T += dt;

        // Print exact vs approximate solution 
        auto u_exact = exactSolution( a, T );
        auto u_num   = sol_new->getValueByGlobalID( 0 );
        auto e       = std::abs( u_exact - u_num );
        AMP::pout << "t=" << T << ": exact-sol=" << u_exact << ", num-sol=" << u_num << ", error=" << e << std::endl;

        // Drop out if we've exceeded max steps
        if ( !timeIntegrator->stepsRemaining() ) {
            AMP_WARNING( "max_integrator_steps has been reached, dropping out of loop now..." );
            break;
        }
    }
}
// end of driver()


int main( int argc, char **argv )
{

    AMP::AMPManager::startup( argc, argv );

    double dt = 0.1;
    if (argc > 1) {
        dt = std::stod(argv[1]);
    }

    // Create a global communicator
    AMP::AMP_MPI comm( AMP_COMM_WORLD );

    
    // Driver
    driver( comm, dt, "ExplicitEuler", "" );
    //driver( comm, dt, "BDF3", "PetscSNESSolver" );
    //driver( comm, "BDF1", "CGSolver" );
    //driver( comm, "BDF1", "GMRESSolver" );

    AMP::AMPManager::shutdown();

    return 0;
}