#include "AMP/IO/PIO.h"
#include "AMP/operators/NullOperator.h"
#include "AMP/operators/Operator.h"
#include "AMP/time_integrators/TimeIntegrator.h"
#include "AMP/time_integrators/TimeIntegratorFactory.h"
#include "AMP/time_integrators/TimeIntegratorParameters.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"


class FunctionOperator : public AMP::Operator::Operator
{
public:
    FunctionOperator( std::function<double( double )> f ) : d_f( f ) {}
    std::string type() const override { return "FunctionOperator"; }
    void apply( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                std::shared_ptr<AMP::LinearAlgebra::Vector> r ) override
    {
        auto it_f = f->begin();
        auto it_r = r->begin();
        for ( size_t i = 0; i < it_f.size(); ++i, ++it_f, ++it_r )
            *it_r = d_f( *it_f );
    }

private:
    std::function<double( double )> d_f;
};


void testIntegrator( const std::string &name,
                     const std::string &test,
                     std::shared_ptr<AMP::TimeIntegrator::TimeIntegratorParameters> params,
                     double ans,
                     AMP::UnitTest &ut )
{
    // Create the time integrator
    auto var            = std::make_shared<AMP::LinearAlgebra::Variable>( "x" );
    auto solution       = AMP::LinearAlgebra::createSimpleVector<double>( 1, var, AMP_COMM_WORLD );
    auto timeIntegrator = AMP::TimeIntegrator::TimeIntegratorFactory::create( params );
    solution->zero();

    // Advance the solution
    double T         = 0;
    double dt        = 0.1;
    double finalTime = timeIntegrator->getFinalTime();
    timeIntegrator->setInitialDt( 0.1 );
    while ( T < finalTime ) {
        timeIntegrator->advanceSolution( dt, T == 0, solution, solution );
        timeIntegrator->setCurrentTime( T );
        T += dt;
    }

    // Check the answer
    double ans2 = static_cast<double>( solution->max() );
    if ( fabs( ans2 - ans ) < 1e-6 )
        ut.passes( name + " - " + test );
    else if ( test == "fixed source" )
        ut.expected_failure(
            AMP::Utilities::stringf( "%s - %s (%f)", name.data(), test.data(), ans2 ) );
    else
        ut.failure( AMP::Utilities::stringf( "%s - %s (%f)", name.data(), test.data(), ans2 ) );
}

void runBasicIntegratorTests( const std::string &name, AMP::UnitTest &ut )
{
    double finalTime = 2.0;

    // Create the vectors
    auto var = std::make_shared<AMP::LinearAlgebra::Variable>( "x" );
    auto ic  = AMP::LinearAlgebra::createSimpleVector<double>( 1, var, AMP_COMM_WORLD );
    ic->zero();

    // Test creating Create the time integrator
    std::shared_ptr<AMP::Database> db = AMP::Database::create(
        "name", name, "initial_time", 0.0, "final_time", finalTime, "max_integrator_steps", 10000 );
    auto params         = std::make_shared<AMP::TimeIntegrator::TimeIntegratorParameters>( db );
    params->d_ic_vector = ic;
    params->d_operator  = std::make_shared<AMP::Operator::NullOperator>();
    try {
        auto timeIntegrator = AMP::TimeIntegrator::TimeIntegratorFactory::create( params );
        ut.passes( name + " - created" );
    } catch ( ... ) {
        ut.failure( name + " - created" );
        return;
    }

    // Test with a fixed source and null operator
    auto source = AMP::LinearAlgebra::createSimpleVector<double>( 1, var, AMP_COMM_WORLD );
    source->setToScalar( 3.0 );
    params->d_pSourceTerm = source;
    params->d_operator    = std::make_shared<AMP::Operator::NullOperator>();
    testIntegrator( name, "fixed source", params, 3 * finalTime, ut );

    // Test with no source and constant operator
    params->d_pSourceTerm = nullptr;
    params->d_operator    = std::make_shared<FunctionOperator>( []( double ) { return 3.0; } );
    testIntegrator( name, "constant", params, 3 * finalTime, ut );
}


int testSimpleTimeIntegration( int argc, char *argv[] )
{

    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    // List of integrators
    auto integrators = { "ExplicitEuler", "RK12", "RK23", "RK34", "RK45", "RK2", "RK4" };
    //   "Backward Euler", "BDF1", "BDF2","BDF3", "BDF4","BDF5","BDF6"

    // Run the tests
    AMP::TimeIntegrator::registerTimeIntegratorFactories();
    for ( auto tmp : integrators )
        runBasicIntegratorTests( tmp, ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
