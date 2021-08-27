#include "AMP/time_integrators/TimeIntegratorFactory.h"
#include "AMP/time_integrators/ExplicitEuler.h"
#include "AMP/time_integrators/RK23TimeIntegrator.h"
#include "AMP/time_integrators/RK2TimeIntegrator.h"
#include "AMP/time_integrators/RK4TimeIntegrator.h"
#include "AMP/utils/Utilities.h"


namespace AMP {
namespace TimeIntegrator {

TimeIntegratorFactory::TimeIntegratorFactory() = default;

TimeIntegratorFactory::~TimeIntegratorFactory() = default;

std::shared_ptr<TimeIntegrator> TimeIntegratorFactory::createTimeIntegrator(
    std::shared_ptr<TimeIntegratorParameters> timeIntegratorParameters )
{
    AMP_ASSERT( timeIntegratorParameters );

    std::shared_ptr<TimeIntegrator> timeIntegrator;

    std::string timeIntegratorName = "";

    std::shared_ptr<AMP::Database> db( timeIntegratorParameters->d_db );

    if ( db->keyExists( "timeIntegrator_name" ) ) {
        timeIntegratorName = db->getString( "timeIntegrator_name" );
    } else {
        AMP_ERROR( "TimeIntegratorFactory"
                   << " -- Required key `timeIntegrator_name'"
                   << " missing in input." );
    }

    if ( timeIntegratorName == "ExplicitEuler" ) {
        timeIntegrator.reset( new ExplicitEuler( timeIntegratorParameters ) );
    } else if ( timeIntegratorName == "RK2" ) {
        timeIntegrator.reset( new RK2TimeIntegrator( timeIntegratorParameters ) );
    } else if ( timeIntegratorName == "RK4" ) {
        timeIntegrator.reset( new RK4TimeIntegrator( timeIntegratorParameters ) );
    } else if ( timeIntegratorName == "RK23" ) {
        timeIntegrator.reset( new RK23TimeIntegrator( timeIntegratorParameters ) );
    } else {
        AMP_ERROR(
            std::string( "TimeIntegratorFactory does not currently create timeIntegrator " ) +
            timeIntegratorName );
    }

    return timeIntegrator;
}
} // namespace TimeIntegrator
} // namespace AMP
