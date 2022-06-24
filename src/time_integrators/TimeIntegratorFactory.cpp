/*
Copyright 2005, The Regents of the University
of California. This software was produced under
a U.S. Government contract (W-7405-ENG-36)
by Los Alamos National Laboratory, which is
operated by the University of California for the
U.S. Department of Energy. The U.S.
Government is licensed to use, reproduce, and
distribute this software. Permission is granted
to the public to copy and use this software
without charge, provided that this Notice and
any statement of authorship are reproduced on
all copies. Neither the Government nor the
University makes any warranty, express or
implied, or assumes any liability or
responsibility for the use of this software.
*/
#include "AMP/time_integrators/TimeIntegratorFactory.h"
#include "AMP/time_integrators/BDFIntegrator.h"
#include "AMP/time_integrators/ExplicitEuler.h"
#include "AMP/time_integrators/RK12TimeIntegrator.h"
#include "AMP/time_integrators/RK23TimeIntegrator.h"
#include "AMP/time_integrators/RK2TimeIntegrator.h"
#include "AMP/time_integrators/RK34TimeIntegrator.h"
#include "AMP/time_integrators/RK45TimeIntegrator.h"
#include "AMP/time_integrators/RK4TimeIntegrator.h"
#include "AMP/time_integrators/TimeIntegrator.h"
#include "AMP/time_integrators/TimeIntegratorParameters.h"


namespace AMP::TimeIntegrator {


// Create the operator
std::unique_ptr<TimeIntegrator>
TimeIntegratorFactory::create( std::shared_ptr<TimeIntegratorParameters> parameters )
{
    AMP_ASSERT( parameters != nullptr );
    auto inputDatabase = parameters->d_db;
    AMP_ASSERT( inputDatabase );
    auto objectName = inputDatabase->getString( "name" );
    return FactoryStrategy<TimeIntegrator, std::shared_ptr<TimeIntegratorParameters>>::create(
        objectName, parameters );
}


// register all known time integrator factories
void registerTimeIntegratorFactories()
{
    auto &timeIntegratorFactory = TimeIntegratorFactory::getFactory();

    timeIntegratorFactory.registerFactory(
        "ExplicitEuler", AMP::TimeIntegrator::ExplicitEuler::createTimeIntegrator );
    timeIntegratorFactory.registerFactory( "ImplicitIntegrator",
                                           BDFIntegrator::createTimeIntegrator );
    timeIntegratorFactory.registerFactory( "Backward Euler", BDFIntegrator::createTimeIntegrator );
    timeIntegratorFactory.registerFactory( "BDF1", BDFIntegrator::createTimeIntegrator );
    timeIntegratorFactory.registerFactory( "BDF2", BDFIntegrator::createTimeIntegrator );
    timeIntegratorFactory.registerFactory( "BDF3", BDFIntegrator::createTimeIntegrator );
    timeIntegratorFactory.registerFactory( "BDF4", BDFIntegrator::createTimeIntegrator );
    timeIntegratorFactory.registerFactory( "BDF5", BDFIntegrator::createTimeIntegrator );
    timeIntegratorFactory.registerFactory( "BDF6", BDFIntegrator::createTimeIntegrator );
    timeIntegratorFactory.registerFactory(
        "RK2", AMP::TimeIntegrator::RK2TimeIntegrator::createTimeIntegrator );
    timeIntegratorFactory.registerFactory(
        "RK4", AMP::TimeIntegrator::RK4TimeIntegrator::createTimeIntegrator );
    timeIntegratorFactory.registerFactory(
        "RK12", AMP::TimeIntegrator::RK12TimeIntegrator::createTimeIntegrator );
    timeIntegratorFactory.registerFactory(
        "RK23", AMP::TimeIntegrator::RK23TimeIntegrator::createTimeIntegrator );
    timeIntegratorFactory.registerFactory(
        "RK34", AMP::TimeIntegrator::RK34TimeIntegrator::createTimeIntegrator );
    timeIntegratorFactory.registerFactory(
        "RK45", AMP::TimeIntegrator::RK45TimeIntegrator::createTimeIntegrator );
}


} // namespace AMP::TimeIntegrator
