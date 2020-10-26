#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/vectors/MultiVariable.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/testHelpers/VectorTests.h"
#include "AMP/vectors/testHelpers/generateVectorFactories.h"

#include "test_ArrayVector.h"

#ifdef USE_TRILINOS_THYRA
#include "AMP/vectors/testHelpers/trilinos/thyra/ThyraVectorFactory.h"
#endif

#include <iostream>


using namespace AMP::unit_test;
using namespace AMP::LinearAlgebra;
using AMP::LinearAlgebra::generateVectorFactory;
using AMP::LinearAlgebra::VectorTests;


std::string SimpleFactory1 = "SimpleVectorFactory<15,false,double>";
std::string SimpleFactory2 = "SimpleVectorFactory<45,true,double>";
std::string SNPVFactory    = "NativePetscVectorFactory";
std::string SMEVFactory    = "ManagedEpetraVectorFactory<SimpleVectorFactory<45,true>>";
#ifdef USE_EXT_PETSC
std::string MVFactory1 = "MultiVectorFactory<" + SMEVFactory + ", 1, " + SNPVFactory + ", 1>";
#else
std::string MVFactory1 =
    "MultiVectorFactory<SimpleVectorFactory<15,false>,1," + SMEVFactory + ",1>";
#endif


int main( int argc, char **argv )
{

    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    // Test ArrayVector dimensions
    std::vector<size_t> dims{ 3, 3, 3, 3 };
    testArrayVectorDimensions<double>( dims, ut );

    // Print the total nummber of vector factories to test
    AMP::pout << "Testing " << getAllFactories().size() << " vector factories\n";
    for ( auto name : getManagedVectorFactories() ) {
        auto factory = generateVectorFactory( name );
        auto name2   = factory->name();
        if ( name != name2 ) {
            ut.failure( "Factory name " + name2 + " does not match input name " + name );
        }
    }

    // Run null vector tests
    AMP::pout << std::endl << "Testing NullVector" << std::endl;
    {
        auto factory = generateVectorFactory( SimpleFactory1 );
        VectorTests tests( factory );
        tests.testNullVector( &ut );
    }

    // Run the basic vector tests
    AMP::pout << std::endl << "Running basic vector tests:" << std::endl;
    for ( auto name : getAllFactories() ) {
        auto factory = generateVectorFactory( name );
        VectorTests tests( factory );
        tests.testBasicVector( &ut );
    }

    // Run the managed vector tests
    AMP::pout << std::endl << "Running managed vector tests:" << std::endl;
    for ( auto name : getManagedVectorFactories() ) {
        auto factory = generateVectorFactory( name );
        VectorTests tests( factory );
        tests.testManagedVector( &ut );
    }

    // Run the petsc vector tests
    AMP::pout << std::endl << "Running petsc vector tests:" << std::endl;
    for ( auto name : getManagedVectorFactories() ) {
        auto factory = generateVectorFactory( name );
        VectorTests tests( factory );
        tests.testPetsc( &ut );
    }

    // Run the sundials vector tests
    AMP::pout << std::endl << "Running sundials vector tests:" << std::endl;
    for ( auto name : getManagedVectorFactories() ) {
        auto factory = generateVectorFactory( name );
        VectorTests tests( factory );
        tests.testSundials( &ut );
    }

    // Run the vector selector tests
    AMP::pout << std::endl << "Running vector selector tests:" << std::endl;
    for ( auto name : getAllFactories() ) {
        auto factory = generateVectorFactory( name );
        VectorTests tests( factory );
        tests.testVectorSelector( &ut );
    }

// Run Belos tests of thyra vectors
#if defined( USE_TRILINOS_THYRA ) && defined( USE_TRILINOS_BELOS )
    AMP::pout << std::endl << "Testing Belos interface to Thyra vectors" << std::endl;
    testBelosThyraVector( ut, NativeThyraFactory() );
    testBelosThyraVector( ut, ManagedThyraFactory( generateVectorFactory( SimpleFactory2 ) ) );
    testBelosThyraVector( ut,
                          ManagedNativeThyraFactory( generateVectorFactory( SimpleFactory2 ) ) );
    testBelosThyraVector( ut, ManagedNativeThyraFactory( generateVectorFactory( MVFactory1 ) ) );
#endif

    ut.report();
    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
