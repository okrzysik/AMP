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
std::string SNPVFactory    = "SimplePetscNativeFactory";
std::string SMEVFactory    = "ManagedEpetraVectorFactory<SimpleVectorFactory<45,true>>";
#ifdef USE_EXT_PETSC
std::string MVFactory1 = "MultiVectorFactory<" + SMEVFactory + ", 1, " + SNPVFactory + ", 1>";
#else
std::string MVFactory1 =
    "MultiVectorFactory<SimpleVectorFactory<15,false>,1," + SMEVFactory + ",1>";
#endif


// Set the tests
inline void testNullVector( AMP::UnitTest &ut, const std::string &factoryName )
{
    auto factory = generateVectorFactory( factoryName );
    VectorTests tests( factory );
    tests.testNullVector( &ut );
    AMP::AMP_MPI( AMP_COMM_WORLD ).barrier();
}
inline void testBasicVector( AMP::UnitTest &ut, const std::string &factoryName )
{
    auto factory = generateVectorFactory( factoryName );
    VectorTests tests( factory );
    tests.testBasicVector( &ut );
    AMP::AMP_MPI( AMP_COMM_WORLD ).barrier();
}
inline void testManagedVector( AMP::UnitTest &ut, const std::string &factoryName )
{
    auto factory = generateVectorFactory( factoryName );
    VectorTests tests( factory );
    tests.testManagedVector( &ut );
    AMP::AMP_MPI( AMP_COMM_WORLD ).barrier();
}
inline void testVectorSelector( AMP::UnitTest &ut, const std::string &factoryName )
{
    auto factory = generateVectorFactory( factoryName );
    VectorTests tests( factory );
    tests.testVectorSelector( &ut );
    AMP::AMP_MPI( AMP_COMM_WORLD ).barrier();
}
inline void checkFactoryName( AMP::UnitTest &ut, const std::string &factoryName )
{
    auto factory = generateVectorFactory( factoryName );
    auto name    = factory->name();
    if ( name != factoryName ) {
        ut.failure( "Factory name " + name + " does not match input name " + factoryName );
    }
}


int main( int argc, char **argv )
{

    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    // Print the total nummber of vector factories to test
    std::cout << "Testing " << getAllFactories().size() << " vector factories\n";
    for ( auto factory : getManagedVectorFactories() )
        checkFactoryName( ut, factory );
    AMP::pout << std::endl;

    // Test ArrayVector dimensions
    std::vector<size_t> dims{ 3, 3, 3, 3 };
    testArrayVectorDimensions<double>( dims, ut );

    // Run null vector tests
    AMP::pout << "Testing NullVector" << std::endl;
    testNullVector( ut, SimpleFactory1 );

    // Run the basic vector tests
    AMP::pout << "Running basic vector tests:" << std::endl;
    for ( auto factory : getAllFactories() ) {
        // AMP::pout << "    " << factory << std::endl;
        testBasicVector( ut, factory );
    }
    AMP::pout << std::endl;


    // Run the managed vector tests
    AMP::pout << "Running managed vector tests:" << std::endl;
    for ( auto factory : getManagedVectorFactories() ) {
        // AMP::pout << "    " << factory << std::endl;
        testManagedVector( ut, factory );
    }
    AMP::pout << std::endl;

    // Run the vector selector tests
    AMP::pout << "Running vector selector tests:" << std::endl;
    for ( auto factory : getAllFactories() ) {
        // AMP::pout << "    " << factory << std::endl;
        testVectorSelector( ut, factory );
    }
    AMP::pout << std::endl;

// Run Belos tests of thyra vectors
#if defined( USE_TRILINOS_THYRA ) && defined( USE_TRILINOS_BELOS )
    AMP::pout << "Testing Belos interface to Thyra vectors" << std::endl;
    testBelosThyraVector( ut, NativeThyraFactory() );
    testBelosThyraVector( ut, ManagedThyraFactory( generateVectorFactory( SimpleFactory2 ) ) );
    testBelosThyraVector( ut,
                          ManagedNativeThyraFactory( generateVectorFactory( SimpleFactory2 ) ) );
    testBelosThyraVector( ut, ManagedNativeThyraFactory( generateVectorFactory( MVFactory1 ) ) );
    AMP::pout << std::endl;
#endif

    ut.report();
    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
