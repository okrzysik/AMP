#include "AMP/AMP_TPLs.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/MultiVariable.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/testHelpers/VectorTests.h"
#include "AMP/vectors/testHelpers/generateVectorFactories.h"

#include "test_ArrayVector.h"

#ifdef AMP_USE_TRILINOS_THYRA
    #include "AMP/vectors/testHelpers/trilinos/thyra/ThyraVectorFactory.h"
#endif

#include "ProfilerApp.h"

#include <iostream>

using namespace AMP::unit_test;
using namespace AMP::LinearAlgebra;
using AMP::LinearAlgebra::generateVectorFactory;
using AMP::LinearAlgebra::VectorTests;

std::string SimpleFactory1 = "SimpleVectorFactory<15,false,double>";
std::string SimpleFactory2 = "SimpleVectorFactory<45,true,double>";
std::string SNPVFactory    = "NativePetscVectorFactory";
std::string SNEVFactory    = "NativeEpetraFactory";
#ifdef AMP_USE_PETSC
std::string MVFactory1 = "MultiVectorFactory<" + SNEVFactory + ", 1, " + SNPVFactory + ", 1>";
#else
std::string MVFactory1 =
    "MultiVectorFactory<SimpleVectorFactory<15,false>,1," + SNEVFactory + ",1>";
#endif


template<class FUN, class... Args>
int testVectorSelectorPerformance( std::shared_ptr<AMP::LinearAlgebra::Vector> vec, Args... ts )
{
    int N_it = 10;
    auto t1  = AMP::Utilities::time();
    for ( int i = 0; i < N_it; i++ ) {
        FUN variable( ts... );
        AMP_ASSERT( variable.subset( vec ) );
    }
    auto t2 = AMP::Utilities::time();
    return 1e9 * ( t2 - t1 ) / N_it;
}
void testVectorSelectorPerformance()
{
    PROFILE( "testVectorSelectorPerformance" );
    AMP::AMP_MPI world_comm( AMP_COMM_WORLD );
    AMP::AMP_MPI self_comm( AMP_COMM_SELF );
    if ( world_comm.getRank() == 0 )
        printf( "       subset performance (ns):              "
                " Variable  Stride  VecComm    World    Self\n" );
    for ( auto name : getAllFactories() ) {
        auto factory  = generateVectorFactory( name );
        auto vec      = factory->getVector();
        auto vec_comm = vec->getComm();
        world_comm.barrier();
        auto t1 = testVectorSelectorPerformance<VS_ByVariableName>( vec, vec->getName() );
        auto t2 = testVectorSelectorPerformance<VS_Stride>( vec, 0, 1 );
        auto t3 = testVectorSelectorPerformance<VS_Comm>( vec, vec_comm );
        auto t4 = testVectorSelectorPerformance<VS_Comm>( vec, world_comm );
        auto t5 = testVectorSelectorPerformance<VS_Comm>( vec, self_comm );
        if ( world_comm.getRank() == 0 ) {
            auto name2 = name.substr( 0, 40 );
            printf( "  %40s  %8i %8i %8i %8i %8i\n", name2.data(), t1, t2, t3, t4, t5 );
        }
    }
}


int main( int argc, char **argv )
{

    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    PROFILE_ENABLE( 3 );

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
    for ( auto name : getAllFactories() ) {
        auto factory = generateVectorFactory( name );
        VectorTests tests( factory );
        tests.testPetsc( &ut );
    }

    // Run the epetra vector tests
    AMP::pout << std::endl << "Running epetra vector tests:" << std::endl;
    for ( auto name : getAllFactories() ) {
        auto factory = generateVectorFactory( name );
        VectorTests tests( factory );
        tests.testEpetra( &ut );
    }

    // Run the sundials vector tests
    AMP::pout << std::endl << "Running sundials vector tests:" << std::endl;
    for ( auto name : getAllFactories() ) {
        auto factory = generateVectorFactory( name );
        VectorTests tests( factory );
        tests.testSundials( &ut );
    }

// Run Belos tests of thyra vectors
#if defined( AMP_USE_TRILINOS_THYRA ) && defined( AMP_USE_TRILINOS_BELOS )
    AMP::pout << std::endl << "Testing Belos interface to Thyra vectors" << std::endl;
    testBelosThyraVector( ut, NativeThyraFactory() );
    testBelosThyraVector( ut, ManagedThyraFactory( generateVectorFactory( SimpleFactory2 ) ) );
    testBelosThyraVector( ut,
                          ManagedNativeThyraFactory( generateVectorFactory( SimpleFactory2 ) ) );
    testBelosThyraVector( ut, ManagedNativeThyraFactory( generateVectorFactory( MVFactory1 ) ) );
#endif

    // Run the vector selector tests
    AMP::pout << std::endl << "Running vector selector tests:" << std::endl;
    for ( auto name : getAllFactories() ) {
        auto factory = generateVectorFactory( name );
        VectorTests tests( factory );
        tests.testVectorSelector( &ut );
    }

    // Test vector selector performance
    testVectorSelectorPerformance();

    PROFILE_SAVE( "testVector" );
    ut.report();
    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
