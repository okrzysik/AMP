#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "vectors/MultiVariable.h"
#include "vectors/MultiVector.h"
#include <iostream>

#include "test_VectorHelpers.h"


int main( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    AMP::pout << "Testing NullVector" << std::endl;
    testNullVector( ut, SimpleFactory1 );

    AMP::pout << "Testing SimpleVector" << std::endl;
    testBasicVector( ut, SimpleFactory1 );
    testBasicVector( ut, SimpleFactory2 );
    //testBasicVector( ut, "SimpleVectorFactory<15,false,float>" );
    AMP::pout << std::endl;

    AMP::pout << "Testing ArrayVector" << std::endl;
    testBasicVector( ut, ArrayFactory1 );
    testBasicVector( ut, ArrayFactory2 );
    std::vector<size_t> dims{ 3, 3, 3, 3 };
    testArrayVectorDimensions<double>( dims, ut );
    AMP::pout << std::endl;

#if defined( USE_EXT_PETSC )
    AMP::pout << "Testing NativePetscVector" << std::endl;
    testBasicVector( ut, NPVFactory );
    AMP::pout << std::endl;
#endif

#if defined( USE_EXT_PETSC ) && defined( USE_EXT_TRILINOS )
    AMP::pout << "Testing ManagedPetscVector<NativePetscVector>" << std::endl;
    testManagedVector( ut, SNPVFactory );
    AMP::pout << std::endl;
#else
    ut.expected_failure( "Compiled without petsc" );
#endif

#ifdef USE_EXT_TRILINOS

#ifdef USE_TRILINOS_THYRA
    std::string ManagedThyraFactory1 = "ManagedThyraFactory<"+SimpleFactory1+">";
    std::string ManagedThyraFactory2 = "ManagedThyraFactory<"+SimpleFactory2+">";
    std::string ManagedNativeThyraFactory1 = "ManagedNativeThyraFactory<"+SimpleFactory1+">";
    std::string ManagedNativeThyraFactory2 = "ManagedNativeThyraFactory<"+SimpleFactory2+">";

    AMP::pout << "Testing NativeThyraVector" << std::endl;
    testBasicVector( ut, "NativeThyraFactory" );
    AMP::pout << std::endl;

    AMP::pout << "Testing ManagedThyraVector" << std::endl;
    testBasicVector( ut, ManagedThyraFactory2 );
    AMP::pout << std::endl;

    AMP::pout << "Testing NativeThyraVector of a ManagedThyraVector" << std::endl;
    testBasicVector( ut, ManagedNativeThyraFactory2 );
    AMP::pout << std::endl;
#endif

    AMP::pout << "Testing Iterator" << std::endl;
    VectorIteratorTests( ut, MVFactory1 );
    AMP::pout << std::endl;

    AMP::pout << "Testing ManagedEpetraVector" << std::endl;
    testManagedVector( ut, SMEVFactory );
    AMP::pout << std::endl;

    AMP::pout << "Testing simple multivector" << std::endl;
    testManagedVector( ut, MVFactory1 );
    AMP::pout << std::endl;

    AMP::pout << "Testing bigger multivector" << std::endl;
    testManagedVector( ut, MVFactory2 );
    AMP::pout << std::endl;

    AMP::pout << "Testing multivector of multivector" << std::endl;
    testManagedVector( ut, MVFactory3 );
    AMP::pout << std::endl;

#ifdef USE_TRILINOS_THYRA
    AMP::pout << "Testing NativeThyraVector of a ManagedThyraVector of a MultVector" << std::endl;
    testBasicVector( ut, ManagedNativeThyraFactory1 );
    AMP::pout << std::endl;
#endif

// Run Belos tests of thyra vectors
#ifdef USE_TRILINOS_THYRA
    AMP::pout << "Testing Belos interface to Thyra vectors" << std::endl;
    testBelosThyraVector( ut, NativeThyraFactory() );
    testBelosThyraVector( ut, ManagedThyraFactory( generateVectorFactory(SimpleFactory2) ) );
    testBelosThyraVector( ut, ManagedNativeThyraFactory( generateVectorFactory(SimpleFactory2) ) );
    testBelosThyraVector( ut, ManagedNativeThyraFactory( generateVectorFactory(MVFactory1) ) );
    AMP::pout << std::endl;
#endif

#else
    ut.expected_failure( "Compiled without trilinos" );
#endif

    ut.report();

    int num_failed = ut.NumFailGlobal();
    if ( num_failed == 0 )
        AMP::pout << "No errors detected" << std::endl;
    AMP::AMPManager::shutdown();
    return num_failed;
}
