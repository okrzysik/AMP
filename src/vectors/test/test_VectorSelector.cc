#include "test_VectorSelector.h"

#include "vectors/MultiVector.h"
#include "vectors/VectorSelector.h"

#include "utils/AMPManager.h"
#include "utils/UnitTest.h"

#ifdef USE_EXT_TRILINOS
#include "vectors/trilinos/ManagedEpetraVector.h"
#endif

#include "vectors/testHelpers/VectorTests.h"
#include "vectors/testHelpers/testVectorFactory.h"


// Typedef some factories
// clang-format off
typedef AMP::LinearAlgebra::SimpleVectorFactory<15,false,double>  SimpleFactory1;
typedef AMP::LinearAlgebra::SimpleVectorFactory<45, true,double>  SimpleFactory2;
typedef AMP::LinearAlgebra::ArrayVectorFactory<4,10,false,double> ArrayFactory1;
typedef AMP::LinearAlgebra::ArrayVectorFactory<4,10, true,double> ArrayFactory2;
#if defined( USE_EXT_PETSC ) && defined( USE_EXT_TRILINOS )
    typedef AMP::LinearAlgebra::SimplePetscNativeFactory<AMP::LinearAlgebra::NativePetscVector> SNPVFactory;
#endif
#ifdef USE_EXT_TRILINOS
    typedef AMP::LinearAlgebra::SimpleManagedVectorFactory<AMP::LinearAlgebra::ManagedEpetraVector> SMEVFactory;
#ifdef USE_EXT_PETSC
    typedef AMP::LinearAlgebra::MultiVectorFactory<SMEVFactory,1,SNPVFactory,1> MVFactory1;
    typedef AMP::LinearAlgebra::MultiVectorFactory<SMEVFactory,3,SNPVFactory,2> MVFactory2;
    typedef AMP::LinearAlgebra::MultiVectorFactory<MVFactory1, 2,MVFactory2, 2> MVFactory3;
#else
    typedef AMP::LinearAlgebra::MultiVectorFactory<SimpleFactory1,1,SMEVFactory,1> MVFactory1;
    typedef AMP::LinearAlgebra::MultiVectorFactory<SimpleFactory1,3,SMEVFactory,2> MVFactory2;
    typedef AMP::LinearAlgebra::MultiVectorFactory<MVFactory1,2,MVFactory2,2> MVFactory3;
#endif
#endif
// clang-format on

using AMP::LinearAlgebra::vectorTests;


int main( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    // Run the vector selector tests on different vectors
    vectorTests<SimpleFactory1>::testVectorSelector( &ut );
    vectorTests<SimpleFactory2>::testVectorSelector( &ut );

    vectorTests<ArrayFactory1>::testVectorSelector( &ut );
    vectorTests<ArrayFactory2>::testVectorSelector( &ut );

#if defined( USE_EXT_PETSC ) && defined( USE_EXT_TRILINOS )
    vectorTests<SNPVFactory>::testVectorSelector( &ut );
#endif
#ifdef USE_EXT_TRILINOS
#ifdef USE_TRILINOS_THYRA
    vectorTests<NativeThyraFactory>::testVectorSelector( &ut );
    vectorTests<ManagedThyraFactory<SimpleFactory2>>::testVectorSelector( &ut );
    vectorTests<ManagedNativeThyraFactory<SimpleFactory2>>::testVectorSelector( &ut );
#endif
    vectorTests<MVFactory1>::testVectorSelector( &ut );
    vectorTests<MVFactory2>::testVectorSelector( &ut );
    vectorTests<MVFactory3>::testVectorSelector( &ut );
    vectorTests<SMEVFactory>::testVectorSelector( &ut );
#ifdef USE_TRILINOS_THYRA
    vectorTests<NativeThyraFactory>( &ut );
    vectorTests<ManagedThyraFactory<SimpleFactory2>>::testVectorSelector( &ut );
    vectorTests<ManagedNativeThyraFactory<SimpleFactory2>>::testVectorSelector( &ut );
    vectorTests<ManagedNativeThyraFactory<MVFactory1>>::testVectorSelector( &ut );
#endif
#endif

// Run the tests on a subsetted vector
#ifdef USE_EXT_TRILINOS
    vectorTests<StridedVectorFactory<SMEVFactory>>::testManagedVector( &ut );
#else
    ut.expected_failure( "Compiled without trilinos" );
#endif

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
