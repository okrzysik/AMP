#include "AMP/utils/AMPManager.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/VectorSelector.h"

#include "test_VectorHelpers.h"


int main( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    // Run the vector selector tests on different vectors
    testVectorSelector( ut, SimpleFactory1 );
    testVectorSelector( ut, SimpleFactory2 );

    testVectorSelector( ut, ArrayFactory1 );
    testVectorSelector( ut, ArrayFactory2 );

#if defined( USE_EXT_PETSC ) && defined( USE_EXT_TRILINOS )
    testVectorSelector( ut, SNPVFactory );
#endif
#ifdef USE_EXT_TRILINOS
#ifdef USE_TRILINOS_THYRA
    testVectorSelector( ut, NTFactory );
    testVectorSelector( ut, MTFactory );
    testVectorSelector( ut, MNTFactory );
    testVectorSelector( ut, MNT_MVFactory );
#endif
    testVectorSelector( ut, MVFactory1 );
    testVectorSelector( ut, MVFactory2 );
    testVectorSelector( ut, MVFactory3 );
    testVectorSelector( ut, SMEVFactory );
#endif

// Run the tests on a subsetted vector
#ifdef USE_EXT_TRILINOS
    testManagedVector( ut, StridedFactory );
#else
    ut.expected_failure( "Compiled without trilinos" );
#endif

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
