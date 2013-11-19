#include "vectors/VectorSelector.h"
#include "vectors/MultiVector.h"
#include "vectors/VectorSelector.h"

#include "utils/AMPManager.h"
#include "utils/UnitTest.h"

#ifdef USE_EXT_TRILINOS
    #include "vectors/trilinos/ManagedEpetraVector.h"
#endif

#include "test_VectorLoops.h"
#include "test_VectorSelector.h"

using namespace AMP::unit_test;


#ifdef USE_EXT_PETSC
    typedef AMP::unit_test::SimplePetscNativeFactory<AMP::LinearAlgebra::NativePetscVector>         SNPVFactory;
#endif
#ifdef USE_EXT_TRILINOS
    typedef AMP::unit_test::SimpleManagedVectorFactory<AMP::LinearAlgebra::ManagedEpetraVector>     SMEVFactory;
    #ifdef USE_EXT_PETSC
        typedef AMP::unit_test::MultiVectorFactory<SMEVFactory,1,SNPVFactory,1>                     MVFactory1;
        typedef AMP::unit_test::MultiVectorFactory<SMEVFactory,3,SNPVFactory,2>                     MVFactory2;
        typedef AMP::unit_test::MultiVectorFactory<MVFactory1,2,MVFactory2,2>                       MVFactory3;
    #else
        typedef AMP::unit_test::MultiVectorFactory<SimpleVectorFactory<15,false>,1,SMEVFactory,1>   MVFactory1;
        typedef AMP::unit_test::MultiVectorFactory<SimpleVectorFactory<15,false>,3,SMEVFactory,2>   MVFactory2;
        typedef AMP::unit_test::MultiVectorFactory<MVFactory1,2,MVFactory2,2>                       MVFactory3;
    #endif
#endif



int main ( int argc , char **argv )
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    // Run the vector selector tests on different vectors
    test_vector_selector_loop<SimpleVectorFactory<15,false> >( &ut );
    test_vector_selector_loop<SimpleVectorFactory<45,true> >( &ut );
    #ifdef USE_EXT_PETSC
        test_vector_selector_loop<SNPVFactory>( &ut );
    #endif
    #ifdef USE_EXT_TRILINOS
        #ifdef USE_TRILINOS_THYRA
            test_vector_selector_loop<NativeThyraFactory>( &ut );
            test_vector_selector_loop<ManagedThyraFactory<SimpleVectorFactory<45,true> > >( &ut );
            test_vector_selector_loop<ManagedNativeThyraFactory<SimpleVectorFactory<45,true> > >( &ut );
        #endif
        test_vector_selector_loop<MVFactory1>( &ut );
        test_vector_selector_loop<MVFactory2>( &ut );
        test_vector_selector_loop<MVFactory3>( &ut );
        test_vector_selector_loop<SMEVFactory>( &ut );
        #ifdef USE_TRILINOS_THYRA
            test_vector_selector_loop<NativeThyraFactory>( &ut );
            test_vector_selector_loop<ManagedThyraFactory<SimpleVectorFactory<45,true> > >( &ut );
            test_vector_selector_loop<ManagedNativeThyraFactory<SimpleVectorFactory<45,true> > >( &ut );
            test_vector_selector_loop<ManagedNativeThyraFactory<MVFactory1> >( &ut );
        #endif
    #endif

    // Run the tests on a subsetted vector
    #ifdef USE_EXT_TRILINOS
        testManagedVector<StridedVectorFactory<SMEVFactory> > ( &ut );
    #else
        ut.expected_failure("Compiled without trilinos");
    #endif

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}


