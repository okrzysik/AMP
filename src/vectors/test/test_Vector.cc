#include <iostream>
#include "vectors/MultiVector.h"
#include "vectors/MultiVariable.h"
#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"

#include "test_VectorLoops.h"
#include "test_Vector.h"

#ifdef USE_EXT_PETSC
    #include "petsc/NativePetscVector.h"
#endif
#ifdef USE_EXT_TRILINOS
    #include "trilinos/ManagedEpetraVector.h"
#endif


#ifdef USE_EXT_PETSC
    typedef AMP::unit_test::SimplePetscNativeFactory<AMP::LinearAlgebra::NativePetscVector>         SNPVFactory;
#endif
#ifdef USE_EXT_TRILINOS
    typedef AMP::unit_test::SimpleManagedVectorFactory<AMP::LinearAlgebra::ManagedEpetraVector>     SMEVFactory;
    typedef AMP::unit_test::MultiVectorFactory<SMEVFactory,1,SNPVFactory,1>                         MVFactory1;
    typedef AMP::unit_test::MultiVectorFactory<SMEVFactory,3,SNPVFactory,2>                         MVFactory2;
    typedef AMP::unit_test::MultiVectorFactory<MVFactory1,2,MVFactory2,2>                           MVFactory3;
#endif


using namespace AMP::unit_test;


int main ( int argc , char **argv )
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
    
    AMP::pout << "Testing NullVector" << std::endl;
    testNullVector( &ut );

    AMP::pout << "Testing SimpleVector" << std::endl;
    testBasicVector<SimpleVectorFactory<15,false> >( &ut );
    testBasicVector<SimpleVectorFactory<45,true> >( &ut );
    AMP::pout << std::endl;
    globalComm.barrier();

    #ifdef USE_EXT_PETSC
        AMP::pout << "Testing NativePetscVector" << std::endl;
        testManagedVector<SNPVFactory> ( &ut );
        AMP::pout << std::endl;
        globalComm.barrier();
    #else
        ut.expected_failure("Compiled without petsc");
    #endif

    #ifdef USE_EXT_TRILINOS

        AMP::pout << "Testing NativeThyraVector" << std::endl;
        testBasicVector<NativeThyraFactory>( &ut );
        AMP::pout << std::endl;
        globalComm.barrier();

        AMP::pout << "Testing ManagedThyraVector" << std::endl;
        testBasicVector<ManagedThyraFactory<SimpleVectorFactory<45,true> > >( &ut );
        AMP::pout << std::endl;
        globalComm.barrier();

        AMP::pout << "Testing NativeThyraVector of a ManagedThyraVector" << std::endl;
        testBasicVector<ManagedNativeThyraFactory<SimpleVectorFactory<45,true> > >( &ut );
        AMP::pout << std::endl;
        globalComm.barrier();

        AMP::pout << "Testing Iterator" << std::endl;
        VectorIteratorTests<MVFactory1> ( &ut );
        AMP::pout << std::endl;
        globalComm.barrier();

        AMP::pout << "Testing ManagedEpetraVector" << std::endl;
        testManagedVector<SMEVFactory> ( &ut );
        AMP::pout << std::endl;
        globalComm.barrier();

        AMP::pout << "Testing simple multivector" << std::endl;
        testManagedVector<MVFactory1> ( &ut );
        AMP::pout << std::endl;
        globalComm.barrier();

        AMP::pout << "Testing bigger multivector" << std::endl;
        testManagedVector<MVFactory2> ( &ut );
        AMP::pout << std::endl;
        globalComm.barrier();

        AMP::pout << "Testing multivector of multivector" << std::endl;
        testManagedVector<MVFactory3> ( &ut );
        AMP::pout << std::endl;
        globalComm.barrier();

        AMP::pout << "Testing NativeThyraVector of a ManagedThyraVector of a MultVector" << std::endl;
        testBasicVector<ManagedNativeThyraFactory<MVFactory1> >( &ut );
        AMP::pout << std::endl;
        globalComm.barrier();

    #else
        ut.expected_failure("Compiled without trilinos");
    #endif

    ut.report();

    int num_failed = ut.NumFailGlobal();
    if ( num_failed==0 )
        AMP::pout << "No errors detected" << std::endl;
    AMP::AMPManager::shutdown();
    return num_failed;

}


