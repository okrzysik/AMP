#include <iostream>
#include "vectors/MultiVector.h"
#include "vectors/MultiVariable.h"
#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"

#include "test_VectorLoops.h"

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
    
    std::cout << "Testing SimpleVector" << std::endl;
    testSimpleVector<15> ( &ut );
    testSimpleVector<45> ( &ut );
    std::cout << std::endl;
    globalComm.barrier();

#ifdef USE_EXT_PETSC
    std::cout << "Testing NativePetscVector" << std::endl;
    test_managed_vectors_loop<SNPVFactory> ( &ut );
    std::cout << std::endl;
    globalComm.barrier();
#else
    ut.expected_failure("Compiled without petsc");
#endif

#ifdef USE_EXT_TRILINOS
    std::cout << "Testing Iterator" << std::endl;
    VectorIteratorTests<MVFactory1> ( &ut );
    std::cout << std::endl;
    globalComm.barrier();

    std::cout << "Testing ManagedEpetraVector" << std::endl;
    test_managed_vectors_loop<SMEVFactory> ( &ut );
    std::cout << std::endl;
    globalComm.barrier();

    std::cout << "Testing simple multivector" << std::endl;
    test_managed_vectors_loop<MVFactory1> ( &ut );
    std::cout << std::endl;
    globalComm.barrier();

    std::cout << "Testing bigger multivector" << std::endl;
    test_managed_vectors_loop<MVFactory2> ( &ut );
    std::cout << std::endl;
    globalComm.barrier();

    std::cout << "Testing multivector of multivector" << std::endl;
    test_managed_vectors_loop<MVFactory3> ( &ut );
    std::cout << std::endl;
    globalComm.barrier();
#else
    ut.expected_failure("Compiled without epetra");
#endif

    ut.report();

    int num_failed = ut.NumFailGlobal();
    if ( num_failed==0 )
        std::cout << "No errors detected" << std::endl;
    AMP::AMPManager::shutdown();
    return num_failed;

}


