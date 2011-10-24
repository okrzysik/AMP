#include <iostream>
#include "vectors/ManagedPetscVector.h"
#include "vectors/ManagedSundialsVector.h"
#include "vectors/trilinos/ManagedEpetraVector.h"
#include "vectors/DualVector.h"
#include "vectors/DualVariable.h"
#include "vectors/MultiVector.h"
#include "vectors/MultiVariable.h"
#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"

#include "test_VectorLoops.h"

using namespace AMP::unit_test;

typedef SimpleManagedVectorFactory<AMP::LinearAlgebra::ManagedEpetraVector>     SMEVFactory;
typedef SimplePetscNativeFactory<AMP::LinearAlgebra::NativePetscVector>         SNPVFactory;
typedef MultiVectorFactory<SMEVFactory,1,SNPVFactory,1>                         MVFactory1;
typedef MultiVectorFactory<SMEVFactory,3,SNPVFactory,2>                         MVFactory2;
typedef MultiVectorFactory<MVFactory1,2,MVFactory2,2>                           MVFactory3;




int main ( int argc , char **argv )
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
    
    std::cout << "Testing Iterator" << std::endl;
    VectorIteratorTests<MVFactory1>::run_test ( &ut );
    std::cout << std::endl;
    globalComm.barrier();

    std::cout << "Testing SimpleVector" << std::endl;
    testSimpleVector<15> ( &ut );
    testSimpleVector<45> ( &ut );
    std::cout << std::endl;
    globalComm.barrier();

    std::cout << "Testing ManagedEpetraVector" << std::endl;
    test_managed_vectors_loop<SMEVFactory> ( &ut );
    std::cout << std::endl;
    globalComm.barrier();

    std::cout << "Testing NativePetscVector" << std::endl;
    test_managed_vectors_loop<SNPVFactory> ( &ut );
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

    ut.report();

    int num_failed = ut.NumFailGlobal();
    if ( num_failed==0 )
        std::cout << "No errors detected" << std::endl;
    AMP::AMPManager::shutdown();
    return num_failed;

}


