#include "vectors/MultiVector.h"
#include "vectors/MultiVariable.h"

#include "test_VectorLoops.h"
#include "utils/AMPManager.h"

#ifdef USES_PETSC
    #include "vectors/petsc/ManagedPetscVector.h"
#endif
#ifdef USE_TRILINOS
    #include "vectors/trilinos/ManagedEpetraVector.h"
#endif
#ifdef USE_SUNDIALS
    #include "vectors/sundials/ManagedSundialsVector.h"
#endif

using namespace AMP::unit_test;

#if defined(USES_PETSC) && defined(USE_TRILINOS)
    typedef CloneFactory <ViewFactory <AMP::LinearAlgebra::PetscVector , SimpleManagedVectorFactory<AMP::LinearAlgebra::ManagedEpetraVector> > >         SMEVFactory;
    typedef CloneFactory <ViewFactory <AMP::LinearAlgebra::PetscVector , SimplePetscNativeFactory<AMP::LinearAlgebra::NativePetscVector> > >             SNPVFactory;
    typedef CloneFactory <ViewFactory <AMP::LinearAlgebra::PetscVector , MultiVectorFactory<SMEVFactory,1,SNPVFactory,1> > >              MVFactory1;
    typedef CloneFactory <ViewFactory <AMP::LinearAlgebra::PetscVector , MultiVectorFactory<SMEVFactory,3,SNPVFactory,2> > >              MVFactory2;
    typedef CloneFactory <ViewFactory <AMP::LinearAlgebra::PetscVector , MultiVectorFactory<MVFactory1,2,MVFactory2,2> > >                MVFactory3;
#endif


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

#if defined(USES_PETSC) && defined(USE_TRILINOS)
    std::cout << "Testing Iterator" << std::endl;
    VectorIteratorTests<MVFactory1> ( &ut );
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
#else
    ut.expected_failure("Compiled without petsc or trilinos");
#endif

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;

}
