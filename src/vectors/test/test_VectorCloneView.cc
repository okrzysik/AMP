#include "vectors/MultiVector.h"
#include "vectors/MultiVariable.h"

#include "test_VectorLoops.h"
#include "utils/AMPManager.h"

#ifdef USE_EXT_PETSC
    #include "vectors/petsc/ManagedPetscVector.h"
#endif
#ifdef USE_EXT_TRILINOS
    #include "vectors/trilinos/ManagedEpetraVector.h"
#endif
#ifdef USE_EXT_SUNDIALS
    #include "vectors/sundials/ManagedSundialsVector.h"
#endif

using namespace AMP::unit_test;

#if defined(USE_EXT_PETSC) && defined(USE_EXT_TRILINOS)
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

    AMP::pout << "Testing SimpleVector" << std::endl;
    testBasicVector<SimpleVectorFactory<15,false> >( &ut );
    testBasicVector<SimpleVectorFactory<45,true> >( &ut );
    AMP::pout << std::endl;
    globalComm.barrier();

    AMP::pout << "Testing ArrayVector" << std::endl;
    testBasicVector<ArrayVectorFactory<4,15,false> >( &ut );
    testBasicVector<ArrayVectorFactory<4,45,true> >( &ut );
    AMP::pout << std::endl;
    globalComm.barrier();

#if defined(USE_EXT_PETSC) && defined(USE_EXT_TRILINOS)
    AMP::pout << "Testing Iterator" << std::endl;
    VectorIteratorTests<MVFactory1> ( &ut );
    AMP::pout << std::endl;
    globalComm.barrier();
  
    AMP::pout << "Testing ManagedEpetraVector" << std::endl;
    testManagedVector<SMEVFactory> ( &ut );
    AMP::pout << std::endl;
    globalComm.barrier();

    AMP::pout << "Testing NativePetscVector" << std::endl;
    testManagedVector<SNPVFactory> ( &ut );
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
#else
    ut.expected_failure("Compiled without petsc or trilinos");
#endif

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;

}
