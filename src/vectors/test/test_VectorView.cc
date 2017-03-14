#include "vectors/MultiVector.h"
#include "vectors/MultiVariable.h"

#include "vectors/testHelpers/VectorTests.h"
#include "vectors/testHelpers/testVectorFactory.h"

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


// Typedef some factories
// clang-format off
typedef AMP::LinearAlgebra::SimpleVectorFactory<15,false,double>  SimpleFactory1;
typedef AMP::LinearAlgebra::SimpleVectorFactory<45, true,double>  SimpleFactory2;
typedef AMP::LinearAlgebra::ArrayVectorFactory<4,10,false,double> ArrayFactory1;
typedef AMP::LinearAlgebra::ArrayVectorFactory<4,10, true,double> ArrayFactory2;
#if defined( USE_EXT_PETSC ) && defined( USE_EXT_TRILINOS )
    typedef AMP::LinearAlgebra::ViewFactory<AMP::LinearAlgebra::PetscVector,AMP::LinearAlgebra::SimpleManagedVectorFactory<AMP::LinearAlgebra::ManagedEpetraVector>> SMEVFactory;
    typedef AMP::LinearAlgebra::ViewFactory<AMP::LinearAlgebra::PetscVector,AMP::LinearAlgebra::SimplePetscNativeFactory<AMP::LinearAlgebra::NativePetscVector>>     SNPVFactory;
    typedef AMP::LinearAlgebra::ViewFactory<AMP::LinearAlgebra::PetscVector,AMP::LinearAlgebra::MultiVectorFactory<SMEVFactory,1,SNPVFactory,1>> MVFactory1;
    typedef AMP::LinearAlgebra::ViewFactory<AMP::LinearAlgebra::PetscVector,AMP::LinearAlgebra::MultiVectorFactory<SMEVFactory,3,SNPVFactory,2>> MVFactory2;
    typedef AMP::LinearAlgebra::ViewFactory<AMP::LinearAlgebra::PetscVector,AMP::LinearAlgebra::MultiVectorFactory<MVFactory1,2,MVFactory2,2>>   MVFactory3;
#endif
// clang-format on

using AMP::LinearAlgebra::vectorTests;


int main( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    AMP::pout << "Testing SimpleVector" << std::endl;
    vectorTests<SimpleFactory1>::testBasicVector( &ut );
    vectorTests<SimpleFactory2>::testBasicVector( &ut );
    AMP::pout << std::endl;
    globalComm.barrier();

    AMP::pout << "Testing ArrayVector" << std::endl;
    vectorTests<ArrayFactory1>::testBasicVector( &ut );
    vectorTests<ArrayFactory2>::testBasicVector( &ut );
    AMP::pout << std::endl;
    globalComm.barrier();

#if defined( USE_EXT_PETSC ) && defined( USE_EXT_TRILINOS )
    AMP::pout << "Testing Iterator" << std::endl;
    vectorTests<MVFactory1>::VectorIteratorTests( &ut );
    AMP::pout << std::endl;
    globalComm.barrier();

    AMP::pout << "Testing ManagedEpetraVector" << std::endl;
    vectorTests<SMEVFactory>::VectorIteratorTests( &ut );
    AMP::pout << std::endl;
    globalComm.barrier();

    AMP::pout << "Testing NativePetscVector" << std::endl;
    vectorTests<SNPVFactory>::VectorIteratorTests( &ut );
    AMP::pout << std::endl;
    globalComm.barrier();

    AMP::pout << "Testing simple multivector" << std::endl;
    vectorTests<MVFactory1>::VectorIteratorTests( &ut );
    AMP::pout << std::endl;
    globalComm.barrier();

    AMP::pout << "Testing bigger multivector" << std::endl;
    vectorTests<MVFactory2>::VectorIteratorTests( &ut );
    AMP::pout << std::endl;
    globalComm.barrier();

    AMP::pout << "Testing multivector of multivector" << std::endl;
    vectorTests<MVFactory3>::VectorIteratorTests( &ut );
    AMP::pout << std::endl;
    globalComm.barrier();
#else
    ut.expected_failure( "Compiled without petsc or trilinos" );
#endif

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
