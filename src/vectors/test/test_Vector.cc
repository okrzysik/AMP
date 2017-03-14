#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "vectors/MultiVariable.h"
#include "vectors/MultiVector.h"
#include <iostream>

#include "vectors/testHelpers/VectorTests.h"
#include "vectors/testHelpers/testVectorFactory.h"

#ifdef USE_EXT_PETSC
#include "vectors/petsc/NativePetscVector.h"
#endif
#ifdef USE_EXT_TRILINOS
#include "vectors/trilinos/ManagedEpetraVector.h"
#endif

#include "test_ArrayVector.h"


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
    typedef AMP::LinearAlgebra::MultiVectorFactory<SMEVFactory, 1, SNPVFactory, 1> MVFactory1;
    typedef AMP::LinearAlgebra::MultiVectorFactory<SMEVFactory, 3, SNPVFactory, 2> MVFactory2;
    typedef AMP::LinearAlgebra::MultiVectorFactory<MVFactory1, 2, MVFactory2, 2> MVFactory3;
#else
    typedef AMP::LinearAlgebra::MultiVectorFactory<SimpleVectorFactory<15,false>,1,SMEVFactory,1> MVFactory1;
    typedef AMP::LinearAlgebra::MultiVectorFactory<SimpleVectorFactory<15,false>,3,SMEVFactory,2> MVFactory2;
    typedef AMP::LinearAlgebra::MultiVectorFactory<MVFactory1, 2, MVFactory2, 2> MVFactory3;
#endif
#endif
// clang-format on


using namespace AMP::unit_test;
using AMP::LinearAlgebra::vectorTests;


int main( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    AMP::pout << "Testing NullVector" << std::endl;
    vectorTests<SimpleFactory1>::testNullVector( &ut );

    AMP::pout << "Testing SimpleVector" << std::endl;
    vectorTests<SimpleFactory1>::testBasicVector( &ut );
    vectorTests<SimpleFactory2>::testBasicVector( &ut );
    AMP::pout << std::endl;
    globalComm.barrier();

    AMP::pout << "Testing ArrayVector" << std::endl;
    vectorTests<ArrayFactory1>::testBasicVector( &ut );
    vectorTests<ArrayFactory2>::testBasicVector( &ut );
    std::vector<size_t> dims{ 3, 3, 3, 3 };
    testArrayVectorDimensions<double>( dims, ut );
    AMP::pout << std::endl;
    globalComm.barrier();

#if defined( USE_EXT_PETSC ) && defined( USE_EXT_TRILINOS )
    AMP::pout << "Testing NativePetscVector" << std::endl;
    vectorTests<SNPVFactory>::testManagedVector( &ut );
    AMP::pout << std::endl;
    globalComm.barrier();
#else
    ut.expected_failure( "Compiled without petsc" );
#endif

#ifdef USE_EXT_TRILINOS

#ifdef USE_TRILINOS_THYRA
    AMP::pout << "Testing NativeThyraVector" << std::endl;
    vectorTests<NativeThyraFactory>::testBasicVector( &ut );
    AMP::pout << std::endl;
    globalComm.barrier();

    AMP::pout << "Testing ManagedThyraVector" << std::endl;
    vectorTests<ManagedThyraFactory<SimpleFactory2>>::testBasicVector( &ut );
    AMP::pout << std::endl;
    globalComm.barrier();

    AMP::pout << "Testing NativeThyraVector of a ManagedThyraVector" << std::endl;
    vectorTests<ManagedNativeThyraFactory<SimpleFactory2>>::testBasicVector( &ut );
    AMP::pout << std::endl;
    globalComm.barrier();
#endif

    AMP::pout << "Testing Iterator" << std::endl;
    vectorTests<MVFactory1>::VectorIteratorTests( &ut );
    AMP::pout << std::endl;
    globalComm.barrier();

    AMP::pout << "Testing ManagedEpetraVector" << std::endl;
    vectorTests<SMEVFactory>::testManagedVector( &ut );
    AMP::pout << std::endl;
    globalComm.barrier();

    AMP::pout << "Testing simple multivector" << std::endl;
    vectorTests<MVFactory1>::testManagedVector( &ut );
    AMP::pout << std::endl;
    globalComm.barrier();

    AMP::pout << "Testing bigger multivector" << std::endl;
    vectorTests<MVFactory2>::testManagedVector( &ut );
    AMP::pout << std::endl;
    globalComm.barrier();

    AMP::pout << "Testing multivector of multivector" << std::endl;
    vectorTests<MVFactory3>::testManagedVector( &ut );
    AMP::pout << std::endl;
    globalComm.barrier();

#ifdef USE_TRILINOS_THYRA
    AMP::pout << "Testing NativeThyraVector of a ManagedThyraVector of a MultVector" << std::endl;
    vectorTests<ManagedNativeThyraFactory<MVFactory1>>::testBasicVector( &ut );
    AMP::pout << std::endl;
    globalComm.barrier();
#endif

// Run Belos tests of thyra vectors
#ifdef USE_TRILINOS_THYRA
    AMP::pout << "Testing Belos interface to Thyra vectors" << std::endl;
    testBelosThyraVector<NativeThyraFactory>( &ut );
    testBelosThyraVector<ManagedThyraFactory<SimpleVectorFactory<45, true>>>( &ut );
    testBelosThyraVector<ManagedNativeThyraFactory<SimpleVectorFactory<45, true>>>( &ut );
    testBelosThyraVector<ManagedNativeThyraFactory<MVFactory1>>( &ut );
    AMP::pout << std::endl;
    globalComm.barrier();
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
