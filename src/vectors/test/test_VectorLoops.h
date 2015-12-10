#ifndef included_AMP_UnitTest_test_VectorLoops
#define included_AMP_UnitTest_test_VectorLoops

#include "utils/UnitTest.h"

#include "test_Vector.h"
#include "test_VectorTests.h"
#include "test_VectorSelectorTests.h"
#ifdef USE_EXT_SUNDIALS
    #include "test_SundialsVectorTests.h"
#endif
#ifdef USE_EXT_PETSC
    #include "test_PetscVector.h"
#endif

/// \cond UNDOCUMENTED

using namespace AMP::unit_test;


#ifdef USE_EXT_PETSC
template <class FACTORY>
void  testPetscVector( AMP::UnitTest *ut )
{
    InstantiatePetscVectors<FACTORY>::run_test( ut );
    DuplicatePetscVector<FACTORY>::run_test( ut );
    StaticCopyPetscVector<FACTORY>::run_test( ut );
    StaticDuplicatePetscVector<FACTORY>::run_test( ut );
    CopyPetscVector<FACTORY>::run_test( ut );
    VerifyDotPetscVector<FACTORY>::run_test( ut );
    VerifyNormsPetscVector<FACTORY>::run_test( ut );
    VerifyScalePetscVector<FACTORY>::run_test( ut );
    VerifyAXPYPetscVector<FACTORY>::run_test( ut );
    VerifyAbsPetscVector<FACTORY>::run_test( ut );
    VerifyMaxPointwiseDividePetscVector<FACTORY>::run_test( ut );
    VerifyGetSizePetscVector<FACTORY>::run_test( ut );
    VerifySwapPetscVector<FACTORY>::run_test( ut );
    VerifyAXPBYPetscVector<FACTORY>::run_test( ut );
    VerifySetPetscVector<FACTORY>::run_test( ut );
    VerifySetRandomPetscVector<FACTORY>::run_test( ut );
    VerifySqrtPetscVector<FACTORY>::run_test( ut );
    VerifyPointwiseMultPetscVector<FACTORY>::run_test( ut );
    VerifyPointwiseDividePetscVector<FACTORY>::run_test( ut );
    VerifyPointwiseMaxPetscVector<FACTORY>::run_test( ut );
    VerifyPointwiseMinPetscVector<FACTORY>::run_test( ut );
    VerifyPointwiseMaxAbsPetscVector<FACTORY>::run_test( ut );
    VerifyLogPetscVector<FACTORY>::run_test( ut );
    VerifyExpPetscVector<FACTORY>::run_test( ut );
    VerifyAYPXPetscVector<FACTORY>::run_test( ut );
    VerifyAXPBYPCZPetscVector<FACTORY>::run_test( ut );
}
#endif


template <class FACTORY>
void testBasicVector( AMP::UnitTest *ut )
{
    InstantiateVector<FACTORY>( ut );
    SetToScalarVector<FACTORY>( ut );
    SetRandomValuesVector<FACTORY>::run_test( ut );
    CloneVector<FACTORY>( ut );
    DotProductVector<FACTORY>( ut );
    AbsVector<FACTORY>( ut );
    L1NormVector<FACTORY>( ut );
    L2NormVector<FACTORY>( ut );
    MaxNormVector<FACTORY>( ut );
    ScaleVector<FACTORY>( ut );
    AddVector<FACTORY>( ut );
    SubtractVector<FACTORY>( ut );
    MultiplyVector<FACTORY>( ut );
    DivideVector<FACTORY>( ut );
    ReciprocalVector<FACTORY>( ut );
    LinearSumVector<FACTORY>::run_test( ut );
    AxpyVector<FACTORY>::run_test( ut );
    AxpbyVector<FACTORY>::run_test( ut );
    CopyVector<FACTORY>::run_test( ut );
    VerifyVectorMin<FACTORY>( ut );
    VerifyVectorMax<FACTORY>( ut );
    VerifyVectorMaxMin<FACTORY>( ut );
    #ifdef USE_EXT_PETSC
        DeepCloneOfView<FACTORY,AMP::LinearAlgebra::PetscVector>( ut );
        Bug_491<FACTORY>( ut );
    #endif
    #ifdef USE_EXT_SUNDIALS
        DeepCloneOfView<FACTORY,AMP::LinearAlgebra::SundialsVector>( ut );
    #endif
    VectorIteratorLengthTest<FACTORY>( ut );
    Bug_728<FACTORY>( ut );
//    VectorIteratorTests<FACTORY>( ut );
    TestMultivectorDuplicate<FACTORY>( ut );
}


#ifdef USE_EXT_SUNDIALS
template <class FACTORY>
void testSundialsVector( AMP::UnitTest *ut )
{
    CloneSundialsVector<FACTORY>::run_test( ut );
    LinearSumSundialsVector<FACTORY>::run_test( ut );
    ConstSundialsVector<FACTORY>::run_test( ut );
    ProdSundialsVector<FACTORY>::run_test( ut );
    DivSundialsVector<FACTORY>::run_test( ut );
    ScaleSundialsVector<FACTORY>::run_test( ut );
    AbsSundialsVector<FACTORY>::run_test( ut );
    InvSundialsVector<FACTORY>::run_test( ut );
    AddConstSundialsVector<FACTORY>::run_test( ut );
    DotProdSundialsVector<FACTORY>::run_test( ut );
    MaxNormSundialsVector<FACTORY>::run_test( ut );
    WRMSNormSundialsVector<FACTORY>::run_test( ut );
    L1NormSundialsVector<FACTORY>::run_test( ut );
}
#endif


template <class FACTORY>
void testManagedVector( AMP::UnitTest *ut )
{
  testBasicVector<FACTORY> ( ut );

  #ifdef USE_EXT_PETSC
    typedef SimplePetscVectorFactory<FACTORY>   PETSC_FACTORY;
    testPetscVector<PetscViewFactory<PETSC_FACTORY> > ( ut );
    testPetscVector<PetscCloneFactory<PetscViewFactory<PETSC_FACTORY> > > ( ut );
  #endif

  #ifdef USE_EXT_SUNDIALS
    testBasicVector<ViewFactory<AMP::LinearAlgebra::SundialsVector , FACTORY> >( ut );
    testBasicVector<CloneFactory<ViewFactory<AMP::LinearAlgebra::SundialsVector , FACTORY> > > ( ut );
    testSundialsVector<ViewFactory<AMP::LinearAlgebra::SundialsVector , FACTORY> >( ut );
    testSundialsVector<CloneFactory<ViewFactory<AMP::LinearAlgebra::SundialsVector , FACTORY> > > ( ut );
  #endif

}


void testNullVector ( AMP::UnitTest *ut )
{
    InstantiateVector<NullVectorFactory>( ut );
}


template <class FACTORY>
void test_parallel_vectors_loop ( AMP::UnitTest *ut )
{
    InstantiateVector<FACTORY>( ut );
    //VerifyVectorGhostCreate<FACTORY>( ut );
    VerifyVectorMakeConsistentSet<FACTORY>( ut );
    VerifyVectorMakeConsistentAdd<FACTORY>( ut );
    CopyVectorConsistency<FACTORY>( ut );
}


template <class FACTORY>
void test_vector_selector_loop ( AMP::UnitTest *ut )
{
    testAllSelectors<FACTORY>( ut );
    test_VS_ByVariableName<FACTORY>( ut );
    test_VS_Comm<FACTORY>( ut );
}


/// \endcond

#endif
