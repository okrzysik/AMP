#ifndef included_AMP_UnitTest_test_VectorLoops
#define included_AMP_UnitTest_test_VectorLoops

#include "utils/UnitTest.h"

#include "test_Vector.h"
#include "test_VectorTests.h"
#ifdef USE_SUNDIALS
    #include "test_SundialsVectorTests.h"
#endif
#ifdef USE_PETSC
    #include "test_PetscVector.h"
#endif

/// \cond UNDOCUMENTED

using namespace AMP::unit_test;


#ifdef USE_PETSC
template <class FACTORY>
void  test_petsc_bottom ( AMP::UnitTest *ut )
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
void test_managed_vectors_bottom ( AMP::UnitTest *ut )
{
  #ifdef USE_PETSC
    DeepCloneOfView<FACTORY,AMP::LinearAlgebra::PetscVector>( ut );
    Bug_491<FACTORY>( ut );
  #endif
  #ifdef USE_SUNDIALS
    DeepCloneOfView<FACTORY,AMP::LinearAlgebra::SundialsVector>( ut );
  #endif
    VectorIteratorLengthTest<FACTORY>( ut );
    Bug_728<FACTORY>( ut );
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
//    VectorIteratorTests<FACTORY>( ut );
}


#ifdef USE_SUNDIALS
template <class FACTORY>
void test_sundials_bottom ( AMP::UnitTest *ut )
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
void test_managed_vectors_loop ( AMP::UnitTest *ut )
{
  test_managed_vectors_bottom<FACTORY> ( ut );

  #ifdef USE_PETSC
    typedef SimplePetscVectorFactory<FACTORY>   PETSC_FACTORY;
    test_petsc_bottom<PetscViewFactory<PETSC_FACTORY> > ( ut );
    test_petsc_bottom<PetscCloneFactory<PetscViewFactory<PETSC_FACTORY> > > ( ut );
  #endif

  #ifdef USE_SUNDIALS
    test_managed_vectors_bottom<ViewFactory<AMP::LinearAlgebra::SundialsVector , FACTORY> >( ut );
    test_managed_vectors_bottom<CloneFactory<ViewFactory<AMP::LinearAlgebra::SundialsVector , FACTORY> > > ( ut );
    test_sundials_bottom<ViewFactory<AMP::LinearAlgebra::SundialsVector , FACTORY> >( ut );
    test_sundials_bottom<CloneFactory<ViewFactory<AMP::LinearAlgebra::SundialsVector , FACTORY> > > ( ut );
  #endif

}


template <int I>
void testSimpleVector ( AMP::UnitTest *ut )
{
    InstantiateVector<SimpleVectorFactory<I> >( ut );
    SetToScalarVector<SimpleVectorFactory<I> >( ut );
    CloneVector<SimpleVectorFactory<I> >( ut );
    ScaleVector<SimpleVectorFactory<I> >( ut );
    AddVector<SimpleVectorFactory<I> >( ut );
    SubtractVector<SimpleVectorFactory<I> >( ut );
    MultiplyVector<SimpleVectorFactory<I> >( ut );
    DivideVector<SimpleVectorFactory<I> >( ut );
    ReciprocalVector<SimpleVectorFactory<I> >( ut );
    LinearSumVector<SimpleVectorFactory<I> >::run_test( ut );
    AxpyVector<SimpleVectorFactory<I> >::run_test( ut );
    AxpbyVector<SimpleVectorFactory<I> >::run_test( ut );
    CopyVector<SimpleVectorFactory<I> >::run_test( ut );
    VectorIteratorTests<SimpleVectorFactory<I> >( ut );
    L1NormVector<SimpleVectorFactory<I> >( ut );
    L2NormVector<SimpleVectorFactory<I> >( ut );
    DotProductVector<SimpleVectorFactory<I> >( ut );
    AbsVector<SimpleVectorFactory<I> >( ut );
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

/// \endcond

#endif
