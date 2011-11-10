#ifndef included_AMP_UnitTest_test_VectorLoops
#define included_AMP_UnitTest_test_VectorLoops

#include "vectors/DualVector.h"
#include "utils/UnitTest.h"

#ifdef USE_SUNDIALS
    #include "test_SundialsVectorTests.h"
    #include "vectors/sundials/ManagedSundialsVector.h"
    #include "vectors/sundials/SundialsVector.h"
#endif
#ifdef USE_PETSC
    #include "test_PetscVector.h"
    #include "vectors/petsc/ManagedPetscVector.h"
    #include "vectors/petsc/PetscVector.h"
#endif
#include "vectors/MultiVector.h"
#include "test_Vector.h"

/// \cond UNDOCUMENTED

using namespace AMP::unit_test;


template <class FACTORY>
void  test_petsc_bottom ( AMP::UnitTest *ut )
{
    InstantiatePetscVectors<FACTORY>::run_test ( ut );
    DuplicatePetscVector<FACTORY>::run_test ( ut );
    StaticCopyPetscVector<FACTORY>::run_test ( ut );
    StaticDuplicatePetscVector<FACTORY>::run_test ( ut );
    CopyPetscVector<FACTORY>::run_test ( ut );
    VerifyDotPetscVector<FACTORY>::run_test ( ut );
    VerifyNormsPetscVector<FACTORY>::run_test ( ut );
    VerifyScalePetscVector<FACTORY>::run_test ( ut );
    VerifyAXPYPetscVector<FACTORY>::run_test ( ut );
    VerifyAbsPetscVector<FACTORY>::run_test ( ut );
    VerifyMaxPointwiseDividePetscVector<FACTORY>::run_test ( ut );
    VerifyGetSizePetscVector<FACTORY>::run_test ( ut );
    VerifySwapPetscVector<FACTORY>::run_test ( ut );
    VerifyAXPBYPetscVector<FACTORY>::run_test ( ut );
    VerifySetPetscVector<FACTORY>::run_test ( ut );
    VerifySetRandomPetscVector<FACTORY>::run_test ( ut );
    VerifySqrtPetscVector<FACTORY>::run_test ( ut );
    VerifyPointwiseMultPetscVector<FACTORY>::run_test ( ut );
    VerifyPointwiseDividePetscVector<FACTORY>::run_test ( ut );
    VerifyPointwiseMaxPetscVector<FACTORY>::run_test ( ut );
    VerifyPointwiseMinPetscVector<FACTORY>::run_test ( ut );
    VerifyPointwiseMaxAbsPetscVector<FACTORY>::run_test ( ut );
    VerifyLogPetscVector<FACTORY>::run_test ( ut );
    VerifyExpPetscVector<FACTORY>::run_test ( ut );
    VerifyAYPXPetscVector<FACTORY>::run_test ( ut );
    VerifyAXPBYPCZPetscVector<FACTORY>::run_test ( ut );
}


template <class FACTORY>
void test_managed_vectors_bottom ( AMP::UnitTest *ut )
{
    DeepCloneOfView<FACTORY,AMP::LinearAlgebra::PetscVector>::run_test ( ut );
#ifdef USE_SUNDIALS
    DeepCloneOfView<FACTORY,AMP::LinearAlgebra::SundialsVector>::run_test ( ut );
#endif
    VectorIteratorLengthTest<FACTORY>::run_test ( ut );
    Bug_491<FACTORY>::run_test ( ut );
    Bug_728<FACTORY>::run_test ( ut );
    InstantiateVector<FACTORY>::run_test ( ut );
    SetToScalarVector<FACTORY>::run_test ( ut );
    SetRandomValuesVector<FACTORY>::run_test ( ut );
    CloneVector<FACTORY>::run_test ( ut );
    DotProductVector<FACTORY>::run_test ( ut );
    AbsVector<FACTORY>::run_test ( ut );
    L1NormVector<FACTORY>::run_test ( ut );
    L2NormVector<FACTORY>::run_test ( ut );
    MaxNormVector<FACTORY>::run_test ( ut );
    ScaleVector<FACTORY>::run_test ( ut );
    AddVector<FACTORY>::run_test ( ut );
    SubtractVector<FACTORY>::run_test ( ut );
    MultiplyVector<FACTORY>::run_test ( ut );
    DivideVector<FACTORY>::run_test ( ut );
    ReciprocalVector<FACTORY>::run_test ( ut );
    LinearSumVector<FACTORY>::run_test ( ut );
    AxpyVector<FACTORY>::run_test ( ut );
    AxpbyVector<FACTORY>::run_test ( ut );
    CopyVector<FACTORY>::run_test ( ut );
    VerifyVectorMin<FACTORY>::run_test ( ut );
    VerifyVectorMax<FACTORY>::run_test ( ut );
    VerifyVectorMaxMin<FACTORY>::run_test ( ut );
//    VectorIteratorTests<FACTORY>::run_test ( ut );
}


#ifdef USE_SUNDIALS
template <class FACTORY>
void test_sundials_bottom ( AMP::UnitTest *ut )
{
    CloneSundialsVector<FACTORY>::run_test ( ut );
    LinearSumSundialsVector<FACTORY>::run_test ( ut );
    ConstSundialsVector<FACTORY>::run_test ( ut );
    ProdSundialsVector<FACTORY>::run_test ( ut );
    DivSundialsVector<FACTORY>::run_test ( ut );
    ScaleSundialsVector<FACTORY>::run_test ( ut );
    AbsSundialsVector<FACTORY>::run_test ( ut );
    InvSundialsVector<FACTORY>::run_test ( ut );
    AddConstSundialsVector<FACTORY>::run_test ( ut );
    DotProdSundialsVector<FACTORY>::run_test ( ut );
    MaxNormSundialsVector<FACTORY>::run_test ( ut );
    WRMSNormSundialsVector<FACTORY>::run_test ( ut );
    L1NormSundialsVector<FACTORY>::run_test ( ut );
}
#endif


template <class FACTORY>
void test_managed_vectors_loop ( AMP::UnitTest *ut )
{
  test_managed_vectors_bottom<FACTORY> ( ut );

  typedef SimplePetscVectorFactory<FACTORY>   PETSC_FACTORY;
  test_petsc_bottom<PetscViewFactory<PETSC_FACTORY> > ( ut );
  test_petsc_bottom<PetscCloneFactory<PetscViewFactory<PETSC_FACTORY> > > ( ut );

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
    InstantiateVector<SimpleVectorFactory<I> >::run_test ( ut );
    SetToScalarVector<SimpleVectorFactory<I> >::run_test ( ut );
    CloneVector<SimpleVectorFactory<I> >::run_test ( ut );
    ScaleVector<SimpleVectorFactory<I> >::run_test ( ut );
    AddVector<SimpleVectorFactory<I> >::run_test ( ut );
    SubtractVector<SimpleVectorFactory<I> >::run_test ( ut );
    MultiplyVector<SimpleVectorFactory<I> >::run_test ( ut );
    DivideVector<SimpleVectorFactory<I> >::run_test ( ut );
    ReciprocalVector<SimpleVectorFactory<I> >::run_test ( ut );
    LinearSumVector<SimpleVectorFactory<I> >::run_test ( ut );
    AxpyVector<SimpleVectorFactory<I> >::run_test ( ut );
    AxpbyVector<SimpleVectorFactory<I> >::run_test ( ut );
    CopyVector<SimpleVectorFactory<I> >::run_test ( ut );
    VectorIteratorTests<SimpleVectorFactory<I> >::run_test ( ut );
    L1NormVector<SimpleVectorFactory<I> >::run_test ( ut );
    L2NormVector<SimpleVectorFactory<I> >::run_test ( ut );
    DotProductVector<SimpleVectorFactory<I> >::run_test ( ut );
    AbsVector<SimpleVectorFactory<I> >::run_test ( ut );
}


template <class FACTORY>
void test_parallel_vectors_loop ( AMP::UnitTest *ut )
{
    InstantiateVector<FACTORY>::run_test ( ut );
    VerifyVectorGhostCreate<FACTORY>::run_test ( ut );
    VerifyVectorMakeConsistentSet<FACTORY>::run_test ( ut );
    VerifyVectorMakeConsistentAdd<FACTORY>::run_test ( ut );
    CopyVectorConsistency<FACTORY>::run_test ( ut );
}

/// \endcond

#endif
