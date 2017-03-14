#ifndef included_AMP_UnitTest_test_VectorLoops
#define included_AMP_UnitTest_test_VectorLoops

#include "vectors/testHelpers/VectorTests.h"
#include "vectors/testHelpers/test_VectorSelectorTests.h"

#include "utils/UnitTest.h"

#ifdef USE_EXT_SUNDIALS
#include "vectors/testHelpers/SundialsVectorTests.inline.h"
#endif
#ifdef USE_EXT_PETSC
#include "vectors/testHelpers/PetscVector.inline.h"
#include "vectors/testHelpers/PetscVectorTests.inline.h"
#endif


namespace AMP {
namespace LinearAlgebra {


template <class FACTORY>
void vectorTests<FACTORY>::testPetscVector( AMP::UnitTest *ut )
{
#ifdef USE_EXT_PETSC
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
#endif
}


template <class FACTORY>
void vectorTests<FACTORY>::testBasicVector( AMP::UnitTest *ut )
{
    InstantiateVector( ut );
    SetToScalarVector( ut );
    SetRandomValuesVector::run_test( ut );
    CloneVector( ut );
    DotProductVector( ut );
    AbsVector( ut );
    L1NormVector( ut );
    L2NormVector( ut );
    MaxNormVector( ut );
    ScaleVector( ut );
    AddVector( ut );
    SubtractVector( ut );
    MultiplyVector( ut );
    DivideVector( ut );
    ReciprocalVector( ut );
    LinearSumVector::run_test( ut );
    AxpyVector::run_test( ut );
    AxpbyVector::run_test( ut );
    CopyVector::run_test( ut );
    CopyRawDataBlockVector::run_test( ut );
    VerifyVectorMin( ut );
    VerifyVectorMax( ut );
    VerifyVectorMaxMin( ut );
#ifdef USE_EXT_PETSC
    DeepCloneOfView<AMP::LinearAlgebra::PetscVector>( ut );
    Bug_491( ut );
#endif
#ifdef USE_EXT_SUNDIALS
    DeepCloneOfView<AMP::LinearAlgebra::SundialsVector>( ut );
#endif
    VectorIteratorLengthTest( ut );
    Bug_728( ut );
    //    VectorIteratorTests( ut );
    TestMultivectorDuplicate( ut );
}


template <class FACTORY>
void vectorTests<FACTORY>::testSundialsVector( AMP::UnitTest *ut )
{
#ifdef USE_EXT_SUNDIALS
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
#endif
}


template <class FACTORY>
void vectorTests<FACTORY>::testManagedVector( AMP::UnitTest *ut )
{
    testBasicVector( ut );

#ifdef USE_EXT_PETSC
    typedef SimplePetscVectorFactory<FACTORY> PETSC_FACTORY;
    vectorTests<PetscViewFactory<PETSC_FACTORY>>::testPetscVector( ut );
    vectorTests<PetscCloneFactory<PetscViewFactory<PETSC_FACTORY>>>::testPetscVector( ut );
#endif

#ifdef USE_EXT_SUNDIALS
    vectorTests<ViewFactory<AMP::LinearAlgebra::SundialsVector, FACTORY>>::testBasicVector( ut );
    vectorTests<CloneFactory<ViewFactory<AMP::LinearAlgebra::SundialsVector, FACTORY>>>::testBasicVector( ut );
    vectorTests<ViewFactory<AMP::LinearAlgebra::SundialsVector, FACTORY>>::testSundialsVector( ut );
    vectorTests<CloneFactory<ViewFactory<AMP::LinearAlgebra::SundialsVector, FACTORY>>>::testSundialsVector( ut );
#endif
}


template <class FACTORY>
void vectorTests<FACTORY>::testNullVector( AMP::UnitTest *ut )
{
    vectorTests<NullVectorFactory>::InstantiateVector( ut );
}


template <class FACTORY>
void vectorTests<FACTORY>::testParallelVectors( AMP::UnitTest *ut )
{
    InstantiateVector( ut );
    // VerifyVectorGhostCreate( ut );
    VerifyVectorMakeConsistentSet( ut );
    VerifyVectorMakeConsistentAdd( ut );
    CopyVectorConsistency( ut );
}


template <class FACTORY>
void vectorTests<FACTORY>::testVectorSelector( AMP::UnitTest *ut )
{
    testAllSelectors<FACTORY>( ut );
    test_VS_ByVariableName<FACTORY>( ut );
    test_VS_Comm<FACTORY>( ut );
}



} // namespace LinearAlgebra
} // namespace AMP

#endif
