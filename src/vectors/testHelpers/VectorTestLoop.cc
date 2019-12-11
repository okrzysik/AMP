#include "AMP/vectors/testHelpers/VectorTests.h"

#include "AMP/utils/UnitTest.h"

#include "AMP/vectors/testHelpers/VectorFactory.h"

#ifdef USE_EXT_SUNDIALS
#include "AMP/vectors/sundials/SundialsVector.h"
#include "AMP/vectors/testHelpers/sundials/SundialsVectorTests.h"
#endif
#ifdef USE_EXT_PETSC
#include "AMP/vectors/petsc/PetscVector.h"
#include "AMP/vectors/testHelpers/petsc/PetscVectorFactory.h"
#include "AMP/vectors/testHelpers/petsc/PetscVectorTests.h"
#endif


namespace AMP {
namespace LinearAlgebra {


void VectorTests::testBasicVector( AMP::UnitTest *ut )
{
    InstantiateVector( ut );
    SetToScalarVector( ut );
    SetRandomValuesVector( ut );
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
    LinearSumVector( ut );
    AxpyVector( ut );
    AxpbyVector( ut );
    CopyVector( ut );
    CopyRawDataBlockVector( ut );
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


void VectorTests::testManagedVector( AMP::UnitTest *ut )
{
    testBasicVector( ut );

#ifdef USE_EXT_PETSC
    {
        auto simplePetscFactory = std::make_shared<SimplePetscVectorFactory>( d_factory );
        auto petscViewFactory   = std::make_shared<PetscViewFactory>( simplePetscFactory );
        auto petscCloneFactory  = std::make_shared<PetscCloneFactory>( petscViewFactory );
        PetscVectorTests test1( petscViewFactory );
        PetscVectorTests test2( petscCloneFactory );
        test1.testPetscVector( ut );
        test2.testPetscVector( ut );
    }
#endif

#ifdef USE_EXT_SUNDIALS
    {
        auto viewFactory =
            std::make_shared<ViewFactory<AMP::LinearAlgebra::SundialsVector>>( d_factory );
        auto cloneFactory = std::make_shared<CloneFactory>( viewFactory );
        VectorTests test1( viewFactory );
        VectorTests test2( cloneFactory );
        SundialsVectorTests test3( viewFactory );
        SundialsVectorTests test4( cloneFactory );
        test1.testBasicVector( ut );
        test2.testBasicVector( ut );
        test3.testSundialsVector( ut );
        test4.testSundialsVector( ut );
    }
#endif
}


void VectorTests::testNullVector( AMP::UnitTest *ut )
{
    auto viewFactory = std::make_shared<NullVectorFactory>();
    VectorTests test( viewFactory );
    test.InstantiateVector( ut );
}


void VectorTests::testParallelVectors( AMP::UnitTest *ut )
{
    InstantiateVector( ut );
    // VerifyVectorGhostCreate( ut );
    VerifyVectorMakeConsistentSet( ut );
    VerifyVectorMakeConsistentAdd( ut );
    CopyVectorConsistency( ut );
}


void VectorTests::testVectorSelector( AMP::UnitTest *ut )
{
    testAllSelectors( ut );
    test_VS_ByVariableName( ut );
    test_VS_Comm( ut );
}


} // namespace LinearAlgebra
} // namespace AMP
