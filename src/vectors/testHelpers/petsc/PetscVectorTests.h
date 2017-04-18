#ifdef USE_EXT_PETSC
#ifndef included_test_PetscVectorTests
#define included_test_PetscVectorTests

#include "string"
#include <algorithm>

#include "utils/UnitTest.h"
#include "vectors/petsc/PetscVector.h"
#include "vectors/petsc/NativePetscVector.h"
#include "vectors/testHelpers/VectorTests.h"

namespace AMP {
namespace LinearAlgebra {

class PetscVectorFactory;

/**
 * \class PetscVectorTests
 * \brief A helper class to store/run tests for a vector
 */
class PetscVectorTests
{
public:
    explicit PetscVectorTests( AMP::shared_ptr<const PetscVectorFactory> factory ): d_factory(factory) {}

    void testPetscVector( AMP::UnitTest *ut );

    void InstantiatePetscVectors( AMP::UnitTest *utils );

    void Bug_612( AMP::UnitTest *utils );

    void DuplicatePetscVector( AMP::UnitTest *utils );

    void StaticDuplicatePetscVector( AMP::UnitTest *utils );

    void StaticCopyPetscVector( AMP::UnitTest *utils );

    void CopyPetscVector( AMP::UnitTest *utils );

    void VerifyPointwiseMaxAbsPetscVector( AMP::UnitTest *utils );

    void VerifyPointwiseMaxPetscVector( AMP::UnitTest *utils );

    void VerifyPointwiseMinPetscVector( AMP::UnitTest *utils );

    void VerifyAXPBYPCZPetscVector( AMP::UnitTest *utils );

    void VerifyAYPXPetscVector( AMP::UnitTest *utils );

    void VerifyExpPetscVector( AMP::UnitTest *utils );

    void VerifyLogPetscVector( AMP::UnitTest *utils );

    void VerifyNormsPetscVector( AMP::UnitTest *utils );

    void VerifyAXPBYPetscVector( AMP::UnitTest *utils );

    void VerifySwapPetscVector( AMP::UnitTest *utils );

    void VerifyGetSizePetscVector( AMP::UnitTest *utils );

    void VerifyMaxPointwiseDividePetscVector( AMP::UnitTest *utils );

    void VerifyAbsPetscVector( AMP::UnitTest *utils );

    void VerifyPointwiseMultPetscVector( AMP::UnitTest *utils );

    void VerifyPointwiseDividePetscVector( AMP::UnitTest *utils );

    void VerifySqrtPetscVector( AMP::UnitTest *utils );

    void VerifySetRandomPetscVector( AMP::UnitTest *utils );

    void VerifySetPetscVector( AMP::UnitTest *utils );

    void VerifyAXPYPetscVector( AMP::UnitTest *utils );

    void VerifyScalePetscVector( AMP::UnitTest *utils );

    void VerifyDotPetscVector( AMP::UnitTest *utils );

private:
    AMP::shared_ptr<const PetscVectorFactory> d_factory;
};


} // LinearAlgebra namespace
} // AMP namespace

/// \endcond

#endif
#endif
