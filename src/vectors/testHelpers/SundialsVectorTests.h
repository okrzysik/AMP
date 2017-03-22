#ifdef USE_EXT_SUNDIALS
#ifndef included_test_SundialsVectorTests
#define included_test_SundialsVectorTests

#include "string"
#include <algorithm>

#include "utils/UnitTest.h"
#include "vectors/petsc/PetscVector.h"
#include "vectors/petsc/NativePetscVector.h"
#include "vectors/testHelpers/VectorTests.h"
#include "vectors/testHelpers/PetscVectorTests.h"
#include "vectors/sundials/ManagedSundialsVector.h"


namespace AMP {
namespace LinearAlgebra {



/**
 * \class SundialsVectorTests
 * \brief A helper class to store/run tests for a vector
 */
class SundialsVectorTests
{
public:
    SundialsVectorTests( AMP::shared_ptr<const VectorFactory> factory ): d_factory(factory) {}

    void testSundialsVector( AMP::UnitTest *ut );

    void CloneSundialsVector( AMP::UnitTest *utils );

    void LinearSumSundialsVector( AMP::UnitTest *utils );

    void ConstSundialsVector( AMP::UnitTest *utils );

    void ProdSundialsVector( AMP::UnitTest *utils );

    void DivSundialsVector( AMP::UnitTest *utils );

    void ScaleSundialsVector( AMP::UnitTest *utils );

    void AbsSundialsVector( AMP::UnitTest *utils );

    void InvSundialsVector( AMP::UnitTest *utils );

    void AddConstSundialsVector( AMP::UnitTest *utils );

    void DotProdSundialsVector( AMP::UnitTest *utils );

    void MaxNormSundialsVector( AMP::UnitTest *utils );

    void WRMSNormSundialsVector( AMP::UnitTest *utils );

    void MinSundialsVector( AMP::UnitTest *utils );

    void L1NormSundialsVector( AMP::UnitTest *utils );

    void MinQuotientSundialsVector( AMP::UnitTest *utils );


private:
    AMP::shared_ptr<const VectorFactory> d_factory;

    static inline AMP::LinearAlgebra::ManagedSundialsVector *getVector( N_Vector &a )
    {
        return reinterpret_cast<AMP::LinearAlgebra::ManagedSundialsVector *>( a->content );
    }

};


} // LinearAlgebra namespace
} // AMP namespace

/// \endcond

#endif
#endif
