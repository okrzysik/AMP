#ifndef included_test_VectorTests
#define included_test_VectorTests

#include "utils/UnitTest.h"


namespace AMP {
namespace LinearAlgebra {


/**
 * \class vectorTests
 * \brief A helper class to store/run tests for a vector
 */
template <typename VECTOR_FACTORY>
class vectorTests
{
public:

    static void testPetscVector( AMP::UnitTest *ut );

    static void testBasicVector( AMP::UnitTest *ut );

    static void testSundialsVector( AMP::UnitTest *ut );

    static void testManagedVector( AMP::UnitTest *ut );

    static void testNullVector( AMP::UnitTest *ut );

    static void testParallelVectors( AMP::UnitTest *ut );

    static void testVectorSelector( AMP::UnitTest *ut );


public:

    static void InstantiateVector( AMP::UnitTest *utils );


    static void CopyVectorConsistency( AMP::UnitTest *utils );


    template <typename VIEWER>
    static void DeepCloneOfView( AMP::UnitTest *utils );


    static void Bug_728( AMP::UnitTest *utils );


    static void SetToScalarVector( AMP::UnitTest *utils );


    static void CloneVector( AMP::UnitTest *utils );


    static void DotProductVector( AMP::UnitTest *utils );


    static void L2NormVector( AMP::UnitTest *utils );


    static void AbsVector( AMP::UnitTest *utils );


    static void L1NormVector( AMP::UnitTest *utils );


    static void MaxNormVector( AMP::UnitTest *utils );


    static void ScaleVector( AMP::UnitTest *utils );


#ifdef USE_EXT_PETSC
    static void Bug_491( AMP::UnitTest *utils );
#endif


    static void AddVector( AMP::UnitTest *utils );


    static void SubtractVector( AMP::UnitTest *utils );


    static void MultiplyVector( AMP::UnitTest *utils );


    static void DivideVector( AMP::UnitTest *utils );


    static void VectorIteratorLengthTest( AMP::UnitTest *utils );


    template <typename ITERATOR>
    static void both_VectorIteratorTests( AMP::LinearAlgebra::Vector::shared_ptr p, AMP::UnitTest *utils );


    static void VectorIteratorTests( AMP::UnitTest *utils );


    static void VerifyVectorMin( AMP::UnitTest *utils );


    static void VerifyVectorMax( AMP::UnitTest *utils );


    static void VerifyVectorMaxMin( AMP::UnitTest *utils );


    class SetRandomValuesVector {
      public:
        static inline const char *get_test_name() { return "vector::setRandomValues"; }
        static void verify_vector( AMP::UnitTest *utils, AMP::LinearAlgebra::Vector::shared_ptr v );
        static void run_test( AMP::UnitTest *utils );
    };


    static void ReciprocalVector( AMP::UnitTest *utils );


    class LinearSumVector {
      public:
        static inline const char *get_test_name() { return "vector::linearSum"; }
        static void do_instance( AMP::UnitTest *utils, double alpha, double beta, const char *msg );
        static void run_test( AMP::UnitTest *utils );
    };


    class AxpyVector {
      public:
        static inline const char *get_test_name() { return "vector::axpy"; }
        static void do_instance( AMP::UnitTest *utils, double alpha, const char *msg );
        static void run_test( AMP::UnitTest *utils );
    };


    class AxpbyVector {
      public:
        static inline const char *get_test_name() { return "vector::axpby"; }
        static void do_instance( AMP::UnitTest *utils, double alpha, double beta, const char *msg );
        static void run_test( AMP::UnitTest *utils );
    };


    class CopyVector {
      public:
        static inline const char *get_test_name() { return "vector::copyVector"; }
        static void run_test( AMP::UnitTest *utils );
    };

    class CopyRawDataBlockVector {
      public:
        static inline const char *get_test_name() { return "vector::copyVector"; }
        static void run_test( AMP::UnitTest *utils );
    };


    static void VerifyVectorGhostCreate( AMP::UnitTest *utils );


    static void VerifyVectorMakeConsistentAdd( AMP::UnitTest *utils );


    static void VerifyVectorMakeConsistentSet( AMP::UnitTest *utils );


    // Test creating a multivector with multiple copies of the data
    // This should always return one copy of the superset of the data
    static void TestMultivectorDuplicate( AMP::UnitTest *utils );

};


} // namespace LinearAlgebra
} // namespace AMP


// Extra includes
#include "vectors/testHelpers/VectorTests.inline.h"
#include "vectors/testHelpers/vectorTestLoop.inline.h"


#endif
