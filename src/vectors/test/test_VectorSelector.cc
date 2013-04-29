#include "vectors/VectorSelector.h"
#include "vectors/MultiVector.h"
#include "vectors/VectorSelector.h"

#include "test_VectorLoops.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"

#ifdef USE_EXT_TRILINOS
    #include "vectors/trilinos/ManagedEpetraVector.h"
#endif

using namespace AMP::unit_test;


template <typename VECTOR_FACTORY>
class SelectTester
{
public:
    static const char * get_test_name () { return "Vector selection mechanism"; }

    static  void run_test ( AMP::UnitTest *utils )
    {
        AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
        AMP::LinearAlgebra::Vector::shared_ptr  vec1 = VECTOR_FACTORY::getVector();
        AMP::LinearAlgebra::Vector::shared_ptr  vec2 = vec1->cloneVector ( "vec2" );
        AMP::LinearAlgebra::Vector::shared_ptr  vec3a = vec1->cloneVector ( "vec3" );
        AMP::LinearAlgebra::Vector::shared_ptr  vec3b = vec1->cloneVector ( "vec3" );
        AMP::LinearAlgebra::Vector::shared_ptr  vec3 = AMP::LinearAlgebra::MultiVector::create ( "vec3" , globalComm );
        vec3->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( vec2 );
        vec3->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( vec3a );
        vec3->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( vec3b );

        // utils.failure ( "This is a test failure" );
        AMP::LinearAlgebra::Vector::shared_ptr  selection1 = vec2->select ( AMP::LinearAlgebra::VS_ByVariableName ( "None" ) , "None" );
        if ( selection1 )
            utils->failure ( "Found vector where there should be none" );
        else
            utils->passes ( "Excluded all vectors" );

        selection1 = vec2->select ( AMP::LinearAlgebra::VS_ByVariableName ( "vec2" ) , "subset" );
        if ( selection1 ) {
            if ( selection1->castTo<AMP::LinearAlgebra::MultiVector>().getVector(0).get() != vec2.get() )
                utils->failure ( "Could not find vector" );
            else
                utils->passes ( "Found vector" );
        } else {
            utils->failure ( "Did not find a vector" );
            return;
        }

        selection1 = vec3->select ( AMP::LinearAlgebra::VS_ByVariableName ( "vec3" ) , "subset" );
        if ( selection1 ) {
            if ( selection1->castTo<AMP::LinearAlgebra::MultiVector>().getVector(0).get() != vec3a.get() )
                utils->failure ( "Could not find vector" );
            else
                utils->passes ( "Found vector" );
            if ( selection1->castTo<AMP::LinearAlgebra::MultiVector>().getVector(1).get() != vec3b.get() )
                utils->failure ( "Could not find vector" );
            else
                utils->passes ( "Found vector" );
        } else {
            utils->failure ( "Did not find a vector" );
            return;
        }

    }

};


template <typename T>
class StridedVectorFactory
{
public:

    static AMP::LinearAlgebra::Variable::shared_ptr  getVariable()
    {
        return AMP::LinearAlgebra::Variable::shared_ptr ();
    }

    static AMP::LinearAlgebra::Vector::shared_ptr   getVector()
    {
        boost::shared_ptr<typename T::vector>  vec = T::getVector();
        AMP::LinearAlgebra::VS_Stride criterion = AMP::LinearAlgebra::VS_Stride(1,3);
        AMP::LinearAlgebra::Vector::shared_ptr  vec_select = vec->select( criterion, "thirds" );
        size_t N1 = vec->getGlobalSize();
        size_t N2 = vec_select->getGlobalSize();
        AMP_ASSERT(N1/3==N2);
        return vec_select;
    }
};


#ifdef USE_EXT_TRILINOS
typedef SimpleManagedVectorFactory<AMP::LinearAlgebra::ManagedEpetraVector>         SMEVFactory;
#endif


int main ( int argc , char **argv )
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    #ifdef USE_EXT_TRILINOS
        SelectTester<SMEVFactory>::run_test ( &ut );
        testManagedVector<StridedVectorFactory<SMEVFactory> > ( &ut );
    #else
        ut.expected_failure("Compiled without trilinos");
    #endif

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}


