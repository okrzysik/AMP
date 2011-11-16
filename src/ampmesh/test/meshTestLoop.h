#include "string.h"
#include "meshTests.h"
#include "utils/UnitTest.h"
#include "ampmesh/Mesh.h"
#include "meshTests.h"
#ifdef USE_AMP_VECTORS
    #include "meshVectorTests.h"
#endif
#ifdef USE_AMP_MATRICIES
    #include "meshMatrixTests.h"
#endif



void MeshTestLoop( AMP::UnitTest *ut, boost::shared_ptr<AMP::Mesh::Mesh> mesh )
{
    // Run some basic sanity checks
    MeshBasicTest( ut, mesh );
    // Test the number of elements
    MeshCountTest( ut, mesh );
    // Test the iterators
    MeshIteratorTest( ut, mesh );
    //VerifyBoundaryNodeIterator::run_test( ut, mesh );
    //VerifyOwnedBoundaryNodeIterator::run_test( ut, mesh );
    //VerifyOwnedNodeIterator::run_test( ut, mesh );
    //VerifyNodeElemMapIteratorTest::run_test( ut, mesh );
    // Test the elements
    //VerifyBoundaryIteratorTest::run_test( ut, mesh );
    // Test Interface
    //VerifyProcAndIsOwnedInterface<ElementHelper>::run_test( ut, mesh );
    //VerifyProcAndIsOwnedInterface<NodeHelper>::run_test( ut, mesh );
    // Test element for node
    //VerifyElementForNode::run_test( ut, mesh );
    // Test displacement
    //DisplaceNodes::run_test( ut, mesh );
    // Bug tests
    //Bug_758::run_test( ut, mesh );
    //Bug_761<1>::run_test( ut, mesh );
    //Bug_761<2>::run_test( ut, mesh );
    //Bug_761<7>::run_test( ut, mesh );
    //Bug_761<8>::run_test( ut, mesh );
    //MeshAdapterTest<AllPassTest>::run_test( ut, mesh );


    // Run the vector tests
    #ifdef USE_AMP_VECTORS
        VerifyGetVectorTest<1>( ut, mesh );
        VerifyGetVectorTest<3>( ut, mesh );
    #endif

    // Run the matrix tests
    #ifdef USE_AMP_MATRICIES
        ut->failure("Matricies are not implimented yet");
        //VerifyGetMatrixTrivialTest<1>( ut, mesh );
        //VerifyGetMatrixTrivialTest<3>( ut, mesh );
    #endif

};

