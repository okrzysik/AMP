#include "string.h"
#include "meshTests.h"
#include "utils/UnitTest.h"
#include "ampmesh/Mesh.h"
#include "meshTests.h"


void MeshTestLoop( AMP::UnitTest *ut, boost::shared_ptr<AMP::Mesh::Mesh> mesh )
{
    // Test the iterators
    MeshIteratorTest( ut, mesh );
    //VerifyBoundaryNodeIterator::run_test( ut, mesh );
    //VerifyOwnedBoundaryNodeIterator::run_test( ut, mesh );
    //VerifyOwnedNodeIterator::run_test( ut, mesh );
    //VerifyNodeElemMapIteratorTest::run_test( ut, mesh );
    //ElementTest::run_test( ut, mesh );
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
};


#ifdef USE_VECTORS
template <typename GENERATOR>
class MeshAdapterVectorLoop
{
public:
    static void  test_mesh( AMP::UnitTest *ut )
    {
        // Create the mesh
        GENERATOR generator;
        AMP::Mesh::MeshAdapter::shared_ptr mesh = generator.getMesh();
        // Run the vector tests
        VerifyGetVectorTest<1>::run_test ( ut, mesh );
        VerifyGetVectorTest<3>::run_test ( ut, mesh );
    }
};
#endif


#ifdef USE_MATRICIES
template <typename GENERATOR>
class MeshAdapterMatrixLoop
{
public:
    static void  test_mesh( AMP::UnitTest *ut )
    {
        // Create the mesh
        GENERATOR generator;
        AMP::Mesh::MeshAdapter::shared_ptr mesh = generator.getMesh();
        // Run the matrix tests
        VerifyGetMatrixTrivialTest<1>::run_test ( ut, mesh );
        VerifyGetMatrixTrivialTest<3>::run_test ( ut, mesh );
    }
};
#endif


