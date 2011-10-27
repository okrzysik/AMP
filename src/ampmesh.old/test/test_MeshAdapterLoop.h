#include "string.h"
#include "test_MeshTests.h"
#include "utils/UnitTest.h"

using namespace AMP::unit_test;

template <typename GENERATOR>
class MeshAdapterLoop
{
public:
    static void  test_mesh( AMP::UnitTest *ut )
    {
        // Create the mesh
        GENERATOR generator;
        AMP::Mesh::MeshAdapter::shared_ptr mesh = generator.getMesh();
        // Test iterators
        VerifyBoundaryNodeIterator::run_test( ut, mesh );
        VerifyOwnedBoundaryNodeIterator::run_test( ut, mesh );
        VerifyOwnedNodeIterator::run_test( ut, mesh );
        ElementIteratorTest::run_test( ut, mesh );
        NodeIteratorTest::run_test( ut, mesh );
        VerifyNodeElemMapIteratorTest::run_test( ut, mesh );
        ElementTest::run_test( ut, mesh );
        // Test elements
        VerifyBoundaryIteratorTest::run_test( ut, mesh );
        // Test Interface
        VerifyProcAndIsOwnedInterface<ElementHelper>::run_test( ut, mesh );
        VerifyProcAndIsOwnedInterface<NodeHelper>::run_test( ut, mesh );
        // Test element for node
        VerifyElementForNode::run_test( ut, mesh );
        // Test displacement
        DisplaceNodes::run_test( ut, mesh );
        // Bug tests
        Bug_758::run_test( ut, mesh );
        Bug_761<1>::run_test( ut, mesh );
        Bug_761<2>::run_test( ut, mesh );
        Bug_761<7>::run_test( ut, mesh );
        Bug_761<8>::run_test( ut, mesh );
        //MeshAdapterTest<AllPassTest>::run_test( ut, mesh );
    }
};


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

