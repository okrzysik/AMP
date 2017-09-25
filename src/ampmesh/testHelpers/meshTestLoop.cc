#include "ampmesh/testHelpers/meshTests.h"

#include "ProfilerApp.h"


namespace AMP {
namespace Mesh {


void meshTests::MeshTestLoop( AMP::UnitTest *ut, AMP::shared_ptr<AMP::Mesh::Mesh> mesh )
{
    PROFILE_START( "MeshTestLoop" );
    // Run some basic sanity checks
    MeshBasicTest( ut, mesh );
    // Test the number of elements
    MeshCountTest( ut, mesh );
    // Test the iterators
    MeshIteratorTest( ut, mesh );
    MeshIteratorOperationTest( ut, mesh );
    VerifyBoundaryIDNodeIterator( ut, mesh );
    VerifyBoundaryIterator( ut, mesh );
    testBlockIDs( ut, mesh );
    MeshIteratorSetOPTest( ut, mesh );
    VerifyGhostIsOwned( ut, mesh );
    // Test the node neighbors
    getNodeNeighbors( ut, mesh );
    // Test displacement
    if ( mesh->isMeshMovable() >= 1 )
        DisplaceMeshScalar( ut, mesh );
    if ( mesh->isMeshMovable() >= 2 )
        DisplaceMeshVector( ut, mesh );
    // VerifyNodeElemMapIteratorTest::run_test( ut, mesh );
    // Test the elements
    // VerifyBoundaryIteratorTest::run_test( ut, mesh );
    // Test Interface
    // VerifyProcAndIsOwnedInterface<ElementHelper>::run_test( ut, mesh );
    // VerifyProcAndIsOwnedInterface<NodeHelper>::run_test( ut, mesh );
    // Test element for node
    // VerifyElementForNode::run_test( ut, mesh );
    // Bug tests
    // Bug_758::run_test( ut, mesh );
    // Bug_761<1>::run_test( ut, mesh );
    // Bug_761<2>::run_test( ut, mesh );
    // Bug_761<7>::run_test( ut, mesh );
    // Bug_761<8>::run_test( ut, mesh );
    // MeshAdapterTest<AllPassTest>::run_test( ut, mesh );
    MeshPerformance( ut, mesh );
    ;
    PROFILE_STOP( "MeshTestLoop" );
}


void meshTests::MeshVectorTestLoop(
    AMP::UnitTest *ut, AMP::shared_ptr<AMP::Mesh::Mesh> mesh, bool fast )
{
// Run the vector tests
#ifdef USE_AMP_VECTORS
    PROFILE_START( "MeshVectorTestLoop" );
    VerifyGetVectorTest<1, true>( ut, mesh );
    if ( !fast ) {
        VerifyGetVectorTest<3, true>( ut, mesh );
        VerifyGetVectorTest<1, false>( ut, mesh );
        VerifyGetVectorTest<3, false>( ut, mesh );
    }
    PROFILE_STOP( "MeshVectorTestLoop" );
#endif
}


void meshTests::MeshMatrixTestLoop(
    AMP::UnitTest *ut, AMP::shared_ptr<AMP::Mesh::Mesh> mesh, bool fast )
{
// Run the matrix tests
#ifdef USE_AMP_MATRICES
    bool run_tests = true;
#if !defined( USE_EXT_PETSC ) && !defined( USE_EXT_TRILINOS )
    if ( AMP::AMP_MPI( AMP_COMM_WORLD ).getSize() > 1 ) {
        ut->expected_failure( "No parallel matrix to test" );
        run_tests = false;
    }
#endif
    if ( run_tests ) {
        PROFILE_START( "MeshMatrixTestLoop" );
        VerifyGetMatrixTrivialTest<1, true>( ut, mesh );
        GhostWriteTest<1, true>( ut, mesh );
        if ( !fast ) {
            VerifyGetMatrixTrivialTest<3, true>( ut, mesh );
            VerifyGetMatrixTrivialTest<1, false>( ut, mesh );
            VerifyGetMatrixTrivialTest<3, false>( ut, mesh );
            GhostWriteTest<3, true>( ut, mesh );
            GhostWriteTest<1, false>( ut, mesh );
            GhostWriteTest<3, false>( ut, mesh );
        }
        PROFILE_STOP( "MeshMatrixTestLoop" );
    }
#endif
}


} // namespace Mesh
} // namespace AMP
