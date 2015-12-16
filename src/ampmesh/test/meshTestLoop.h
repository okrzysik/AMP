#include "ampmesh/Mesh.h"
#include "ampmesh/SubsetMesh.h"
#include "meshTests.h"
#include "meshTests.h"
#include "string.h"
#include "utils/UnitTest.h"
#ifdef USE_AMP_VECTORS
#include "meshVectorTests.h"
#endif
#ifdef USE_AMP_MATRICES
#include "meshMatrixTests.h"
#endif
#include "ProfilerApp.h"


void MeshTestLoop( AMP::UnitTest *ut, AMP::shared_ptr<AMP::Mesh::Mesh> mesh )
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
    if ( AMP::dynamic_pointer_cast<AMP::Mesh::SubsetMesh>( mesh ).get() != nullptr )
        ut->expected_failure( "Displace mesh tests are not valid for sub-meshes" );
    else
        DisplaceMesh( ut, mesh );

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
    PROFILE_STOP( "MeshTestLoop" );
}


void MeshVectorTestLoop( AMP::UnitTest *ut, AMP::shared_ptr<AMP::Mesh::Mesh> mesh )
{
// Run the vector tests
#ifdef USE_AMP_VECTORS
    PROFILE_START( "MeshVectorTestLoop" );
    VerifyGetVectorTest<1, false>( ut, mesh );
    VerifyGetVectorTest<3, false>( ut, mesh );
    VerifyGetVectorTest<1, true>( ut, mesh );
    VerifyGetVectorTest<3, true>( ut, mesh );
    PROFILE_STOP( "MeshVectorTestLoop" );
#endif
}


void MeshMatrixTestLoop( AMP::UnitTest *ut, AMP::shared_ptr<AMP::Mesh::Mesh> mesh )
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
        VerifyGetMatrixTrivialTest<1, false>( ut, mesh );
        VerifyGetMatrixTrivialTest<3, false>( ut, mesh );
        VerifyGetMatrixTrivialTest<1, true>( ut, mesh );
        VerifyGetMatrixTrivialTest<3, true>( ut, mesh );
        GhostWriteTest<1, false>( ut, mesh );
        GhostWriteTest<3, false>( ut, mesh );
        GhostWriteTest<1, true>( ut, mesh );
        GhostWriteTest<3, true>( ut, mesh );
        PROFILE_STOP( "MeshMatrixTestLoop" );
    }
#endif
}
