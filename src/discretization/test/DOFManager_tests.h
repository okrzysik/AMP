#ifndef included_AMP_DOFManager_tests
#define included_AMP_DOFManager_tests

#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/MultiDOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/discretization/subsetDOFManager.h"
#include "AMP/utils/UnitTest.h"

#ifdef USE_AMP_VECTORS
#include "AMP/vectors/MeshVariable.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/VectorSelector.h"
#endif

#include "../../ampmesh/test/meshGenerators.h"


// Function to test getting the DOFs for a mesh iterator
void testGetDOFIterator( AMP::UnitTest *ut,
                         const AMP::Mesh::MeshIterator &iterator,
                         AMP::Discretization::DOFManager::shared_ptr DOF )
{
    bool pass1 = true;
    bool pass2 = true;
    auto cur   = iterator.begin();
    std::vector<size_t> dofs;
    for ( size_t i = 0; i < cur.size(); i++ ) {
        const auto id = cur->globalID();
        DOF->getDOFs( id, dofs );
        if ( dofs.empty() )
            pass1 = false;
        for ( size_t dof : dofs ) {
            auto id2 = DOF->getElement( dof ).globalID();
            if ( id2 != id )
                pass2 = false;
        }
        ++cur;
    }
    if ( pass1 )
        ut->passes( "Got the DOFs for every element in iterator" );
    else
        ut->failure( "Got the DOFs for every element in iterator" );
    if ( pass2 )
        ut->passes( "getElementID" );
    else
        ut->failure( "getElementID" );
}


// Function to test some very basic properites
void testBasics( AMP::Discretization::DOFManager::shared_ptr DOF, AMP::UnitTest *ut )
{
    bool passAll = true;

    // Get the comm
    AMP::AMP_MPI comm = DOF->getComm();
    if ( comm.isNull() ) {
        passAll = false;
        ut->failure( "Comm is not valid" );
    }

    // Get the dof ranges
    size_t start    = DOF->beginDOF();
    size_t end      = DOF->endDOF();
    size_t N_local  = DOF->numLocalDOF();
    size_t N_global = DOF->numGlobalDOF();
    if ( N_local != end - start ) {
        passAll = false;
        ut->failure( "N_local is the wrong size" );
    }
    if ( N_global != comm.sumReduce( N_local ) ) {
        passAll = false;
        ut->failure( "N_global is the wrong size" );
    }

    // Trivial check of == and !=
    if ( *DOF != *DOF || !( *DOF == *DOF ) ) {
        passAll = false;
        ut->failure( "Failed operator==" );
    }

    // Test getting the DOFs for the DOF iterator
    AMP::UnitTest ut2;
    testGetDOFIterator( &ut2, DOF->getIterator(), DOF );
    if ( ut2.NumFailLocal() != 0 ) {
        passAll = false;
        ut->failure( "Failed checking local iterator" );
    }

    // Check the results of the basic tests
    if ( passAll )
        ut->passes( "Basic DOF tests" );
    else
        ut->failure( "Basic DOF tests" );
}


// Test subsetting for different comms
void testSubsetComm( AMP::Discretization::DOFManager::shared_ptr DOF, AMP::UnitTest *ut )
{
    auto subsetDOF = DOF->subset( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    if ( *DOF == *subsetDOF )
        ut->passes( "Subset DOF for COMM_WORLD" );
    else
        ut->failure( "Subset DOF on COMM_WORLD" );
    subsetDOF = DOF->subset( DOF->getComm() );
    if ( *DOF == *subsetDOF )
        ut->passes( "Subset DOF for DOF comm" );
    else
        ut->failure( "Subset DOF on DOF comm" );
    subsetDOF = DOF->subset( AMP::AMP_MPI( AMP_COMM_SELF ) );
    if ( subsetDOF->numLocalDOF() == DOF->numLocalDOF() &&
         subsetDOF->numGlobalDOF() == DOF->numLocalDOF() )
        ut->passes( "Subset DOF for COMM_SELF" );
    else
        ut->failure( "Subset DOF on COMM_SELF" );
}


// Test subsetting for different meshes
void testSubsetMesh( AMP::Mesh::Mesh::shared_ptr mesh,
                     AMP::Discretization::DOFManager::shared_ptr DOF,
                     bool is_nodal,
                     int DOFsPerNode,
                     int gcw,
                     AMP::UnitTest *ut )
{
    auto subsetDOF = DOF->subset( mesh );
    if ( *DOF == *subsetDOF )
        ut->passes( "Subset DOF on full mesh" );
    else
        ut->failure( "Subset DOF on full mesh" );
    auto meshIDs = mesh->getBaseMeshIDs();
    if ( meshIDs.size() > 1 && is_nodal ) {
        bool passes = true;
        for ( size_t i = 0; i < meshIDs.size(); i++ ) {
            auto subsetMesh = mesh->Subset( meshIDs[i] );
            if ( subsetMesh.get() != nullptr ) {
                subsetDOF = DOF->subset( subsetMesh );
                testGetDOFIterator(
                    ut, subsetMesh->getIterator( AMP::Mesh::GeomType::Vertex, gcw ), subsetDOF );
                auto mesh_DOF = AMP::Discretization::simpleDOFManager::create(
                    subsetMesh, AMP::Mesh::GeomType::Vertex, gcw, DOFsPerNode, false );
                if ( *mesh_DOF != *subsetDOF )
                    passes = false;
            }
        }
        if ( passes )
            ut->passes( "Subset DOF for each mesh" );
        else
            ut->failure( "Subset DOF for each mesh" );
    }
}


// Function to test the index conversion given a multDOFManager
void testMultiDOFMap( AMP::UnitTest *ut,
                      std::shared_ptr<AMP::Discretization::multiDOFManager> multiDOF )
{
    // First create a global DOF list
    size_t N_global = multiDOF->numGlobalDOF();
    std::vector<size_t> globalDOFList( N_global );
    for ( size_t i = 0; i < N_global; i++ )
        globalDOFList[i] = i;

    // Loop through the DOFManagers
    auto managers = multiDOF->getDOFManagers();
    for ( size_t i = 0; i < managers.size(); i++ ) {
        AMP::AMP_MPI comm = managers[i]->getComm();
        comm.barrier();
        // Convert the global list to a local list
        std::vector<size_t> tmp = multiDOF->getSubDOF( (int) i, globalDOFList );
        std::vector<size_t> localDOFList, localToGlobal;
        localDOFList.reserve( N_global );
        localToGlobal.reserve( N_global );
        size_t neg_one = ~( (size_t) 0 );
        for ( size_t j = 0; j < N_global; j++ ) {
            if ( tmp[j] != neg_one ) {
                localDOFList.push_back( tmp[j] );
                localToGlobal.push_back( j );
            }
        }
        // Check that we created the proper list
        bool passes = localDOFList.size() == managers[i]->numGlobalDOF();
        for ( size_t j = 0; j < localDOFList.size(); j++ ) {
            if ( localDOFList[j] != j )
                passes = false;
        }
        if ( passes )
            ut->passes( "Conversion from global to sub DOFs in multiDOFManager" );
        else
            ut->failure( "Conversion from global to sub DOFs in multiDOFManager" );
        // Check that we can convert back to the global ids
        auto globalDOFList2 = multiDOF->getGlobalDOF( (int) i, localDOFList );
        passes              = globalDOFList2.size() == localDOFList.size();
        for ( size_t j = 0; j < globalDOFList2.size(); j++ ) {
            if ( globalDOFList2[j] != localToGlobal[j] )
                passes = false;
        }
        if ( passes )
            ut->passes( "Conversion from sub to global DOFs in multiDOFManager" );
        else
            ut->failure( "Conversion from sub to global DOFs in multiDOFManager" );
    }
}


// Function to test that a multivector with a DOFManager repeated correctly sets the values
#ifdef USE_AMP_VECTORS
void testMultiDOFVector( AMP::UnitTest *ut, AMP::Discretization::DOFManager::shared_ptr DOF )
{
    // Create the individual vectors
    auto var1 = std::make_shared<AMP::LinearAlgebra::Variable>( "a" );
    auto var2 = std::make_shared<AMP::LinearAlgebra::Variable>( "b" );
    auto var3 = std::make_shared<AMP::LinearAlgebra::Variable>( "c" );
    auto vec1 = AMP::LinearAlgebra::createVector( DOF, var1, true );
    auto vec2 = AMP::LinearAlgebra::createVector( DOF, var2, true );
    // Create the multivector
    auto multiVector = AMP::LinearAlgebra::MultiVector::create( var3, DOF->getComm() );
    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> vectors( 2 );
    vectors[0] = vec1;
    vectors[1] = vec2;
    multiVector->addVector( vectors );
    auto multiDOF = multiVector->getDOFManager();
    // Check that we can set each value correctly
    multiVector->zero(*multiVector);
    multiVector->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
    AMP::Mesh::MeshIterator it = DOF->getIterator();
    std::vector<size_t> dof1, dof2;
    bool uniqueMultiDOFs = true;
    for ( size_t i = 0; i < it.size(); i++ ) {
        DOF->getDOFs( it->globalID(), dof1 );
        multiDOF->getDOFs( it->globalID(), dof2 );
        AMP_ASSERT( dof1.size() == 1 );
        AMP_ASSERT( dof2.size() == 2 );
        if ( dof2[0] == dof2[1] )
            uniqueMultiDOFs = false;
        multiVector->setValueByGlobalID( dof2[0], dof1[0] );
        multiVector->setValueByGlobalID( dof2[1], 2 * dof1[0] );
        ++it;
    }
    if ( uniqueMultiDOFs )
        ut->passes( "MultiDOFManger with duplicate subDOFManagers returns unique DOFs" );
    else
        ut->failure( "MultiDOFManger with duplicate subDOFManagers returns unique DOFs" );
    multiVector->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
    double vec1norm        = vec1->L1Norm(vec1);
    double vec2norm        = vec2->L1Norm(vec2);
    double multiVectorNorm = multiVector->L1Norm(*multiVector);
    double N_tot           = ( DOF->numGlobalDOF() * ( DOF->numGlobalDOF() - 1 ) ) / 2;
    if ( vec1norm == N_tot && vec2norm == 2 * N_tot && multiVectorNorm == 3 * N_tot )
        ut->passes( "MultiVector with repeated DOFs sets values correctly" );
    else
        ut->failure( "MultiVector with repeated DOFs sets values correctly" );
}
#endif


#endif
