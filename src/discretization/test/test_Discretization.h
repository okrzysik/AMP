#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/MultiDOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/discretization/structuredFaceDOFManager.h"
#include "AMP/discretization/subsetDOFManager.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MultiMesh.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Utilities.h"

#include "DOFManager_tests.h"

#include "../../mesh/test/meshGenerators.h"


using namespace AMP::unit_test;


// Function to test subsetting a DOF manager
template<class GENERATOR, bool SPLIT>
void testSubsetDOFManager( AMP::UnitTest *ut )
{
    // Get the mesh
    GENERATOR mesh_generator;
    mesh_generator.build_mesh();
    auto mesh = mesh_generator.getMesh();

    // Create a simple DOF manager
    auto DOF = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 1, SPLIT );
    testGetDOFIterator( ut, mesh->getIterator( AMP::Mesh::GeomType::Vertex, 1 ), DOF );

    // Run basic tests that should be valid for all DOFManagers
    testBasics( DOF, ut );
    testSubsetComm( DOF, ut );
    testSubsetMesh( mesh, DOF, true, 1, 1, ut );

    // Subset for iterators
    auto iterator  = mesh->getIterator( AMP::Mesh::GeomType::Vertex, 1 );
    auto subsetDOF = DOF->subset( iterator, mesh->getComm() );
    testGetDOFIterator( ut, mesh->getIterator( AMP::Mesh::GeomType::Vertex, 1 ), subsetDOF );
    if ( *DOF == *subsetDOF )
        ut->passes( "Subset DOF on full mesh iterator: " + GENERATOR::name() );
    else
        ut->failure( "Subset DOF on full mesh iterator: " + GENERATOR::name() );
    iterator  = mesh->getSurfaceIterator( AMP::Mesh::GeomType::Vertex, 1 );
    subsetDOF = DOF->subset( iterator, mesh->getComm() );
    testGetDOFIterator( ut, mesh->getSurfaceIterator( AMP::Mesh::GeomType::Vertex, 1 ), subsetDOF );
    iterator = mesh->getSurfaceIterator( AMP::Mesh::GeomType::Vertex, 0 );
    if ( subsetDOF->numGlobalDOF() < DOF->numGlobalDOF() &&
         subsetDOF->numLocalDOF() == iterator.size() )
        ut->passes( "Subset DOF on surface mesh iterator: " + GENERATOR::name() );
    else
        ut->failure( "Subset DOF on surface mesh iterator: " + GENERATOR::name() );

    // Subset using VectorSelector
    auto submesh = mesh->Subset( mesh->getSurfaceIterator( AMP::Mesh::GeomType::Face, 1 ) );
    AMP::LinearAlgebra::MeshIteratorVariable meshSubsetVariable(
        "surfaceSubset",
        submesh->getIterator( AMP::Mesh::GeomType::Vertex, 1 ),
        submesh->getComm() );
    subsetDOF = meshSubsetVariable.getSubsetDOF( DOF );
    testGetDOFIterator( ut, submesh->getIterator( AMP::Mesh::GeomType::Vertex, 1 ), subsetDOF );

    // Test subsetting for a subset mesh
    if ( mesh->getGeomType() == AMP::Mesh::GeomType::Volume ) {
        auto surface_mesh =
            mesh->Subset( mesh->getSurfaceIterator( AMP::Mesh::GeomType::Face, 1 ) );
        subsetDOF = DOF->subset( surface_mesh );
        size_t N1 = surface_mesh->numGlobalElements( AMP::Mesh::GeomType::Vertex );
        size_t N2 = subsetDOF->numGlobalDOF();
        if ( N1 == N2 )
            ut->passes( "Subset DOF on surface mesh: " + GENERATOR::name() );
        else
            ut->failure( "Subset DOF on surface mesh: " + GENERATOR::name() );
    }
}


// Function to test the creation/destruction of a simpleDOFManager
template<class GENERATOR>
void testSimpleDOFManager( AMP::UnitTest *ut )
{
    // Get the mesh
    GENERATOR mesh_generator;
    mesh_generator.build_mesh();
    auto mesh = mesh_generator.getMesh();

    // Create some simple DOF managers
    auto DOF1 = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 0, 1, false );
    auto DOF2 = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 1, false );
    auto DOF3 = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 3, false );
    auto DOF4 = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 1, true );

    // Run basic tests that should be valid for all DOFManagers
    testBasics( DOF1, ut );
    testBasics( DOF2, ut );
    testBasics( DOF3, ut );
    testBasics( DOF4, ut );
    testSubsetComm( DOF1, ut );
    testSubsetComm( DOF2, ut );
    testSubsetComm( DOF3, ut );
    testSubsetComm( DOF4, ut );
    testSubsetMesh( mesh, DOF1, true, 1, 0, ut );
    testSubsetMesh( mesh, DOF2, true, 1, 1, ut );
    testSubsetMesh( mesh, DOF3, true, 3, 1, ut );
    testSubsetMesh( mesh, DOF4, true, 1, 1, ut );

    // Check some simple properties
    if ( DOF1->numLocalDOF() > 0 )
        ut->passes( "Non-empty DOFs: " + GENERATOR::name() );
    else
        ut->failure( "Non-empty DOFs: " + GENERATOR::name() );
    if ( DOF1->endDOF() - DOF1->beginDOF() == DOF1->numLocalDOF() )
        ut->passes( "Non-empty DOFs: " + GENERATOR::name() );
    else
        ut->failure( "Non-empty DOFs: " + GENERATOR::name() );
    if ( DOF2->numLocalDOF() == DOF1->numLocalDOF() &&
         DOF3->numLocalDOF() == 3 * DOF1->numLocalDOF() &&
         DOF4->numLocalDOF() == DOF1->numLocalDOF() && DOF2->beginDOF() == DOF1->beginDOF() &&
         DOF3->beginDOF() == 3 * DOF1->beginDOF() && DOF4->beginDOF() == DOF1->beginDOF() &&
         DOF2->endDOF() == DOF1->endDOF() && DOF3->endDOF() == 3 * DOF1->endDOF() &&
         DOF4->endDOF() == DOF1->endDOF() )
        ut->passes( "DOFs agree: " + GENERATOR::name() );
    else
        ut->failure( "DOFs agree: " + GENERATOR::name() );

    // Check that the iterator size matches the mesh
    auto meshIterator = mesh->getIterator( AMP::Mesh::GeomType::Vertex, 0 );
    auto DOF1Iterator = DOF1->getIterator();
    auto DOF2Iterator = DOF2->getIterator();
    auto DOF3Iterator = DOF3->getIterator();
    auto DOF4Iterator = DOF4->getIterator();
    if ( DOF1Iterator.size() == meshIterator.size() )
        ut->passes( "DOF1 has the correct size for the iterator: " + GENERATOR::name() );
    else
        ut->failure( "DOF1 has the correct size for the iterator: " + GENERATOR::name() );
    if ( DOF2Iterator.size() == meshIterator.size() )
        ut->passes( "DOF2 has the correct size for the iterator: " + GENERATOR::name() );
    else
        ut->failure( "DOF2 has the correct size for the iterator: " + GENERATOR::name() );
    if ( DOF3Iterator.size() == meshIterator.size() )
        ut->passes( "DOF3 has the correct size for the iterator: " + GENERATOR::name() );
    else
        ut->failure( "DOF3 has the correct size for the iterator: " + GENERATOR::name() );
    if ( DOF4Iterator.size() == meshIterator.size() )
        ut->passes( "DOF4 has the correct size for the iterator: " + GENERATOR::name() );
    else
        ut->failure( "DOF4 has the correct size for the iterator: " + GENERATOR::name() );

    // Check the iterator based constructor
    DOF1 = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 1, false );
    DOF2 = AMP::Discretization::simpleDOFManager::create(
        mesh,
        mesh->getIterator( AMP::Mesh::GeomType::Vertex, 1 ),
        mesh->getIterator( AMP::Mesh::GeomType::Vertex, 0 ),
        1 );
    if ( DOF1->numGlobalDOF() == DOF2->numGlobalDOF() )
        ut->passes( "iterator-based constructor created: " + GENERATOR::name() );
    else
        ut->failure( "iterator-based constructor created: " + GENERATOR::name() );

    // Check that we can get the DOFs
    testGetDOFIterator( ut, mesh->getIterator( AMP::Mesh::GeomType::Vertex, 0 ), DOF1 );
    testGetDOFIterator( ut, mesh->getIterator( AMP::Mesh::GeomType::Vertex, 1 ), DOF2 );
    testGetDOFIterator( ut, mesh->getIterator( AMP::Mesh::GeomType::Vertex, 1 ), DOF3 );
    testGetDOFIterator( ut, mesh->getIterator( AMP::Mesh::GeomType::Vertex, 1 ), DOF4 );
}


// Function to test the creation/destruction of a multiDOFManager
template<class GENERATOR>
void testMultiDOFManager( AMP::UnitTest *ut )
{
    // Get the mesh
    GENERATOR mesh_generator;
    mesh_generator.build_mesh();
    AMP::Mesh::Mesh::shared_ptr mesh = mesh_generator.getMesh();

    // Create a simple DOF manager and check if it is a multiDOF manager
    auto DOFs = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    if ( std::dynamic_pointer_cast<AMP::Mesh::MultiMesh>( mesh ) ) {
        auto multiDOF = std::dynamic_pointer_cast<AMP::Discretization::multiDOFManager>( DOFs );
        if ( multiDOF ) {
            ut->passes( "Created multiDOFManager from simpleDOFManager: " + GENERATOR::name() );
            testMultiDOFMap( ut, multiDOF );
        } else {
            ut->failure( "Created multiDOFManager from simpleDOFManager: " + GENERATOR::name() );
        }
    }

    // Run basic tests that should be valid for all DOFManagers
    testBasics( DOFs, ut );
    testSubsetComm( DOFs, ut );
    testSubsetMesh( mesh, DOFs, true, 1, 1, ut );

    // Test a multivector that contains multiple vectors with the same DOFManager
    testMultiDOFVector( ut, DOFs );

    // Create a multiDOFManager with repeated mesh elements and make sure the iterator only iterates
    // once through each element
    std::vector<AMP::Discretization::DOFManager::shared_ptr> managers( 2, DOFs );
    auto DOF2 = std::make_shared<AMP::Discretization::multiDOFManager>( DOFs->getComm(), managers );
    AMP::Mesh::MeshIterator iterator1 = DOFs->getIterator();
    AMP::Mesh::MeshIterator iterator2 = DOF2->getIterator();
    if ( iterator1.size() == iterator2.size() )
        ut->passes( "multiDOFManager iterates once through each element: " + GENERATOR::name() );
    else
        ut->failure( "multiDOFManager iterates once through each element: " + GENERATOR::name() );

    // Run basic tests that should be valid for all DOFManagers
    testBasics( DOF2, ut );
    testSubsetComm( DOF2, ut );
    testSubsetMesh( mesh, DOF2, true, 2, 1, ut );
}


// Function to test the creation/destruction of a multiDOFManager
template<class GENERATOR, int Nx, int Ny, int Nz, int GCW>
void testStructureDOFManager( AMP::UnitTest *ut )
{
    // Get the mesh
    GENERATOR mesh_generator;
    mesh_generator.build_mesh();
    AMP::Mesh::Mesh::shared_ptr mesh = mesh_generator.getMesh();

    // Create a simple DOF manager and check if it is a multiDOF manager
    int dofsPerFace[3] = { Nx, Ny, Nz };
    auto DOFs = AMP::Discretization::structuredFaceDOFManager::create( mesh, dofsPerFace, GCW );

    // Run basic tests that should be valid for all DOFManagers
    testBasics( DOFs, ut );
    testSubsetComm( DOFs, ut );
    testSubsetMesh( mesh, DOFs, false, 0, GCW, ut );

    // Check getRowDOFs
    if ( Nx > 0 && Ny > 0 && Nz > 0 ) {
        bool pass = true;
        auto it   = mesh->getIterator( AMP::Mesh::GeomType::Volume, 0 );
        for ( size_t ii = 0; ii < it.size(); ++ii, ++it ) {
            auto faces = it->getElements( AMP::Mesh::GeomType::Face );
            for ( size_t i = 0; i < faces.size(); i++ ) {
                std::vector<size_t> rows = DOFs->getRowDOFs( faces[i] );
                std::vector<size_t> dofs;
                for ( size_t j = 0; j < faces.size(); j++ ) {
                    DOFs->getDOFs( faces[j].globalID(), dofs );
                    for ( size_t k = 0; k < dofs.size(); k++ ) {
                        size_t index = AMP::Utilities::findfirst( rows, dofs[k] );
                        if ( index == rows.size() ) {
                            index--;
                        }
                        if ( rows[index] != dofs[k] )
                            pass = false;
                    }
                }
            }
        }
        if ( pass )
            ut->passes( "getRowDOFs found all faces that share a volume: " + GENERATOR::name() );
        else
            ut->failure( "getRowDOFs found all faces that share a volume: " + GENERATOR::name() );
    }
}
