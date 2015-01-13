#include "ampmesh/Mesh.h"
#include "ampmesh/MultiMesh.h"
#include "discretization/DOF_Manager.h"
#include "discretization/subsetDOFManager.h"
#include "discretization/simpleDOF_Manager.h"
#include "discretization/MultiDOF_Manager.h"
#include "discretization/structuredFaceDOFManager.h"
#include "utils/AMP_MPI.h"
#include "../../ampmesh/test/meshGenerators.h"
#include "DOFManager_tests.h"


using namespace AMP::unit_test;



// Function to test subsetting a DOF manager
template <class GENERATOR, bool SPLIT>
void testSubsetDOFManager( AMP::UnitTest *ut )
{
    // Get the mesh
    GENERATOR mesh_generator;
    mesh_generator.build_mesh();
    AMP::Mesh::Mesh::shared_ptr mesh = mesh_generator.getMesh();

    // Create a simple DOF manager
    AMP::Discretization::DOFManager::shared_ptr DOF = 
        AMP::Discretization::simpleDOFManager::create( mesh, AMP::Mesh::Vertex, 1, 1, SPLIT );
    testGetDOFIterator(  ut, mesh->getIterator(AMP::Mesh::Vertex,1), DOF );

    // Run basic tests that should be valid for all DOFManagers
    testBasics( DOF, ut );
    testSubsetComm( DOF, ut );
    testSubsetMesh( mesh, DOF, true, 1, 1, ut );

    // Subset for iterators
    AMP::Mesh::MeshIterator iterator = mesh->getIterator( AMP::Mesh::Vertex, 1 );
    AMP::Discretization::DOFManager::shared_ptr subsetDOF = DOF->subset( iterator, mesh->getComm() );
    testGetDOFIterator(  ut, mesh->getIterator(AMP::Mesh::Vertex,1), subsetDOF );
    if ( *DOF == *subsetDOF )
        ut->passes("Subset DOF on full mesh iterator: "+GENERATOR::name());
    else
        ut->failure("Subset DOF on full mesh iterator: "+GENERATOR::name());
    iterator = mesh->getSurfaceIterator( AMP::Mesh::Vertex, 1 );
    subsetDOF = DOF->subset( iterator, mesh->getComm() );
    testGetDOFIterator(  ut, mesh->getSurfaceIterator(AMP::Mesh::Vertex,1), subsetDOF );
    iterator = mesh->getSurfaceIterator( AMP::Mesh::Vertex, 0 );
    if ( subsetDOF->numGlobalDOF()<DOF->numGlobalDOF() && subsetDOF->numLocalDOF()==iterator.size() )
        ut->passes("Subset DOF on surface mesh iterator: "+GENERATOR::name());
    else
        ut->failure("Subset DOF on surface mesh iterator: "+GENERATOR::name());

    // Subset using VectorSelector
    #ifdef USE_AMP_VECTORS
        AMP::Mesh::Mesh::shared_ptr submesh = mesh->Subset( mesh->getSurfaceIterator(AMP::Mesh::Face,1) );
        AMP::LinearAlgebra::MeshIteratorVariable meshSubsetVariable( 
            "surfaceSubset", submesh->getIterator(AMP::Mesh::Vertex,1), submesh->getComm() );
        subsetDOF = meshSubsetVariable.getSubsetDOF( DOF );
        testGetDOFIterator(  ut, submesh->getIterator(AMP::Mesh::Vertex,1), subsetDOF );
    #endif

    // Test subsetting for a subset mesh
    if ( mesh->getGeomType() == AMP::Mesh::Volume ) {
        AMP::Mesh::Mesh::shared_ptr surface_mesh = mesh->Subset( mesh->getSurfaceIterator(AMP::Mesh::Face,1) );
        subsetDOF = DOF->subset( surface_mesh );
        size_t N1 = surface_mesh->numGlobalElements(AMP::Mesh::Vertex);
        size_t N2 = subsetDOF->numGlobalDOF();
        if ( N1==N2 )
            ut->passes("Subset DOF on surface mesh: "+GENERATOR::name());
        else
            ut->failure("Subset DOF on surface mesh: "+GENERATOR::name());
    }
}


// Function to test the creation/destruction of a simpleDOFManager
template <class GENERATOR>
void testSimpleDOFManager( AMP::UnitTest *ut )
{
    // Get the mesh
    GENERATOR mesh_generator;
    mesh_generator.build_mesh();
    AMP::Mesh::Mesh::shared_ptr mesh = mesh_generator.getMesh();
    
    // Create some simple DOF managers
    AMP::Discretization::DOFManager::shared_ptr DOF1 =  AMP::Discretization::simpleDOFManager::create( mesh, AMP::Mesh::Vertex, 0, 1, false );
    AMP::Discretization::DOFManager::shared_ptr DOF2 =  AMP::Discretization::simpleDOFManager::create( mesh, AMP::Mesh::Vertex, 1, 1, false );
    AMP::Discretization::DOFManager::shared_ptr DOF3 =  AMP::Discretization::simpleDOFManager::create( mesh, AMP::Mesh::Vertex, 1, 3, false );
    AMP::Discretization::DOFManager::shared_ptr DOF4 =  AMP::Discretization::simpleDOFManager::create( mesh, AMP::Mesh::Vertex, 1, 1, true );
    
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
        ut->passes("Non-empty DOFs: "+GENERATOR::name());
    else
        ut->failure("Non-empty DOFs: "+GENERATOR::name());
    if ( DOF1->endDOF()-DOF1->beginDOF() == DOF1->numLocalDOF() )
        ut->passes("Non-empty DOFs: "+GENERATOR::name());
    else
        ut->failure("Non-empty DOFs: "+GENERATOR::name());
    if ( DOF2->numLocalDOF()==DOF1->numLocalDOF() && DOF3->numLocalDOF()==3*DOF1->numLocalDOF() && DOF4->numLocalDOF()==DOF1->numLocalDOF() &&
            DOF2->beginDOF()==DOF1->beginDOF() && DOF3->beginDOF()==3*DOF1->beginDOF() && DOF4->beginDOF()==DOF1->beginDOF() && 
            DOF2->endDOF()==DOF1->endDOF() && DOF3->endDOF()==3*DOF1->endDOF() && DOF4->endDOF()==DOF1->endDOF() )
        ut->passes("DOFs agree: "+GENERATOR::name());
    else
        ut->failure("DOFs agree: "+GENERATOR::name());

    // Check that the iterator size matches the mesh
    AMP::Mesh::MeshIterator meshIterator = mesh->getIterator(AMP::Mesh::Vertex,0);
    AMP::Mesh::MeshIterator DOF1Iterator = DOF1->getIterator();
    AMP::Mesh::MeshIterator DOF2Iterator = DOF2->getIterator();
    AMP::Mesh::MeshIterator DOF3Iterator = DOF3->getIterator();
    AMP::Mesh::MeshIterator DOF4Iterator = DOF4->getIterator();
    if ( DOF1Iterator.size() == meshIterator.size() )
        ut->passes("DOF1 has the correct size for the iterator: "+GENERATOR::name());
    else
        ut->failure("DOF1 has the correct size for the iterator: "+GENERATOR::name());
    if ( DOF2Iterator.size() == meshIterator.size() )
        ut->passes("DOF2 has the correct size for the iterator: "+GENERATOR::name());
    else
        ut->failure("DOF2 has the correct size for the iterator: "+GENERATOR::name());
    if ( DOF3Iterator.size() == meshIterator.size() )
        ut->passes("DOF3 has the correct size for the iterator: "+GENERATOR::name());
    else
        ut->failure("DOF3 has the correct size for the iterator: "+GENERATOR::name());
    if ( DOF4Iterator.size() == meshIterator.size() )
        ut->passes("DOF4 has the correct size for the iterator: "+GENERATOR::name());
    else
        ut->failure("DOF4 has the correct size for the iterator: "+GENERATOR::name());

    // Check the iterator based constructor
    DOF1 = AMP::Discretization::simpleDOFManager::create( mesh, AMP::Mesh::Vertex, 1, 1, false );
    DOF2 = AMP::Discretization::simpleDOFManager::create( 
        mesh, mesh->getIterator(AMP::Mesh::Vertex,1), mesh->getIterator(AMP::Mesh::Vertex,0), 1 );
    if ( DOF1->numGlobalDOF() == DOF2->numGlobalDOF())
        ut->passes("iterator-based constructor created: "+GENERATOR::name());
    else
        ut->failure("iterator-based constructor created: "+GENERATOR::name());

    // Check that we can get the DOFs
    testGetDOFIterator(  ut, mesh->getIterator(AMP::Mesh::Vertex,0), DOF1 );
    testGetDOFIterator(  ut, mesh->getIterator(AMP::Mesh::Vertex,1), DOF2 );
    testGetDOFIterator(  ut, mesh->getIterator(AMP::Mesh::Vertex,1), DOF3 );
    testGetDOFIterator(  ut, mesh->getIterator(AMP::Mesh::Vertex,1), DOF4 );
}



// Function to test the creation/destruction of a multiDOFManager
template <class GENERATOR>
void testMultiDOFManager( AMP::UnitTest *ut )
{
    // Get the mesh
    GENERATOR mesh_generator;
    mesh_generator.build_mesh();
    AMP::Mesh::Mesh::shared_ptr mesh = mesh_generator.getMesh();
    
    // Create a simple DOF manager and check if it is a multiDOF manager
    AMP::Discretization::DOFManager::shared_ptr DOFs =  AMP::Discretization::simpleDOFManager::create( mesh, AMP::Mesh::Vertex, 1, 1, true );
    if ( AMP::dynamic_pointer_cast<AMP::Mesh::MultiMesh>(mesh).get() != NULL ) {
        AMP::shared_ptr<AMP::Discretization::multiDOFManager> multiDOF = AMP::dynamic_pointer_cast<AMP::Discretization::multiDOFManager>(DOFs);
        if ( multiDOF.get() != NULL ) {
            ut->passes("Created multiDOFManager from simpleDOFManager: "+GENERATOR::name());
            testMultiDOFMap( ut, multiDOF );
        } else {
            ut->failure("Created multiDOFManager from simpleDOFManager: "+GENERATOR::name());
        }
    }

    // Run basic tests that should be valid for all DOFManagers
    testBasics( DOFs, ut );
    testSubsetComm( DOFs, ut );
    testSubsetMesh( mesh, DOFs, true, 1, 1, ut );

    // Test a multivector that contains multiple vectors with the same DOFManager
    #ifdef USE_AMP_VECTORS
        testMultiDOFVector( ut, DOFs );
    #else
        ut->expected_failure("Can't test multivector without vectors: "+GENERATOR::name());
    #endif

    // Create a multiDOFManager with repeated mesh elements and make sure the iterator only iterates once through each element
    std::vector<AMP::Discretization::DOFManager::shared_ptr> managers(2,DOFs);
    AMP::Discretization::DOFManager::shared_ptr DOF2( new AMP::Discretization::multiDOFManager( DOFs->getComm(), managers ) );
    AMP::Mesh::MeshIterator iterator1 = DOFs->getIterator();
    AMP::Mesh::MeshIterator iterator2 = DOF2->getIterator();
    if ( iterator1.size()==iterator2.size() )
        ut->passes("multiDOFManager iterates once through each element: "+GENERATOR::name());
    else
        ut->failure("multiDOFManager iterates once through each element: "+GENERATOR::name());

    // Run basic tests that should be valid for all DOFManagers
    testBasics( DOF2, ut );
    testSubsetComm( DOF2, ut );
    testSubsetMesh( mesh, DOF2, true, 2, 1, ut );
}


// Function to test the creation/destruction of a multiDOFManager
template < class GENERATOR, int Nx, int Ny, int Nz, int GCW >
void testStructureDOFManager( AMP::UnitTest *ut )
{
    // Get the mesh
    GENERATOR mesh_generator;
    mesh_generator.build_mesh();
    AMP::Mesh::Mesh::shared_ptr mesh = mesh_generator.getMesh();
    
    // Create a simple DOF manager and check if it is a multiDOF manager
    int dofsPerFace[3]={Nx,Ny,Nz};
    AMP::Discretization::DOFManager::shared_ptr DOFs = AMP::Discretization::structuredFaceDOFManager::create( mesh, dofsPerFace, GCW );

    // Run basic tests that should be valid for all DOFManagers
    testBasics( DOFs, ut );
    testSubsetComm( DOFs, ut );
    testSubsetMesh( mesh, DOFs, false, 0, GCW, ut );

    // Check getRowDOFs
    if ( Nx>0 && Ny>0 && Nz>0 ) {
        bool pass = true;
        AMP::Mesh::MeshIterator it = mesh->getIterator(AMP::Mesh::Volume,0);
        for (size_t ii=0; ii<it.size(); ++ii, ++it) {
            std::vector<AMP::Mesh::MeshElement> faces = it->getElements(AMP::Mesh::Face);
            for (size_t i=0; i<faces.size(); i++) {
                std::vector<size_t> row = DOFs->getRowDOFs(faces[i]);
                std::vector<size_t> dofs;
                for (size_t j=0; j<faces.size(); j++) {
                    DOFs->getDOFs(faces[j].globalID(),dofs);
                    for (size_t k=0; k<dofs.size(); k++) {
                        size_t index = AMP::Utilities::findfirst(row,dofs[k]);
                        if ( index == row.size() ) { index--; }
                        if ( row[index] != dofs[k] )
                            pass = false;
                    }
                }
            }
        }
        if ( pass )
            ut->passes("getRowDOFs found all faces that share a volume: "+GENERATOR::name());
        else
            ut->failure("getRowDOFs found all faces that share a volume: "+GENERATOR::name());
    }
}



