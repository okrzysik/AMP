#include "ampmesh/Mesh.h"
#include "ampmesh/MultiMesh.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "discretization/MultiDOF_Manager.h"
#include "../../ampmesh/test/meshGenerators.h"

using namespace AMP::unit_test;


// Function to test subsetting a DOF manager
template <class GENERATOR>
void testSubsetDOFManager( AMP::UnitTest *ut )
{
    // Get the mesh
    GENERATOR mesh_generator;
    mesh_generator.build_mesh();
    AMP::Mesh::Mesh::shared_ptr mesh = mesh_generator.getMesh();

    // Create a simple DOF manager
    AMP::Discretization::DOFManager::shared_ptr DOF =  AMP::Discretization::simpleDOFManager::create( mesh, AMP::Mesh::Vertex, 0, 1, false );
    
    // Subset for each mesh
    AMP::Discretization::DOFManager::shared_ptr subsetDOF = DOF->subset( mesh );
    if ( DOF->numGlobalDOF() == subsetDOF->numGlobalDOF() )
        ut->passes("Subset DOF on full mesh");
    else
        ut->failure("Subset DOF on full mesh");
    std::vector<AMP::Mesh::MeshID> meshIDs = mesh->getBaseMeshIDs();
    if ( meshIDs.size() > 1 ) {
        size_t tot_size = 0;
        for (size_t i=0; i<meshIDs.size(); i++) {
            AMP::Mesh::Mesh::shared_ptr subsetMesh = mesh->Subset(meshIDs[i]);
            if ( subsetMesh.get() != NULL ) {
                subsetDOF = DOF->subset(subsetMesh);
                tot_size += subsetDOF->numLocalDOF();
            }
        }
        tot_size = DOF->getComm().sumReduce(tot_size);
        if ( tot_size == DOF->numGlobalDOF() )
            ut->passes("Subset DOF for each mesh");
        else
            ut->failure("Subset DOF for each mesh");
    }

    // Subset for iterators
    AMP::Mesh::MeshIterator iterator = mesh->getIterator( AMP::Mesh::Vertex, 0 );
    subsetDOF = DOF->subset( iterator, mesh->getComm() );
    if ( DOF->numGlobalDOF() == subsetDOF->numGlobalDOF() )
        ut->passes("Subset DOF on full mesh iterator");
    else
        ut->failure("Subset DOF on full mesh iterator");
    iterator = mesh->getSurfaceIterator( AMP::Mesh::Vertex, 0 );
    subsetDOF = DOF->subset( iterator, mesh->getComm() );
    if ( subsetDOF->numGlobalDOF()<DOF->numGlobalDOF() && subsetDOF->numGlobalDOF()>0 )
        ut->passes("Subset DOF on surface mesh iterator");
    else
        ut->failure("Subset DOF on surface mesh iterator");

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
    
    // Check some simple properties
    if ( DOF1->numLocalDOF() > 0 )
        ut->passes("Non-empty DOFs");
    else
        ut->failure("Non-empty DOFs");
    if ( DOF1->endDOF()-DOF1->beginDOF() == DOF1->numLocalDOF() )
        ut->passes("Non-empty DOFs");
    else
        ut->failure("Non-empty DOFs");
    if ( DOF2->numLocalDOF()==DOF1->numLocalDOF() && DOF3->numLocalDOF()==3*DOF1->numLocalDOF() && DOF4->numLocalDOF()==DOF1->numLocalDOF() &&
            DOF2->beginDOF()==DOF1->beginDOF() && DOF3->beginDOF()==3*DOF1->beginDOF() && DOF4->beginDOF()==DOF1->beginDOF() && 
            DOF2->endDOF()==DOF1->endDOF() && DOF3->endDOF()==3*DOF1->endDOF() && DOF4->endDOF()==DOF1->endDOF() )
        ut->passes("DOFs agree");
    else
        ut->failure("DOFs agree");

    // Check the iterator based constructor
    DOF1 =  AMP::Discretization::simpleDOFManager::create( mesh, AMP::Mesh::Vertex, 1, 1, false );
    DOF2 =  AMP::Discretization::simpleDOFManager::create( mesh, mesh->getIterator(AMP::Mesh::Vertex,1), mesh->getIterator(AMP::Mesh::Vertex,0), 1 );
    if ( DOF1->numGlobalDOF() == DOF2->numGlobalDOF())
        ut->passes("iterator-based constructor created");
    else
        ut->failure("iterator-based constructor created");
}


// Function to test the index conversion given a multDOFManager
void testMultiDOFMap( AMP::UnitTest *ut, boost::shared_ptr<AMP::Discretization::multiDOFManager> multiDOF )
{
    // First create a global DOF list
    size_t N_global = multiDOF->numGlobalDOF();
    std::vector<size_t> globalDOFList(N_global);
    for (size_t i=0; i<N_global; i++)
        globalDOFList[i] = i;

    // Loop through the DOFManagers
    std::vector<AMP::Discretization::DOFManager::shared_ptr> managers = multiDOF->getDOFManagers();
    for (size_t i=0; i<managers.size(); i++) {
        AMP::AMP_MPI comm = managers[i]->getComm();
        comm.barrier();
        // Convert the global list to a local list
        std::vector<size_t> tmp = multiDOF->getSubDOF( managers[i], globalDOFList );
        std::vector<size_t> localDOFList, localToGlobal;
        localDOFList.reserve(N_global);
        localToGlobal.reserve(N_global);
        size_t neg_one = ~((size_t)0);
        for (size_t j=0; j<N_global; j++) {
            if ( tmp[j] != neg_one ) {
                localDOFList.push_back(tmp[j]);
                localToGlobal.push_back(j);
            }
        }
        // Check that we created the proper list
        bool passes = localDOFList.size()==managers[i]->numGlobalDOF();
        for (size_t j=0; j<localDOFList.size(); j++) {
            if ( localDOFList[j] != j )
                passes = false;
        }
        if ( passes )
            ut->passes("Conversion from global to sub DOFs in multiDOFManager");
        else
            ut->failure("Conversion from global to sub DOFs in multiDOFManager");
        // Check that we can convert back to the global ids
        std::vector<size_t> globalDOFList2 = multiDOF->getGlobalDOF( managers[i], localDOFList );
        passes = globalDOFList2.size()==localDOFList.size();
        for (size_t j=0; j<globalDOFList2.size(); j++) {
            if ( globalDOFList2[j] != localToGlobal[j] )
                passes = false;
        }
        if ( passes )
            ut->passes("Conversion from sub to global DOFs in multiDOFManager");
        else
            ut->failure("Conversion from sub to global DOFs in multiDOFManager");
    }
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
    if ( boost::dynamic_pointer_cast<AMP::Mesh::MultiMesh>(mesh).get() != NULL ) {
        boost::shared_ptr<AMP::Discretization::multiDOFManager> multiDOF = boost::dynamic_pointer_cast<AMP::Discretization::multiDOFManager>(DOFs);
        if ( multiDOF.get() != NULL ) {
            ut->passes("Created multiDOFManager from simpleDOFManager");
            testMultiDOFMap( ut, multiDOF );
        } else {
            ut->failure("Created multiDOFManager from simpleDOFManager");
        }
    }

}


