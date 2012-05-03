#include "ampmesh/Mesh.h"
#include "ampmesh/MultiMesh.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "discretization/MultiDOF_Manager.h"
#include "utils/AMP_MPI.h"
#include "../../ampmesh/test/meshGenerators.h"

#ifdef USE_AMP_VECTORS
    #include "vectors/Variable.h"
    #include "vectors/Vector.h"
    #include "vectors/MultiVector.h"
    #include "vectors/VectorBuilder.h"
#endif


using namespace AMP::unit_test;


// Function to test getting the DOFs for a mesh iterator
void testGetDOFIterator(  AMP::UnitTest *ut, const AMP::Mesh::MeshIterator &iterator, AMP::Discretization::DOFManager::shared_ptr DOF )
{
    bool passes = true;
    AMP::Mesh::MeshIterator cur = iterator.begin();
    std::vector<size_t> dofs;
    for (size_t i=0; i<cur.size(); i++) {
        DOF->getDOFs( cur->globalID(), dofs );
        if ( dofs.empty() )
            passes = false;
    }
    if ( passes ) 
        ut->passes("Got the DOFs for every element in iterator");
    else
        ut->failure("Got the DOFs for every element in iterator");
}


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
        AMP::Discretization::simpleDOFManager::create( mesh, AMP::Mesh::Vertex, 0, 1, SPLIT );
    testGetDOFIterator(  ut, mesh->getIterator(AMP::Mesh::Vertex,1), DOF );
    AMP::Discretization::DOFManager::shared_ptr subsetDOF;

    // Subset for comms
    subsetDOF = DOF->subset( AMP::AMP_MPI(AMP_COMM_WORLD) );
    if ( *DOF == *subsetDOF )
        ut->passes("Subset DOF for COMM_WORLD");
    else
        ut->failure("Subset DOF on COMM_WORLD");
    subsetDOF = DOF->subset( DOF->getComm() );
    if ( *DOF == *subsetDOF )
        ut->passes("Subset DOF for DOF comm");
    else
        ut->failure("Subset DOF on DOF comm");
    subsetDOF = DOF->subset( AMP::AMP_MPI(AMP_COMM_SELF) );
    if ( subsetDOF->numLocalDOF()==DOF->numLocalDOF() && subsetDOF->numGlobalDOF()==DOF->numLocalDOF() )
        ut->passes("Subset DOF for COMM_SELF");
    else
        ut->failure("Subset DOF on COMM_SELF");
    
    // Subset for each mesh
    subsetDOF = DOF->subset( mesh );
    if ( *DOF == *subsetDOF )
        ut->passes("Subset DOF on full mesh");
    else
        ut->failure("Subset DOF on full mesh");
    std::vector<AMP::Mesh::MeshID> meshIDs = mesh->getBaseMeshIDs();
    if ( meshIDs.size() > 1 ) {
        bool passes = true;
        for (size_t i=0; i<meshIDs.size(); i++) {
            AMP::Mesh::Mesh::shared_ptr subsetMesh = mesh->Subset(meshIDs[i]);
            if ( subsetMesh.get() != NULL ) {
                subsetDOF = DOF->subset(subsetMesh);
                testGetDOFIterator(  ut, subsetMesh->getIterator(AMP::Mesh::Vertex,1), subsetDOF );
                AMP::Discretization::DOFManager::shared_ptr mesh_DOF = 
                    AMP::Discretization::simpleDOFManager::create( subsetMesh, AMP::Mesh::Vertex, 0, 1, false );
                if ( *mesh_DOF != *subsetDOF )
                    passes = false;
            }
        }
        if ( passes )
            ut->passes("Subset DOF for each mesh");
        else
            ut->failure("Subset DOF for each mesh");
    }

    // Subset for iterators
    AMP::Mesh::MeshIterator iterator = mesh->getIterator( AMP::Mesh::Vertex, 0 );
    subsetDOF = DOF->subset( iterator, mesh->getComm() );
    testGetDOFIterator(  ut, mesh->getIterator(AMP::Mesh::Vertex,1), subsetDOF );
    if ( *DOF == *subsetDOF )
        ut->passes("Subset DOF on full mesh iterator");
    else
        ut->failure("Subset DOF on full mesh iterator");
    iterator = mesh->getSurfaceIterator( AMP::Mesh::Vertex, 0 );
    subsetDOF = DOF->subset( iterator, mesh->getComm() );
    testGetDOFIterator(  ut, mesh->getSurfaceIterator(AMP::Mesh::Vertex,1), subsetDOF );
    if ( subsetDOF->numGlobalDOF()<DOF->numGlobalDOF() && subsetDOF->numLocalDOF()==iterator.size() )
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

    // Check that the iterator size matches the mesh
    AMP::Mesh::MeshIterator meshIterator = mesh->getIterator(AMP::Mesh::Vertex,0);
    AMP::Mesh::MeshIterator DOF1Iterator = DOF1->getIterator();
    AMP::Mesh::MeshIterator DOF2Iterator = DOF2->getIterator();
    AMP::Mesh::MeshIterator DOF3Iterator = DOF3->getIterator();
    AMP::Mesh::MeshIterator DOF4Iterator = DOF4->getIterator();
    if ( DOF1Iterator.size() == meshIterator.size() )
        ut->passes("DOF1 has the correct size for the iterator");
    else
        ut->failure("DOF1 has the correct size for the iterator");
    if ( DOF2Iterator.size() == meshIterator.size() )
        ut->passes("DOF2 has the correct size for the iterator");
    else
        ut->failure("DOF2 has the correct size for the iterator");
    if ( DOF3Iterator.size() == meshIterator.size() )
        ut->passes("DOF3 has the correct size for the iterator");
    else
        ut->failure("DOF3 has the correct size for the iterator");
    if ( DOF4Iterator.size() == meshIterator.size() )
        ut->passes("DOF4 has the correct size for the iterator");
    else
        ut->failure("DOF4 has the correct size for the iterator");

    // Check the iterator based constructor
    DOF1 =  AMP::Discretization::simpleDOFManager::create( mesh, AMP::Mesh::Vertex, 1, 1, false );
    DOF2 =  AMP::Discretization::simpleDOFManager::create( mesh, mesh->getIterator(AMP::Mesh::Vertex,1), mesh->getIterator(AMP::Mesh::Vertex,0), 1 );
    if ( DOF1->numGlobalDOF() == DOF2->numGlobalDOF())
        ut->passes("iterator-based constructor created");
    else
        ut->failure("iterator-based constructor created");

    // Check that we can get the DOFs
    testGetDOFIterator(  ut, mesh->getIterator(AMP::Mesh::Vertex,0), DOF1 );
    testGetDOFIterator(  ut, mesh->getIterator(AMP::Mesh::Vertex,1), DOF2 );
    testGetDOFIterator(  ut, mesh->getIterator(AMP::Mesh::Vertex,1), DOF3 );
    testGetDOFIterator(  ut, mesh->getIterator(AMP::Mesh::Vertex,1), DOF4 );
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
        std::vector<size_t> tmp = multiDOF->getSubDOF( i, globalDOFList );
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
        std::vector<size_t> globalDOFList2 = multiDOF->getGlobalDOF( i, localDOFList );
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


// Function to test that a multivector with a DOFManager repeated correctly sets the values
#ifdef USE_AMP_VECTORS
void testMultiDOFVector( AMP::UnitTest *ut, AMP::Discretization::DOFManager::shared_ptr DOF )
{
    // Create the individual vectors
    AMP::LinearAlgebra::Variable::shared_ptr var1( new AMP::LinearAlgebra::Variable("a") );
    AMP::LinearAlgebra::Variable::shared_ptr var2( new AMP::LinearAlgebra::Variable("b") );
    AMP::LinearAlgebra::Variable::shared_ptr var3( new AMP::LinearAlgebra::Variable("c") );
    AMP::LinearAlgebra::Vector::shared_ptr vec1 = AMP::LinearAlgebra::createVector( DOF, var1, true );
    AMP::LinearAlgebra::Vector::shared_ptr vec2 = AMP::LinearAlgebra::createVector( DOF, var2, true );
    // Create the multivector
    boost::shared_ptr<AMP::LinearAlgebra::MultiVector> multiVector = AMP::LinearAlgebra::MultiVector::create( var3, DOF->getComm() );
    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> vectors(2);
    vectors[0] = vec1;
    vectors[1] = vec2;
    multiVector->addVector(vectors);
    AMP::Discretization::DOFManager::shared_ptr multiDOF = multiVector->getDOFManager();
    // Check that we can set each value correctly
    multiVector->zero();
    multiVector->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    AMP::Mesh::MeshIterator it = DOF->getIterator();
    std::vector<size_t> dof1, dof2;
    bool uniqueMultiDOFs = true;
    for (size_t i=0; i<it.size(); i++) {
        DOF->getDOFs( it->globalID(), dof1 );
        multiDOF->getDOFs( it->globalID(), dof2 );
        AMP_ASSERT(dof1.size()==1);
        AMP_ASSERT(dof2.size()==2);
        if ( dof2[0]==dof2[1] )
            uniqueMultiDOFs = false;
        multiVector->setValueByGlobalID( dof2[0], dof1[0] );
        multiVector->setValueByGlobalID( dof2[1], 2*dof1[0] );
        ++it;
    }
    if ( uniqueMultiDOFs )
        ut->passes("MultiDOFManger with duplicate subDOFManagers returns unique DOFs");
    else
        ut->failure("MultiDOFManger with duplicate subDOFManagers returns unique DOFs");
    multiVector->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    double vec1norm = vec1->L1Norm();
    double vec2norm = vec2->L1Norm();
    double multiVectorNorm = multiVector->L1Norm();
    double N_tot = (DOF->numGlobalDOF()*(DOF->numGlobalDOF()-1))/2;
    if ( vec1norm==N_tot && vec2norm==2*N_tot && multiVectorNorm==3*N_tot )
        ut->passes("MultiVector with repeated DOFs sets values correctly");
    else
        ut->failure("MultiVector with repeated DOFs sets values correctly");        
}
#endif


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

    // Test a multivector that contains multiple vectors with the same DOFManager
    #ifdef USE_AMP_VECTORS
        testMultiDOFVector( ut, DOFs );
    #else
        ut->expected_failure("Can't test multivector without vectors");
    #endif

    // Create a multiDOFManager with repeated mesh elements and make sure the iterator only iterates once through each element
    std::vector<AMP::Discretization::DOFManager::shared_ptr> managers(2,DOFs);
    AMP::Discretization::DOFManager::shared_ptr DOF2( new AMP::Discretization::multiDOFManager( DOFs->getComm(), managers ) );
    AMP::Mesh::MeshIterator iterator1 = DOFs->getIterator();
    AMP::Mesh::MeshIterator iterator2 = DOF2->getIterator();
    if ( iterator1.size()==iterator2.size() )
        ut->passes("multiDOFManager iterates once through each element");
    else
        ut->failure("multiDOFManager iterates once through each element");
}



