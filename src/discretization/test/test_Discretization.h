#include "ampmesh/Mesh.h"
#include "ampmesh/MultiMesh.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "discretization/MultiDOF_Manager.h"
#include "../../ampmesh/test/meshGenerators.h"

using namespace AMP::unit_test;


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
        for (size_t j=1; j<localDOFList.size(); j++) {
            if ( localDOFList[j] != localDOFList[j-1]+1 )
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


