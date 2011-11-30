#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/PIO.h"
#include "ampmesh/Mesh.h"
#include "discretization/DOF_Manager.h"
#include "operators/map/NodeToNodeMap.h"
#include "operators/map/AsyncMapColumnOperator.h"



void  setBoundary ( size_t whichBnd , AMP::LinearAlgebra::Vector::shared_ptr &v1, AMP::Mesh::Mesh::shared_ptr mesh )
{
    AMP::LinearAlgebra::Vector::shared_ptr  vv1 = v1->castTo<AMP::LinearAlgebra::MultiVector>().getVector(0);
    AMP::Discretization::DOFManager::shared_ptr  d1 = vv1->getDOFManager();

    //AMP::DOFMap::shared_ptr  d2 = mesh->getDOFMap ( v1->castTo<AMP::LinearAlgebra::MultiVector>().getVector ( 1 )->getVariable() );
    AMP::Mesh::MeshIterator  curBnd = mesh->beginOwnedBoundary ( mesh->getMapBoundaryId ( whichBnd ) );
    AMP::Mesh::MeshIterator  endBnd = curBnd->end();

    while ( curBnd != endBnd ) {
        vv1->setValueByGlobalID( d1->getGlobalID ( curBnd->globalID() , 0 ) , curBnd->x() );
        // v1->castTo<AMP::LinearAlgebra::MultiVector>().getVector ( 1 )->setValueByGlobalID( d2->getGlobalID ( curBnd->globalID() , 0 ) , curBnd->x() );
        ++curBnd;
    }
}


void  runTest ( const std::string &fname , AMP::UnitTest *ut )
{
    boost::shared_ptr<AMP::InputDatabase>  input_db ( new AMP::InputDatabase ( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile ( fname , input_db );
    input_db->printClassData (AMP::plog);
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

    AMP::Mesh::MeshManagerParameters::shared_ptr  meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
    AMP::Mesh::MeshManager::shared_ptr  manager ( new AMP::Mesh::MeshManager ( meshmgrParams ) );
    globalComm.barrier();

    // Test the creation/destruction of NodeToNodeMap (no apply call)
    try { 
        boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator>  n2nmaps;
        n2nmaps = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::NodeToNodeMap> ( manager , input_db  );
        n2nmaps.reset();
        ut->passes("Created / Destroyed NodeToNodeMap ("+fname+")");
    } catch ( ... ) {
        ut->failure("Created / Destroyed NodeToNodeMap ("+fname+")");
    }

    // Perform a complete test of NodeToNodeMap
    boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator>  n2nmaps;
    n2nmaps = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::NodeToNodeMap> ( manager , input_db  );

    AMP::LinearAlgebra::Vector::shared_ptr  dummy;
    AMP::LinearAlgebra::Vector::shared_ptr  v1 = manager->createVector ( n2nmaps->getOutputVariable() );
    AMP::LinearAlgebra::Vector::shared_ptr  v2 = manager->createVector ( n2nmaps->getOutputVariable() );
    n2nmaps->setVector ( v2 );

    for ( size_t i = 0 ; i != manager->getNumMaps() ; i++ ) {
        AMP::Mesh::MeshManager::Adapter::shared_ptr curMesh = manager->getMesh ( manager->getMapMeshName ( i ) );
        AMP::LinearAlgebra::Vector::shared_ptr tv1 = v1->select ( AMP::Mesh::VS_ByMeshTmpl<AMP::Mesh::MeshManager::Adapter> ( curMesh ) , "tt" );
        setBoundary ( i , tv1 , manager , curMesh );
    }

    std::cout << v1->maxNorm() << std::endl;
    n2nmaps->apply ( dummy , v1 , v2 );
    std::cout << v2->maxNorm() << std::endl;
    v1->subtract ( v1 , v2 );
    if ( v1->maxNorm() < 1.e-12 ) {
        ut->passes("Node to node map test ("+fname+")");
    } else {
        ut->failure("Node to node map test ("+fname+")");
    }
    std::cout << v1->maxNorm() << std::endl;
}


int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    if ( argc<=1 ) {
        ut.failure("No input files specified");
    } else {
        for (int i=1; i<argc; i++) {
            std::cout << "Running test with input file: " << argv[i] << std::endl;
            runTest ( argv[i] , &ut );
        }
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}   

