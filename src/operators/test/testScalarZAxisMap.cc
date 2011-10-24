#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/PIO.h"

#include "operators/map/ScalarZAxisMap.h"
#include "operators/map/AsyncMapColumnOperator.h"

#include "petscsys.h"


void  setBoundary ( size_t whichBnd , AMP::LinearAlgebra::Vector::shared_ptr &v1 , AMP::Mesh::MeshManager::shared_ptr &manager )
{
    AMP::Mesh::MeshManager::Adapter::shared_ptr  mesh = manager->getMesh();
    AMP::Mesh::DOFMap::shared_ptr  d1 = mesh->getDOFMap ( v1->castTo<AMP::LinearAlgebra::MultiVector>().getVector ( 0 )->getVariable() );
    // AMP::Mesh::DOFMap::shared_ptr  d2 = mesh->getDOFMap ( v1->castTo<AMP::LinearAlgebra::MultiVector>().getVector ( 1 )->getVariable() );
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator  curBnd = mesh->beginOwnedBoundary ( manager->getMapBoundaryId ( whichBnd ) );
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator  endBnd = mesh->endOwnedBoundary   ( manager->getMapBoundaryId ( whichBnd ) );

    while ( curBnd != endBnd )
    {
        v1->castTo<AMP::LinearAlgebra::MultiVector>().getVector ( 0 )->setValueByGlobalID( d1->getGlobalID ( curBnd->globalID() , 0 ) , curBnd->z() );
        // v1->castTo<AMP::LinearAlgebra::MultiVector>().getVector ( 1 )->setValueByGlobalID( d2->getGlobalID ( curBnd->globalID() , 0 ) , curBnd->x() );
        curBnd++;
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

    AMP::Mesh::MeshManager::Adapter::shared_ptr  mesh = manager->getMesh ();

    // Test the creation/destruction of ScalarZAxisMap (no apply call)
    try { 
        boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator>  gapmaps;
        gapmaps = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::ScalarZAxisMap> ( manager , input_db );
        gapmaps.reset();
        ut->passes("Created / Destroyed ScalarZAxisMap");
    } catch ( ... ) {
        ut->failure("Created / Destroyed ScalarZAxisMap");
    }
    
    // Perform a complete test of ScalarZAxisMap
    boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator>  gapmaps;
    gapmaps = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::ScalarZAxisMap> ( manager , input_db );

    AMP::LinearAlgebra::Vector::shared_ptr  dummy;
    AMP::LinearAlgebra::Vector::shared_ptr  v1 = mesh->createVector ( gapmaps->getOutputVariable() );
    AMP::LinearAlgebra::Vector::shared_ptr  v2 = mesh->createVector ( gapmaps->getOutputVariable() );
    gapmaps->setVector ( v2 );
    for ( size_t i = 0 ; i != manager->getNumMaps() ; i++ ) {
        setBoundary ( i , v1 , manager );
    }

    gapmaps->apply ( dummy , v1 , v2 );
    v1->subtract ( v1 , v2 );
    if ( v1->maxNorm() < 1.e-12 )
        ut->passes ( "Node to node map test" );
    else
        ut->failure ( "Node to node map test" );

    gapmaps.reset();
}


int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    AMP::AMP_MPI globalComm = AMP::AMP_MPI(AMP_COMM_WORLD);
    int  numNodes = globalComm.getSize();
    if ( numNodes > 1 ) {
        runTest ( "inputScalarZAxisMap-1" , &ut );
    } else {
        ut.expected_failure("testGapMap does not work for serial runs (designed for domain decomposition 1)");
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;

}
