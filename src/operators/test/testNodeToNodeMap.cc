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
#include "discretization/simpleDOF_Manager.h"
#include "operators/map/NodeToNodeMap.h"
#include "operators/map/AsyncMapColumnOperator.h"



void  setBoundary ( int id , AMP::LinearAlgebra::Vector::shared_ptr &v1, AMP::Mesh::Mesh::shared_ptr mesh )
{
    AMP::LinearAlgebra::Vector::shared_ptr  vv1 = v1->castTo<AMP::LinearAlgebra::MultiVector>().getVector(0);
    AMP::Discretization::DOFManager::shared_ptr  d1 = vv1->getDOFManager();

    AMP::Mesh::MeshIterator  curBnd = mesh->getIDsetIterator( AMP::Mesh::Vertex, id, 0 );
    AMP::Mesh::MeshIterator  endBnd = curBnd.end();

    std::vector<unsigned int> ids;
    while ( curBnd != endBnd ) {
        d1->getDOFs( *curBnd, ids );
        std::vector<double> x = curBnd->coord();
        vv1->setLocalValuesByGlobalID( ids.size(), (int*)&ids[0], &x[0] );
        ++curBnd;
    }
}


void  runTest ( const std::string &fname , AMP::UnitTest *ut )
{

    // Read the input file
    boost::shared_ptr<AMP::InputDatabase>  input_db ( new AMP::InputDatabase ( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile ( fname , input_db );
    input_db->printClassData (AMP::plog);

    // Get the Mesh database and create the mesh parameters
    boost::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(mesh_db));
    params->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));

    // Create the meshes from the input database
    boost::shared_ptr<AMP::Mesh::Mesh> mesh = AMP::Mesh::Mesh::buildMesh(params);

    // Get the database for the node to node maps
    boost::shared_ptr<AMP::Database> map_db = input_db->getDatabase( "NodeToNodeMaps" );


    // Test the creation/destruction of NodeToNodeMap (no apply call)
    try { 
        boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator>  n2nmaps;
        n2nmaps = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::NodeToNodeMap> ( mesh , map_db  );
        n2nmaps.reset();
        ut->passes("Created / Destroyed NodeToNodeMap ("+fname+")");
    } catch ( ... ) {
        ut->failure("Created / Destroyed NodeToNodeMap ("+fname+")");
    }


    // Perform a complete test of NodeToNodeMap
    boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator>  n2nmaps;
    n2nmaps = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::NodeToNodeMap> ( mesh , map_db  );

    // Create a simple DOFManager and the vectors
    AMP::Discretization::DOFManagerParameters::shared_ptr DOFparams( new AMP::Discretization::DOFManagerParameters(mesh) );
    boost::shared_ptr<AMP::Discretization::simpleDOFManager> DOFs( new AMP::Discretization::simpleDOFManager(mesh,AMP::Mesh::Vertex,1,3) );
    AMP::LinearAlgebra::Vector::shared_ptr v1 = DOFs->createVector ( n2nmaps->getOutputVariable() );
    AMP::LinearAlgebra::Vector::shared_ptr v2 = DOFs->createVector ( n2nmaps->getOutputVariable() );
    n2nmaps->setVector ( v2 );

/*    for (int i=0; i<n2nmaps->getNumberOfOperators(); i++) {
        boost::shared_ptr<AMP::Operator::AsyncMapOperator>  map = getOperator(i);
        AMP::Mesh::Mesh::shared_ptr curMesh = manager->getMesh ( manager->getMapMeshName ( i ) );
        AMP::LinearAlgebra::Vector::shared_ptr tv1 = v1->select ( AMP::Mesh::VS_ByMeshTmpl<AMP::Mesh::MeshManager::Adapter> ( curMesh ) , "tt" );
        setBoundary( id , tv1 , manager , curMesh );
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
    std::cout << v1->maxNorm() << std::endl;*/
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

