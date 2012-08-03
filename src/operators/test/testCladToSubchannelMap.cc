#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/PIO.h"
#include "ampmesh/Mesh.h"
#include "ampmesh/SiloIO.h"
#include "ampmesh/MeshElementVectorIterator.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
//#include "operators/map/testCladToSubChannelMap.h"
#include "operators/map/AsyncMapColumnOperator.h"
#include "vectors/Variable.h"
#include "vectors/VectorBuilder.h"


double getTemp(const std::vector<double> &x) {
    return 500 + x[0]*100 + x[1]*100 + x[2]*100;
}


AMP::Mesh::MeshIterator getZFaceIterator(AMP::Mesh::Mesh::shared_ptr subChannel, int ghostWidth)
{
    std::multimap<double,AMP::Mesh::MeshElement> xyFace;
    AMP::Mesh::MeshIterator iterator = subChannel->getIterator( AMP::Mesh::Face, ghostWidth );
    for(size_t i=0; i<iterator.size(); ++i ) {
        std::vector<AMP::Mesh::MeshElement> nodes = iterator->getElements(AMP::Mesh::Vertex);
        std::vector<double> center = iterator->centroid();
        bool is_valid = true;
        for (size_t j=0; j<nodes.size(); ++j) {
            std::vector<double> coord = nodes[j].coord();
            if ( !AMP::Utilities::approx_equal(coord[2],center[2], 1e-6) )
                is_valid = false;
        }
        if ( is_valid ) {
            xyFace.insert(std::pair<double,AMP::Mesh::MeshElement>(center[2],*iterator));
        }
        ++iterator;
    }
    boost::shared_ptr<std::vector<AMP::Mesh::MeshElement> > elements( new std::vector<AMP::Mesh::MeshElement>() );
    elements->reserve(xyFace.size());
    for (std::multimap<double,AMP::Mesh::MeshElement>::iterator it=xyFace.begin(); it!=xyFace.end(); ++it)
        elements->push_back( it->second );
    return AMP::Mesh::MultiVectorIterator( elements );
}


void  runTest ( const std::string &fname , AMP::UnitTest *ut )
{
    // Read the input file
    boost::shared_ptr<AMP::InputDatabase>  input_db ( new AMP::InputDatabase ( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile ( fname , input_db );
    input_db->printClassData (AMP::plog);

    // Get the Mesh database and create the mesh parameters
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
    boost::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(mesh_db));
    params->setComm(globalComm);

    // Create the meshes from the input database
    AMP::Mesh::Mesh::shared_ptr manager = AMP::Mesh::Mesh::buildMesh(params);
    AMP::Mesh::Mesh::shared_ptr pin_mesh = manager->Subset("MultiPin");
    pin_mesh->setName("MultiPin");
    AMP::Mesh::Mesh::shared_ptr subchannel_mesh = manager->Subset("subchannel");
    subchannel_mesh->setName("subchannel");
    AMP::Mesh::Mesh::shared_ptr clad_mesh = pin_mesh->Subset("clad");
    AMP::Mesh::Mesh::shared_ptr subchannel_face = subchannel_mesh->Subset(getZFaceIterator(subchannel_mesh,1));

    // Get the database for the map
    //boost::shared_ptr<AMP::Database> map_db = input_db->getDatabase( "MeshToMeshMaps" );

    // Create the DOFManagers and the vectors
    //int DOFsPerNode = map_db->getInteger("DOFsPerObject");
    //std::string varName = map_db->getString("VariableName");
    int DOFsPerNode = 1;
    std::string varName = "Temperature";
    AMP::LinearAlgebra::Variable::shared_ptr temperature( new AMP::LinearAlgebra::Variable(varName) );
    AMP::Discretization::DOFManager::shared_ptr  pin_DOFs = 
        AMP::Discretization::simpleDOFManager::create(pin_mesh,AMP::Mesh::Vertex,1,DOFsPerNode);
    AMP::Discretization::DOFManager::shared_ptr  subchannel_DOFs = 
        AMP::Discretization::simpleDOFManager::create(subchannel_face,AMP::Mesh::Face,1,DOFsPerNode);

    // Create the vectors
    AMP::LinearAlgebra::Vector::shared_ptr  dummy;
    AMP::LinearAlgebra::Vector::shared_ptr T1 = AMP::LinearAlgebra::createVector( pin_DOFs, temperature );
    AMP::LinearAlgebra::Vector::shared_ptr T2 = AMP::LinearAlgebra::createVector( subchannel_DOFs, temperature );
    T1->setToScalar(0.0);
    T2->setToScalar(0.0);

    // Initialize the pin temperatures
    AMP::Mesh::MeshIterator it = pin_mesh->getIterator(AMP::Mesh::Vertex,0);
    std::vector<size_t> dofs;
    for (size_t i=0; i<it.size(); i++) {
        pin_DOFs->getDOFs(it->globalID(),dofs);
        T1->setValueByGlobalID(dofs[0],getTemp(it->coord()));
        ++it;
    }

    /*// Test the creation/destruction of ScalarZAxisMap (no apply call)
    try { 
        boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator>  gapmaps;
        gapmaps = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::ScalarZAxisMap> ( mesh, map_db );
        gapmaps.reset();
        ut->passes("Created / Destroyed ScalarZAxisMap");
    } catch ( ... ) {
        ut->failure("Created / Destroyed ScalarZAxisMap");
    }

    // Perform a complete test of ScalarZAxisMap
    boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator>  gapmaps;
    gapmaps = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::ScalarZAxisMap> ( mesh, map_db );
    gapmaps->setVector ( v2 );
    
    // Initialize the vectors
    v1->setToScalar(0.0);
    v2->setToScalar(0.0);
    size_t N_maps = (size_t) map_db->getInteger("N_maps");
    std::vector<std::string> mesh1 = map_db->getStringArray("Mesh1");
    std::vector<std::string> mesh2 = map_db->getStringArray("Mesh2");
    std::vector<int> surface1 = map_db->getIntegerArray("Surface1");
    std::vector<int> surface2 = map_db->getIntegerArray("Surface2");
    AMP_ASSERT(mesh1.size()==N_maps||mesh1.size()==1);
    AMP_ASSERT(mesh2.size()==N_maps||mesh2.size()==1);
    AMP_ASSERT(surface1.size()==N_maps||surface1.size()==1);
    AMP_ASSERT(surface2.size()==N_maps||surface2.size()==1);
    for (size_t i=0; i<N_maps; i++) {
        std::string meshname1,  meshname2;
        if ( mesh1.size() == N_maps ) {
            meshname1 = mesh1[i];
            meshname2 = mesh2[i];
        } else {
            meshname1 = mesh1[0];
            meshname2 = mesh2[0];
        }
        int surface_id1, surface_id2;
        if ( surface1.size() == N_maps ) {
            surface_id1 = surface1[i];
            surface_id2 = surface2[i];
        } else {
            surface_id1 = surface1[0];
            surface_id2 = surface2[0];
        }
        AMP::Mesh::Mesh::shared_ptr curMesh = mesh->Subset( meshname1 );
        setBoundary( surface_id1, v1, curMesh );
        curMesh = mesh->Subset( meshname2 );
        setBoundary( surface_id2, v1, curMesh );
    }

    // Apply the maps
    globalComm.barrier();
    gapmaps->apply ( dummy , v1 , v2 );
    v1->subtract ( v1 , v2 );
    if ( v1->maxNorm() < 1.e-12 )
        ut->passes ( "Node to node map test" );
    else
        ut->failure ( "Node to node map test" );

*/

    // Write the results
    #ifdef USE_SILO
        AMP::Mesh::SiloIO::shared_ptr  siloWriter( new AMP::Mesh::SiloIO);
        siloWriter->registerVector( T1, pin_mesh, AMP::Mesh::Vertex, "Temperature" );
        siloWriter->registerVector( T2, subchannel_face, AMP::Mesh::Face, "Temperature" );
        siloWriter->setDecomposition( 0 );
        siloWriter->writeFile( fname, 0 );
    #endif
}


int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    AMP::AMP_MPI globalComm = AMP::AMP_MPI(AMP_COMM_WORLD);
    runTest ( "inputCladToSubchannelMap-1" , &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;

}
