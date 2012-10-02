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
#include "operators/map/SubchannelToCladMap.h"
#include "operators/map/SubchannelToCladGPMap.h"
#include "operators/map/AsyncMapColumnOperator.h"
#include "vectors/Variable.h"
#include "vectors/VectorBuilder.h"


double getTemp(const std::vector<double> &x) {
    return 500 + x[2]*100;
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
    AMP::Mesh::Mesh::shared_ptr clad_mesh;
    if ( pin_mesh.get()!=NULL ) {
        pin_mesh->setName("MultiPin");
        clad_mesh = pin_mesh->Subset("clad");
    }
    AMP::Mesh::Mesh::shared_ptr subchannel_mesh = manager->Subset("subchannel");
    AMP::Mesh::Mesh::shared_ptr subchannel_face;
    if ( subchannel_mesh.get()!=NULL ) {
        subchannel_mesh->setName("subchannel");
        subchannel_face = subchannel_mesh->Subset(getZFaceIterator(subchannel_mesh,1));
    }

    // Get the database for the map
    boost::shared_ptr<AMP::Database> nodal_map_db = input_db->getDatabase( "SubchannelToNodeMap" );
    boost::shared_ptr<AMP::Database> gauss_map_db = input_db->getDatabase( "SubchannelToGPMap" );

    // Create the DOFManagers and the vectors
    //int DOFsPerNode = map_db->getInteger("DOFsPerObject");
    //std::string varName = map_db->getString("VariableName");
    int DOFsPerNode = 1;
    std::string varName = "Temperature";
    AMP::LinearAlgebra::Variable::shared_ptr temperature( new AMP::LinearAlgebra::Variable(varName) );
    AMP::Discretization::DOFManager::shared_ptr  pin_DOFs;
    AMP::Discretization::DOFManager::shared_ptr  subchannel_DOFs;
    AMP::LinearAlgebra::Vector::shared_ptr T1;
    AMP::LinearAlgebra::Vector::shared_ptr T2;
    AMP::LinearAlgebra::Vector::shared_ptr dummy;
    if ( pin_mesh.get()!=NULL ) {
        pin_DOFs = AMP::Discretization::simpleDOFManager::create(pin_mesh,AMP::Mesh::Vertex,1,DOFsPerNode);
        T1 = AMP::LinearAlgebra::createVector( pin_DOFs, temperature );
        T1->setToScalar(500);
    }
    if ( subchannel_face.get()!=NULL ) {
        subchannel_DOFs = AMP::Discretization::simpleDOFManager::create(subchannel_face,AMP::Mesh::Face,1,DOFsPerNode);
        T2 = AMP::LinearAlgebra::createVector( subchannel_DOFs, temperature );
        T2->setToScalar(0.0);
    }

    // Initialize the subchannel temperatures
    if ( subchannel_face.get()!=NULL ) {
        AMP::Mesh::MeshIterator it = subchannel_face->getIterator(AMP::Mesh::Face,0);
        std::vector<size_t> dofs;
        for (size_t i=0; i<it.size(); i++) {
            subchannel_DOFs->getDOFs(it->globalID(),dofs);
            T2->setValueByGlobalID(dofs[0],getTemp(it->centroid()));
            ++it;
        }
    }

    // Test the creation/destruction of SubchannelToCladMap (no apply call)
    try { 
        boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator>  map;
        map = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::SubchannelToCladMap>( manager, nodal_map_db );
        map.reset();
        ut->passes("Created / Destroyed SubchannelToCladMap");
    } catch ( ... ) {
        ut->failure("Created / Destroyed SubchannelToCladMap");
    }
    try { 
        boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator>  map;
        map = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::SubchannelToCladGPMap>( manager, gauss_map_db );
        map.reset();
        ut->passes("Created / Destroyed SubchannelToCladGPMap");
    } catch ( ... ) {
        ut->failure("Created / Destroyed SubchannelToCladGPMap");
    }

    // Perform a complete test of SubchannelToCladMap
    boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator>  map;
    map = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::SubchannelToCladMap>( manager, nodal_map_db );
    map->setVector( T1 );
    
    // Apply the map
    globalComm.barrier();
    map->apply( dummy, T2, dummy );

    // Check the results
    if ( subchannel_face.get()!=NULL ) {
        bool passes = true;
        AMP::Mesh::MeshIterator it = subchannel_face->getIterator(AMP::Mesh::Face,1);
        std::vector<size_t> dofs;
        for (size_t i=0; i<it.size(); i++) {
            subchannel_DOFs->getDOFs(it->globalID(),dofs);
            AMP_ASSERT(dofs.size()==1);
            std::vector<double> pos = it->centroid();
            double v1 = T2->getValueByGlobalID(dofs[0]);
            double v2 = getTemp(pos);
            if ( !AMP::Utilities::approx_equal(v1,v2) )
                passes = false;
        }
        if ( passes )
            ut->passes("correctly mapped temperature");
        else
            ut->failure("correctly mapped temperature");
    }

    // Write the results
    #ifdef USE_EXT_SILO
        AMP::Mesh::SiloIO::shared_ptr  siloWriter( new AMP::Mesh::SiloIO);
        if ( T1.get()!=NULL )
            siloWriter->registerVector( T1, pin_mesh, AMP::Mesh::Vertex, "Temperature" );
        if ( T2.get()!=NULL )
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
    runTest ( "inputSubchannelToCladMap-1" , &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;

}
