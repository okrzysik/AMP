#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MeshElementVectorIterator.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/operators/map/AsyncMapColumnOperator.h"
#include "AMP/operators/map/CladToSubchannelMap.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/VectorBuilder.h"


double getTemp( const AMP::Mesh::Point &x ) { return 500 + x[0] * 100 + x[1] * 100 + x[2] * 100; }


AMP::Mesh::MeshIterator getZFaceIterator( AMP::Mesh::Mesh::shared_ptr subChannel, int ghostWidth )
{
    std::multimap<double, AMP::Mesh::MeshElement> xyFace;
    AMP::Mesh::MeshIterator iterator =
        subChannel->getIterator( AMP::Mesh::GeomType::Face, ghostWidth );
    for ( size_t i = 0; i < iterator.size(); ++i ) {
        auto nodes    = iterator->getElements( AMP::Mesh::GeomType::Vertex );
        auto center   = iterator->centroid();
        bool is_valid = true;
        for ( auto &node : nodes ) {
            auto coord = node.coord();
            if ( !AMP::Utilities::approx_equal( coord[2], center[2], 1e-6 ) )
                is_valid = false;
        }
        if ( is_valid ) {
            xyFace.insert( std::pair<double, AMP::Mesh::MeshElement>( center[2], *iterator ) );
        }
        ++iterator;
    }
    std::shared_ptr<std::vector<AMP::Mesh::MeshElement>> elements(
        new std::vector<AMP::Mesh::MeshElement>() );
    elements->reserve( xyFace.size() );
    for ( auto &elem : xyFace )
        elements->push_back( elem.second );
    return AMP::Mesh::MultiVectorIterator( elements );
}


static void runTest( const std::string &fname, AMP::UnitTest *ut )
{
    // Read the input file
    auto input_db = AMP::Database::parseInputFile( fname );
    input_db->print( AMP::plog );

    // Get the Mesh database and create the mesh parameters
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    auto mesh_db = input_db->getDatabase( "Mesh" );
    auto params  = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    params->setComm( globalComm );

    // Create the meshes from the input database
    auto manager  = AMP::Mesh::Mesh::buildMesh( params );
    auto pin_mesh = manager->Subset( "MultiPin" );
    AMP::Mesh::Mesh::shared_ptr clad_mesh;
    if ( pin_mesh.get() != nullptr ) {
        pin_mesh->setName( "MultiPin" );
        clad_mesh = pin_mesh->Subset( "clad" );
    }
    auto subchannel_mesh = manager->Subset( "subchannel" );
    AMP::Mesh::Mesh::shared_ptr subchannel_face;
    if ( subchannel_mesh.get() != nullptr ) {
        subchannel_mesh->setName( "subchannel" );
        subchannel_face = subchannel_mesh->Subset( getZFaceIterator( subchannel_mesh, 1 ) );
    }

    // Get the database for the map
    auto map_db = input_db->getDatabase( "MeshToMeshMaps" );

    // Create the DOFManagers and the vectors
    // int DOFsPerNode = map_db->getScalar<int>("DOFsPerObject");
    // std::string varName = map_db->getString("VariableName");
    int DOFsPerNode     = 1;
    std::string varName = "Temperature";
    auto temperature    = std::make_shared<AMP::LinearAlgebra::Variable>( varName );
    AMP::Discretization::DOFManager::shared_ptr pin_DOFs;
    AMP::Discretization::DOFManager::shared_ptr subchannel_DOFs;
    AMP::LinearAlgebra::Vector::shared_ptr T1;
    AMP::LinearAlgebra::Vector::shared_ptr T2;
    AMP::LinearAlgebra::Vector::shared_ptr dummy;
    if ( pin_mesh.get() != nullptr ) {
        pin_DOFs = AMP::Discretization::simpleDOFManager::create(
            pin_mesh, AMP::Mesh::GeomType::Vertex, 1, DOFsPerNode );
        T1 = AMP::LinearAlgebra::createVector( pin_DOFs, temperature );
        T1->setToScalar( 0.0, T1 );
    }
    if ( subchannel_face.get() != nullptr ) {
        subchannel_DOFs = AMP::Discretization::simpleDOFManager::create(
            subchannel_face, AMP::Mesh::GeomType::Face, 1, DOFsPerNode );
        T2 = AMP::LinearAlgebra::createVector( subchannel_DOFs, temperature );
        T2->setToScalar( 0.0, T2 );
    }

    // Initialize the pin temperatures
    if ( pin_mesh.get() != nullptr ) {
        auto it = pin_mesh->getIterator( AMP::Mesh::GeomType::Vertex, 0 );
        std::vector<size_t> dofs;
        for ( size_t i = 0; i < it.size(); i++ ) {
            pin_DOFs->getDOFs( it->globalID(), dofs );
            T1->setValueByGlobalID( dofs[0], getTemp( it->coord() ) );
            ++it;
        }
    }

    // Test the creation/destruction of CladToSubchannelMap (no apply call)
    try {
        auto map = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::CladToSubchannelMap>(
            manager, map_db );
        map.reset();
        ut->passes( "Created / Destroyed CladToSubchannelMap" );
    } catch ( ... ) {
        ut->failure( "Created / Destroyed CladToSubchannelMap" );
    }

    // Perform a complete test of CladToSubchannelMap
    auto map = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::CladToSubchannelMap>(
        manager, map_db );
    map->setVector( T2 );

    // Apply the map
    globalComm.barrier();
    map->apply( T1, T2 );

// Check the results
/*if ( subchannel_face.get()!=NULL ) {
    bool passes = true;
    AMP::Mesh::MeshIterator it = subchannel_face->getIterator(AMP::Mesh::GeomType::Face,1);
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
}*/

// Write the results
#ifdef USE_EXT_SILO
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    if ( T1.get() != nullptr )
        siloWriter->registerVector( T1, pin_mesh, AMP::Mesh::GeomType::Vertex, "Temperature" );
    if ( T2.get() != nullptr )
        siloWriter->registerVector( T2, subchannel_face, AMP::Mesh::GeomType::Face, "Temperature" );
    siloWriter->setDecomposition( 1 );
    siloWriter->writeFile( fname, 0 );
#endif
}


int testCladToSubchannelMap( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    runTest( "inputCladToSubchannelMap-1", &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
