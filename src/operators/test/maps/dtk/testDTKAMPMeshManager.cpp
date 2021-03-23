#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/operators/map/dtk/DTKAMPMeshEntityIterator.h"
#include "AMP/operators/map/dtk/DTKAMPMeshManager.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/VectorBuilder.h"
#include <memory>

#include <DTK_BasicEntityPredicates.hpp>

#include <cstdlib>
#include <iostream>
#include <string>


bool selectAll( DataTransferKit::Entity entity ) { return true; }

static void myTest( AMP::UnitTest *ut )
{
    std::string exeName( "testDTKAMPMeshManager" );
    std::string log_file = "output_" + exeName;
    std::string msgPrefix;
    AMP::logOnlyNodeZero( log_file );

    AMP::pout << "Loading the  mesh" << std::endl;
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    std::string input_file = "input_" + exeName;
    auto input_db          = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto meshDatabase = input_db->getDatabase( "Mesh" );

    auto meshParams = std::make_shared<AMP::Mesh::MeshParameters>( meshDatabase );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    auto mesh = AMP::Mesh::Mesh::buildMesh( meshParams );

    int const ghostWidth  = 0;
    int const dofsPerNode = 1;
    auto dofManager       = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, ghostWidth, dofsPerNode );

    // Create a mesh manager.
    DataTransferKit::SelectAllPredicate predicate;
    AMP::Operator::DTKAMPMeshManager dtk_mesh_manager( mesh, dofManager, predicate.getFunction() );

    // Get the function space.
    auto function_space = dtk_mesh_manager.functionSpace();

    // Test the entity set and entity selector by getting an iterator over the nodes.
    auto node_iterator = mesh->getIterator( AMP::Mesh::GeomType::Vertex, ghostWidth );
    auto dtk_node_iterator =
        function_space->entitySet()->entityIterator( 0, function_space->selectFunction() );
    AMP_ASSERT( dtk_node_iterator.size() == node_iterator.size() );
    AMP_ASSERT( dtk_node_iterator.size() > 0 );

    // Test the function space shape function and local map.
    std::vector<std::size_t> node_dofs;
    Teuchos::Array<std::size_t> dof_ids;
    Teuchos::Array<double> centroid( 3 );
    std::vector<double> coords;
    for ( dtk_node_iterator = dtk_node_iterator.begin(), node_iterator = node_iterator.begin();
          dtk_node_iterator != dtk_node_iterator.end();
          ++dtk_node_iterator, ++node_iterator ) {
        // Local map.
        function_space->localMap()->centroid( *dtk_node_iterator, centroid() );
        coords = node_iterator->coord();
        AMP_ASSERT( centroid[0] == coords[0] );
        AMP_ASSERT( centroid[1] == coords[1] );
        AMP_ASSERT( centroid[2] == coords[2] );

        // Shape function.
        function_space->shapeFunction()->entitySupportIds( *dtk_node_iterator, dof_ids );
        dofManager->getDOFs( node_iterator->globalID(), node_dofs );
        AMP_ASSERT( 1 == dof_ids.size() );
        AMP_ASSERT( 1 == node_dofs.size() );
        AMP_ASSERT( dof_ids[0] == node_dofs[0] );
    }

    ut->passes( exeName );
}


int main( int argc, char *argv[] )
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );
    AMP::UnitTest ut;

    myTest( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
