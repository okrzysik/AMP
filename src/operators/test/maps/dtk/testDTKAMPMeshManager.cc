
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <cstdlib>
#include <iostream>
#include <string>

#include "utils/shared_ptr.h"

#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"

#include "vectors/VectorBuilder.h"

#include "discretization/simpleDOF_Manager.h"

#include "ampmesh/Mesh.h"

#include "operators/map/dtk/DTKAMPMeshEntityIterator.h"
#include "operators/map/dtk/DTKAMPMeshManager.h"

#include <DTK_BasicEntityPredicates.hpp>

bool selectAll( DataTransferKit::Entity entity ) { return true; }

void myTest( AMP::UnitTest *ut )
{
    std::string exeName( "testDTKAMPMeshManager" );
    std::string log_file = "output_" + exeName;
    std::string msgPrefix;
    AMP::PIO::logOnlyNodeZero( log_file );

    AMP::pout << "Loading the  mesh" << std::endl;
    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    std::string input_file = "input_" + exeName;
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    AMP::shared_ptr<AMP::Database> meshDatabase = input_db->getDatabase( "Mesh" );

    AMP::shared_ptr<AMP::Mesh::MeshParameters> meshParams(
        new AMP::Mesh::MeshParameters( meshDatabase ) );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::Mesh::Mesh::shared_ptr mesh = AMP::Mesh::Mesh::buildMesh( meshParams );

    bool const split      = true;
    int const ghostWidth  = 0;
    int const dofsPerNode = 1;
    AMP::Discretization::DOFManager::shared_ptr dofManager =
        AMP::Discretization::simpleDOFManager::create(
            mesh, AMP::Mesh::Vertex, ghostWidth, dofsPerNode );
    AMP::LinearAlgebra::Variable::shared_ptr variable( new AMP::LinearAlgebra::Variable( "var" ) );
    AMP::LinearAlgebra::Vector::shared_ptr vector =
        AMP::LinearAlgebra::createVector( dofManager, variable, split );

    // Create a mesh manager.
    DataTransferKit::SelectAllPredicate predicate;
    AMP::Operator::DTKAMPMeshManager dtk_mesh_manager(
        mesh, dofManager, predicate.getFunction() );

    // Get the function space.
    Teuchos::RCP<DataTransferKit::FunctionSpace> function_space = dtk_mesh_manager.functionSpace();

    // Test the entity set and entity selector by getting an iterator over the nodes.
    AMP::Mesh::MeshIterator node_iterator = mesh->getIterator( AMP::Mesh::Vertex, ghostWidth );
    DataTransferKit::EntityIterator dtk_node_iterator = 
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
