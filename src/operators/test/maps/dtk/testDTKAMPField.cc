
#include <utils/AMPManager.h>
#include <utils/AMP_MPI.h>
#include <utils/Database.h>
#include <utils/InputDatabase.h>
#include <utils/InputManager.h>
#include <utils/PIO.h>
#include <utils/UnitTest.h>
#include <utils/Utilities.h>
#include <utils/shared_ptr.h>

#include <vectors/VectorBuilder.h>

#include <discretization/simpleDOF_Manager.h>

#include <ampmesh/Mesh.h>

#include <operators/map/dtk/DTKAMPField.h>

#include <cstdlib>
#include <iostream>
#include <string>


void myTest( AMP::UnitTest *ut )
{
    std::string exeName( "testDTKAMPMeshManager" );
    std::string log_file = "output_" + exeName;
    std::string msgPrefix;
    AMP::PIO::logOnlyNodeZero( log_file );

    // Load input and build the mesh.
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

    // Create a vector.
    bool const split      = true;
    int const ghostWidth  = 0;
    int const dofsPerNode = 1;
    AMP::Discretization::DOFManager::shared_ptr dofManager =
        AMP::Discretization::simpleDOFManager::create(
            mesh, AMP::Mesh::Vertex, ghostWidth, dofsPerNode );
    AMP::LinearAlgebra::Variable::shared_ptr variable( new AMP::LinearAlgebra::Variable( "var" ) );
    AMP::LinearAlgebra::Vector::shared_ptr ampVector =
        AMP::LinearAlgebra::createVector( dofManager, variable, split );
    std::vector<std::size_t> dofIndices;
    AMP::Mesh::MeshIterator meshIterator = mesh->getIterator( AMP::Mesh::Vertex, ghostWidth );
    for ( meshIterator = meshIterator.begin(); meshIterator != meshIterator.end();
          ++meshIterator ) {
        dofManager->getDOFs( meshIterator->globalID(), dofIndices );
        ampVector->setLocalValueByGlobalID( dofIndices[0], static_cast<double>( dofIndices[0] ) );
    }

    // Create a dtk field around the amp vector
    AMP::Operator::DTKAMPField dtk_field( ampVector );

    // Check the dimension.
    AMP_ASSERT( 1 == dtk_field.dimension() );

    // Check the support ids.
    Teuchos::ArrayView<const DataTransferKit::SupportId> support_ids =
        dtk_field.getLocalSupportIds();
    int counter = 0;
    for ( meshIterator = meshIterator.begin(); meshIterator != meshIterator.end();
          ++meshIterator, ++counter ) {
        dofManager->getDOFs( meshIterator->globalID(), dofIndices );
        AMP_ASSERT( support_ids[counter] == dofIndices[0] );
    }

    // Check reading data.
    for ( meshIterator = meshIterator.begin(); meshIterator != meshIterator.end();
          ++meshIterator, ++counter ) {
        dofManager->getDOFs( meshIterator->globalID(), dofIndices );
        AMP_ASSERT( dtk_field.readFieldData( dofIndices[0], 0 ) ==
                    ampVector->getLocalValueByGlobalID( dofIndices[0] ) );
    }

    // Check setting data.
    for ( meshIterator = meshIterator.begin(); meshIterator != meshIterator.end();
          ++meshIterator, ++counter ) {
        dofManager->getDOFs( meshIterator->globalID(), dofIndices );
        dtk_field.writeFieldData( dofIndices[0], 0, 2.0 * dofIndices[0] );
    }
    for ( meshIterator = meshIterator.begin(); meshIterator != meshIterator.end();
          ++meshIterator, ++counter ) {
        dofManager->getDOFs( meshIterator->globalID(), dofIndices );
        AMP_ASSERT( 2.0 * dofIndices[0] == ampVector->getLocalValueByGlobalID( dofIndices[0] ) );
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
