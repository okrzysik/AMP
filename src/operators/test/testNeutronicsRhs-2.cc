#include "materials/Material.h"
#include "utils/AMPManager.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/shared_ptr.h"
#include "vectors/Variable.h"
#include <string>

#include "operators/ElementOperationFactory.h"
#include "operators/ElementPhysicsModelFactory.h"
#include "operators/NeutronicsRhs.h"
#include "operators/libmesh/SourceNonlinearElement.h"
#include "utils/Writer.h"
#include "vectors/Vector.h"

#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/OperatorBuilder.h"

#include <exception>


void sourceTest( AMP::UnitTest *ut, const std::string &exeName )
{
    // Initialization
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logAllNodes( log_file );

    std::cout << "testing with input file " << input_file << std::endl;


    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh( mgrParams );

    //--------------------------------------------------
    //   CREATE THE VOLUME INTEGRAL OPERATOR -----------
    //--------------------------------------------------

    AMP_INSIST( input_db->keyExists( "NeutronicsRhs" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> unusedModel;

    // Construct with OperatorBuilder
    {
        AMP::shared_ptr<AMP::Operator::NeutronicsRhs> ntxBld =
            AMP::dynamic_pointer_cast<AMP::Operator::NeutronicsRhs>(
                AMP::Operator::OperatorBuilder::createOperator(
                    meshAdapter, "NeutronicsRhs", input_db, unusedModel ) );
        AMP_INSIST( ntxBld.get() != nullptr, "NULL rhs out of OperatorBuilder" );
        ut->passes( "NeutronicsRhs was constructed by OperatorBuilder for: " + input_file );
        // ntxBld->setTimeStep(0);
        // ut->passes( "NeutronicsRhs, constructed by OperatorBuilder, set the time for: " +
        // input_file);
    }
}


// Input and output file names
int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    const int NUMFILES          = 1;
    std::string files[NUMFILES] = { "testNeutronicsRhs-db" };

    for ( auto &file : files )
        sourceTest( &ut, file );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
