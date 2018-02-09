#include "AMP/materials/Material.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/InputDatabase.h"
#include "AMP/utils/InputManager.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/shared_ptr.h"
#include "AMP/vectors/Variable.h"
#include <string>

#include "AMP/operators/ElementOperationFactory.h"
#include "AMP/operators/ElementPhysicsModelFactory.h"
#include "AMP/operators/NeutronicsRhs.h"
#include "AMP/operators/libmesh/SourceNonlinearElement.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/Vector.h"

#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"

#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/vectors/VectorBuilder.h"

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
    AMP::shared_ptr<AMP::Database> ntx_db = input_db->getDatabase( "NeutronicsRhs" );

    // Construct stand-alone.
    {
        // construct it.
        AMP::shared_ptr<AMP::Operator::NeutronicsRhsParameters> ntxPrm(
            new AMP::Operator::NeutronicsRhsParameters( ntx_db ) );
        AMP::shared_ptr<AMP::Operator::NeutronicsRhs> ntxRhs(
            new AMP::Operator::NeutronicsRhs( ntxPrm ) );
        ut->passes( "NeutronicsRhs was constructed stand-alone for: " + input_file );
        // set the time.
        ntxRhs->setTimeStep( 0 );
        ut->passes( "NeutronicsRhs, constructed stand-alone, set the time for: " + input_file );
        // set the time.
        // ntxRhs->setTimeInSeconds(8000000.);
        // ut.passes( "NeutronicsRhs, constructed stand-alone, set the time in seconds for: " +
        // input_file);
        // Create a DOF manager for a gauss point vector
        int DOFsPerNode = 8;
        int ghostWidth  = 1;
        bool split      = true;
        AMP::Discretization::DOFManager::shared_ptr dof_map =
            AMP::Discretization::simpleDOFManager::create(
                meshAdapter, AMP::Mesh::GeomType::Volume, ghostWidth, DOFsPerNode, split );
        // create a variable/vector combo.
        AMP::LinearAlgebra::Vector::shared_ptr nullVec;
        // AMP::Operator::NeutronicsRhs::SP_HexGaussPointVariable outVar(new
        // AMP::Operator::NeutronicsRhs::HexGaussPointVariable("outpower") );
        AMP::LinearAlgebra::Variable::shared_ptr outVar = ntxRhs->getOutputVariable();
        AMP::LinearAlgebra::Vector::shared_ptr outVec =
            AMP::LinearAlgebra::createVector( dof_map, outVar, split );
        ntxRhs->apply( nullVec, outVec );
    }

    // Construct with OperatorBuilder
    {
        AMP::shared_ptr<AMP::Operator::NeutronicsRhs> ntxBld =
            AMP::dynamic_pointer_cast<AMP::Operator::NeutronicsRhs>(
                AMP::Operator::OperatorBuilder::createOperator(
                    meshAdapter, "NeutronicsRhs", input_db, unusedModel ) );
        AMP_INSIST( ntxBld.get() != nullptr, "NULL rhs out of OperatorBuilder" );
        ut->passes( "NeutronicsRhs was constructed by OperatorBuilder for: " + input_file );
        // ntxBld->setTimeStep(0);
        // ut.passes( "NeutronicsRhs, constructed by OperatorBuilder, set the time for: " +
        // input_file);
    }
}
// Input and output file names
int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    const int NUMFILES          = 4;
    std::string files[NUMFILES] = {
        "testNeutronicsRhs-1", "testNeutronicsRhs-2", "testNeutronicsRhs-3", "testNeutronicsRhs-4"
    };

    for ( auto &file : files )
        sourceTest( &ut, file );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
