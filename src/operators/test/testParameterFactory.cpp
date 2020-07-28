#include "AMP/utils/AMPManager.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include <memory>

#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"

#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/vectors/VectorBuilder.h"

#include "AMP/operators/ParameterFactory.h"
#include "AMP/operators/boundary/DirichletMatrixCorrectionParameters.h"


static void ParameterFactoryTest( AMP::UnitTest *ut )
{
    std::string exeName( "testParameterFactory-1" );
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );


    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Get the Mesh database and create the mesh parameters
    std::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
    std::shared_ptr<AMP::Mesh::MeshParameters> params( new AMP::Mesh::MeshParameters( database ) );
    params->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );

    // Create the meshes from the input database
    AMP::Mesh::Mesh::shared_ptr mesh = AMP::Mesh::Mesh::buildMesh( params );

    AMP_INSIST( input_db->keyExists( "Parameter" ), "Key ''Parameter'' is missing!" );
    std::shared_ptr<AMP::Database> elemOp_db = input_db->getDatabase( "Parameter" );
    std::shared_ptr<AMP::Operator::OperatorParameters> operatorParameters =
        AMP::Operator::ParameterFactory::createParameter( elemOp_db, mesh );

    if ( elemOp_db->getString( "name" ) == "DirichletMatrixCorrection" ) {
        std::shared_ptr<AMP::Operator::DirichletMatrixCorrectionParameters> operatorParams =
            std::dynamic_pointer_cast<AMP::Operator::DirichletMatrixCorrectionParameters>(
                operatorParameters );

        if ( operatorParams.get() != nullptr )
            ut->passes( exeName.c_str() );
        else
            ut->failure( exeName.c_str() );
    } else {
        ut->failure( exeName.c_str() );
    }
}


int testParameterFactory( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    ParameterFactoryTest( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
