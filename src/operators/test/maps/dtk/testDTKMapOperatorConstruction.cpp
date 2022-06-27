#include "AMP/IO/PIO.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshFactory.h"
#include "AMP/operators/map/dtk/DTKMapOperator.h"
#include "AMP/operators/map/dtk/MultiDofDTKMapOperator.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/VectorSelector.h"
#include <discretization/simpleDOF_Manager.h>


static void dtkConsruction( AMP::UnitTest *ut, std::string input_file )
{
    std::string log_file = "log_DTK";
    AMP::logOnlyNodeZero( log_file );


    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    // PROFILE_START("SetupDriver");
    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto mesh_db   = input_db->getDatabase( "Mesh" );
    auto mgrParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    auto manager = std::make_shared<AMP::Mesh::MeshFactory::create>( mgrParams );
    AMP::pout << "Finished loading meshes" << std::endl;

    auto nodalDofMap = AMP::Discretization::simpleDOFManager::create(
        manager, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    auto variables = std::make_shared<AMP::LinearAlgebra::Variable>( "Field" );
    auto fieldVec  = AMP::LinearAlgebra::createVector( nodalDofMap, variable );
    auto resVec    = AMP::LinearAlgebra::createVector( nodalDofMap, variable );
    auto dummyVec  = AMP::LinearAlgebra::createVector( nodalDofMap, variable );
    auto mapVec    = AMP::LinearAlgebra::createVector( nodalDofMap, variable, true );

    AMP::pout << "----------------------------\n";
    AMP::pout << "     CREATE MAP OPERATOR    \n";
    AMP::pout << "----------------------------\n";
    std::shared_ptr<AMP::Database> nullDatabase;

    auto Mesh1 = manager->Subset( "Mesh1" );
    auto Mesh2 = manager->Subset( "Mesh2" );

    auto mapOperatorParams =
        std::make_shared<AMP::Operator::MultiDofDTKMapOperatorParameters>( nullDatabase );
    mapOperatorParams->d_globalComm    = AMP_COMM_WORLD;
    mapOperatorParams->d_Mesh1         = Mesh1;
    mapOperatorParams->d_BoundaryID1   = 1;
    mapOperatorParams->d_Variable1     = variable->getName();
    mapOperatorParams->d_StrideOffset1 = 0;
    mapOperatorParams->d_StrideLength1 = 1;
    mapOperatorParams->d_Mesh2         = Mesh2;
    mapOperatorParams->d_BoundaryID2   = 3;
    mapOperatorParams->d_Variable2     = variable->getName();
    mapOperatorParams->d_StrideOffset2 = 0;
    mapOperatorParams->d_StrideLength2 = 1;
    mapOperatorParams->d_SourceVector  = fieldVec;
    mapOperatorParams->d_TargetVector  = mapVec;
    auto mapOperator = std::make_shared<AMP::Operator::MultiDofDTKMapOperator>( mapOperatorParams );
    ut->passes( "DTK Map Operator creation" );

    try {
        mapOperator->apply( fieldVec, resVec );
        ut->passes( "DTK Map Operator apply" );
    } catch ( std::exception ) {
        ut->failure( "DTK Map Operator apply" );
    }
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::string inputFile = "input_testDTKConstruction";
    if ( argc > 1 )
        inputFile = argv[1];
    dtkConsruction( &ut, inputFile );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
