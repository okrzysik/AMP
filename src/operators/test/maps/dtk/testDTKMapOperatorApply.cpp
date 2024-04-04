#include "AMP/IO/PIO.h"
#include "AMP/IO/Writer.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshFactory.h"
#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/map/dtk/DTKMapOperator.h"
#include "AMP/operators/map/dtk/MultiDofDTKMapOperator.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/VectorSelector.h"
#include <discretization/simpleDOF_Manager.h>


static void thermalTest( AMP::UnitTest *ut, const std::string &input_file )
{
    std::string log_file = "log_DTKMapOperatorApply";
    std::string out_file = "out_DTKMapOperatorApply";

    AMP::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto mesh_db   = input_db->getDatabase( "Mesh" );
    auto mgrParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    auto manager = AMP::Mesh::MeshFactory::create( mgrParams );
    AMP::pout << "Finished loading meshes" << std::endl;

    int DOFsPerNode     = 1;
    int nodalGhostWidth = 1;
    bool split          = true;
    auto nodalDofMap    = AMP::Discretization::simpleDOFManager::create(
        manager, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );

    auto thermalVariable = std::make_shared<AMP::LinearAlgebra::Variable>( "Temperature" );

    auto SolutionVec      = AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable );
    auto RightHandSideVec = AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable );
    auto ResidualVec      = AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable );
    auto RankVec          = AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable );

    auto thermalMapVec = AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable, true );

    auto siloWriter = AMP::IO::Writer::buildWriter( "Silo" );
    siloWriter->registerMesh( manager );
    siloWriter->registerVector( SolutionVec, manager, AMP::Mesh::GeomType::Vertex, "SolutionVec" );
    siloWriter->registerVector( thermalMapVec, manager, AMP::Mesh::GeomType::Vertex, "MapVec" );
    siloWriter->registerVector( RankVec, manager, AMP::Mesh::GeomType::Vertex, "RankVec" );

    RightHandSideVec->setToScalar( 0.0 );

    double initialMapValue = input_db->getWithDefault<double>( "IMapValue", 1. );
    SolutionVec->setToScalar( 298. );
    thermalMapVec->setToScalar( initialMapValue );
    RankVec->setToScalar( globalComm.getRank() );
    RightHandSideVec->setToScalar( 0 );

    siloWriter->writeFile( out_file, 0 );

    auto DTKdb        = input_db->getDatabase( "DTKMaps" );
    auto mapColParams = std::make_shared<AMP::Operator::ColumnOperatorParameters>( input_db );
    auto mapsColumn   = std::make_shared<AMP::Operator::ColumnOperator>( mapColParams );

    auto nmaps          = DTKdb->getVector<int>( "N_maps" );
    auto mesh1          = DTKdb->getVector<std::string>( "Mesh1" );
    auto mesh2          = DTKdb->getVector<std::string>( "Mesh2" );
    auto surfaceID1     = DTKdb->getVector<int>( "Surface1" );
    auto surfaceID2     = DTKdb->getVector<int>( "Surface2" );
    auto var1           = DTKdb->getVector<std::string>( "Variable1" );
    auto var2           = DTKdb->getVector<std::string>( "Variable2" );
    auto inputDofsSize  = DTKdb->getVector<int>( "InputDOFsPerObject" );
    auto inputStride    = DTKdb->getVector<int>( "InputStride" );
    auto outputDofsSize = DTKdb->getVector<int>( "OutputDOFsPerObject" );
    auto outputStride   = DTKdb->getVector<int>( "OutputStride" );

    AMP::pout << "----------------------------\n";
    AMP::pout << "     CREATE MAP OPERATOR    \n";
    AMP::pout << "----------------------------\n";
    std::shared_ptr<AMP::Database> nullDatabase;
    for ( size_t i = 0; i < nmaps.size(); i++ ) {
        AMP::pout << "interface b/w" << std::endl;
        auto mapOperatorParams =
            std::make_shared<AMP::Operator::MultiDofDTKMapOperatorParameters>( nullDatabase );
        mapOperatorParams->d_globalComm    = AMP_COMM_WORLD;
        mapOperatorParams->d_Mesh1         = manager->Subset( mesh1[i] );
        mapOperatorParams->d_BoundaryID1   = surfaceID1[i];
        mapOperatorParams->d_Variable1     = var1[i];
        mapOperatorParams->d_StrideOffset1 = inputStride[i];
        mapOperatorParams->d_StrideLength1 = inputDofsSize[i];
        mapOperatorParams->d_Mesh2         = manager->Subset( mesh2[i] );
        mapOperatorParams->d_BoundaryID2   = surfaceID2[i];
        mapOperatorParams->d_Variable2     = var2[i];
        mapOperatorParams->d_StrideOffset2 = outputStride[i];
        mapOperatorParams->d_StrideLength2 = outputDofsSize[i];
        mapOperatorParams->d_SourceVector  = SolutionVec;
        mapOperatorParams->d_TargetVector  = thermalMapVec;
        auto mapOperator =
            std::make_shared<AMP::Operator::MultiDofDTKMapOperator>( mapOperatorParams );
        mapsColumn->append( mapOperator );
    }

    siloWriter->writeFile( out_file, 1 );
    mapsColumn->apply( SolutionVec, ResidualVec );
    siloWriter->writeFile( out_file, 2 );

    AMP::pout << " L2Norm of Map Vec " << std::setprecision( 17 ) << thermalMapVec->L2Norm()
              << std::endl;
    AMP::pout << " Thermal Map Vec Max : " << thermalMapVec->max() << " Min "
              << thermalMapVec->min() << std::endl;

    if ( AMP::Utilities::approx_equal( thermalMapVec->max(), 298. ) ) {
        ut->passes( "Mapped Values are consistent" );
    } else {
        ut->failure( "Mapped Values are inconsistent" );
    }

    ut->passes( input_file );
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::string inputFile = "input_testDTKMapOperatorApply";

    try {
        if ( argc > 1 )
            inputFile = argv[1];
        thermalTest( &ut, inputFile );
    } catch ( std::exception &err ) {
        AMP::pout << "ERROR: While testing " << argv[0] << err.what() << std::endl;
        ut.failure( "ERROR: While testing" );
    } catch ( ... ) {
        AMP::pout << "ERROR: While testing " << argv[0] << "An unknown exception was thrown."
                  << std::endl;
        ut.failure( "ERROR: While testing" );
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();

    AMP::AMPManager::shutdown();

    return num_failed;
}
