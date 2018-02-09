#include "AMP/ampmesh/Mesh.h"
#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/map/dtk/DTKMapOperator.h"
#include "AMP/operators/map/dtk/MultiDofDTKMapOperator.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/InputDatabase.h"
#include "AMP/utils/InputManager.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/VectorSelector.h"
#include <discretization/simpleDOF_Manager.h>


void thermalTest( AMP::UnitTest *ut, const std::string &input_file )
{
    std::string log_file = "log_DTKMapOperatorApply";
    std::string out_file = "out_DTKMapOperatorApply";

    AMP::PIO::logOnlyNodeZero( log_file );

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    // PROFILE_START("SetupDriver");
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::Mesh::Mesh::shared_ptr manager( AMP::Mesh::Mesh::buildMesh( mgrParams ) );
    AMP::pout << "Finished loading meshes" << std::endl;

    int DOFsPerNode     = 1;
    int nodalGhostWidth = 1;
    bool split          = true;
    AMP::Discretization::DOFManager::shared_ptr nodalDofMap =
        AMP::Discretization::simpleDOFManager::create(
            manager, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );

    AMP::LinearAlgebra::Variable::shared_ptr thermalVariable(
        new AMP::LinearAlgebra::Variable( "Temperature" ) );

    AMP::LinearAlgebra::Vector::shared_ptr SolutionVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable );
    AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable );
    AMP::LinearAlgebra::Vector::shared_ptr ResidualVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable );
    AMP::LinearAlgebra::Vector::shared_ptr RankVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable );

    AMP::LinearAlgebra::Vector::shared_ptr thermalMapVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable, true );

    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    siloWriter->registerMesh( manager );
    siloWriter->registerVector( SolutionVec, manager, AMP::Mesh::GeomType::Vertex, "SolutionVec" );
    siloWriter->registerVector( thermalMapVec, manager, AMP::Mesh::GeomType::Vertex, "MapVec" );
    siloWriter->registerVector( RankVec, manager, AMP::Mesh::GeomType::Vertex, "RankVec" );

    RightHandSideVec->setToScalar( 0.0 );


    AMP::Mesh::Mesh::shared_ptr cellSandwichMesh = manager->Subset( "CellSandwich" );
    AMP::Mesh::Mesh::shared_ptr ccMesh           = manager->Subset( "CellCurrentCollectors" );

    double initialMapValue = input_db->getDoubleWithDefault( "IMapValue", 1. );
    SolutionVec->setToScalar( 298. );
    thermalMapVec->setToScalar( initialMapValue );
    RankVec->setToScalar( globalComm.getRank() );
    RightHandSideVec->setToScalar( 0 );

    siloWriter->writeFile( out_file, 0 );

    AMP::shared_ptr<AMP::Database> DTKdb = input_db->getDatabase( "DTKMaps" );
    AMP::shared_ptr<AMP::Operator::ColumnOperatorParameters> mapColParams(
        new AMP::Operator::ColumnOperatorParameters( input_db ) );
    AMP::shared_ptr<AMP::Operator::ColumnOperator> mapsColumn(
        new AMP::Operator::ColumnOperator( mapColParams ) );

    std::vector<int> nmaps          = DTKdb->getIntegerArray( "N_maps" );
    std::vector<std::string> mesh1  = DTKdb->getStringArray( "Mesh1" );
    std::vector<std::string> mesh2  = DTKdb->getStringArray( "Mesh2" );
    std::vector<int> surfaceID1     = DTKdb->getIntegerArray( "Surface1" );
    std::vector<int> surfaceID2     = DTKdb->getIntegerArray( "Surface2" );
    std::vector<std::string> var1   = DTKdb->getStringArray( "Variable1" );
    std::vector<std::string> var2   = DTKdb->getStringArray( "Variable2" );
    std::vector<int> inputDofsSize  = DTKdb->getIntegerArray( "InputDOFsPerObject" );
    std::vector<int> inputStride    = DTKdb->getIntegerArray( "InputStride" );
    std::vector<int> outputDofsSize = DTKdb->getIntegerArray( "OutputDOFsPerObject" );
    std::vector<int> outputStride   = DTKdb->getIntegerArray( "OutputStride" );

    AMP::pout << "----------------------------\n";
    AMP::pout << "     CREATE MAP OPERATOR    \n";
    AMP::pout << "----------------------------\n";
    AMP::shared_ptr<AMP::Database> nullDatabase;
    for ( size_t i = 0; i < nmaps.size(); i++ ) {
        AMP::pout << "interface b/w" << std::endl;
        AMP::shared_ptr<AMP::Operator::MultiDofDTKMapOperatorParameters> mapOperatorParams(
            new AMP::Operator::MultiDofDTKMapOperatorParameters( nullDatabase ) );
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
        AMP::shared_ptr<AMP::Operator::Operator> mapOperator(
            new AMP::Operator::MultiDofDTKMapOperator( mapOperatorParams ) );

        mapsColumn->append( mapOperator );
    }

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

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

    // PROFILE_SAVE(inputFile);

    AMP::AMPManager::shutdown();

    return num_failed;
}
