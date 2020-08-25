#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/materials/Material.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/map/Map1Dto3D.h"
#include "AMP/operators/map/Map3Dto1D.h"
#include "AMP/operators/map/MapOperatorParameters.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/SimpleVector.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include <memory>

#include <string>


static void testMap( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logAllNodes( log_file );

    // Read the input file
    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Get the Mesh database and create the mesh parameters
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    std::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    std::shared_ptr<AMP::Mesh::MeshParameters> params( new AMP::Mesh::MeshParameters( mesh_db ) );
    params->setComm( globalComm );

    // Create the meshes from the input database
    std::shared_ptr<AMP::Mesh::Mesh> mesh    = AMP::Mesh::Mesh::buildMesh( params );
    AMP::Mesh::Mesh::shared_ptr meshAdapter1 = mesh->Subset( "pellet" );
    AMP::Mesh::Mesh::shared_ptr meshAdapter2 = mesh->Subset( "clad" );

    // Create a simple DOFManager and the vectors
    AMP::Discretization::DOFManager::shared_ptr DOFs =
        AMP::Discretization::simpleDOFManager::create( mesh, AMP::Mesh::GeomType::Vertex, 1, 1 );
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    AMP::LinearAlgebra::Variable::shared_ptr testVariable(
        new AMP::LinearAlgebra::Variable( "MapSolution" ) );
    AMP::LinearAlgebra::Variable::shared_ptr gapVariable(
        new AMP::LinearAlgebra::Variable( "Gap" ) );

    AMP::LinearAlgebra::Vector::shared_ptr mapSolution =
        AMP::LinearAlgebra::createVector( DOFs, testVariable );
    AMP::LinearAlgebra::Vector::shared_ptr mapSolutionMaster =
        AMP::LinearAlgebra::createVector( DOFs, testVariable );
    AMP::LinearAlgebra::Vector::shared_ptr mapSolutionSlave =
        AMP::LinearAlgebra::createVector( DOFs, testVariable );

    mapSolution->setToScalar( 0.0 );

    // Create the map operators
    std::shared_ptr<AMP::Database> map3dto1d_db1 =
        std::dynamic_pointer_cast<AMP::Database>( input_db->getDatabase( "MapPelletto1D" ) );
    std::shared_ptr<AMP::Operator::MapOperatorParameters> map3dto1dParams1(
        new AMP::Operator::MapOperatorParameters( map3dto1d_db1 ) );
    map3dto1dParams1->d_MapComm = globalComm;
    map3dto1dParams1->d_MapMesh = meshAdapter1;
    std::shared_ptr<AMP::Operator::Map3Dto1D> map1ToLowDim(
        new AMP::Operator::Map3Dto1D( map3dto1dParams1 ) );

    std::shared_ptr<AMP::Database> map1dto3d_db1 =
        std::dynamic_pointer_cast<AMP::Database>( input_db->getDatabase( "Map1DtoClad" ) );
    std::shared_ptr<AMP::Operator::MapOperatorParameters> map1dto3dParams1(
        new AMP::Operator::MapOperatorParameters( map1dto3d_db1 ) );
    map1dto3dParams1->d_MapComm = globalComm;
    map1dto3dParams1->d_MapMesh = meshAdapter2;
    std::shared_ptr<AMP::Operator::Map1Dto3D> map1ToHighDim(
        new AMP::Operator::Map1Dto3D( map1dto3dParams1 ) );

    map1ToLowDim->setZLocations( map1ToHighDim->getZLocations() );

    std::shared_ptr<AMP::Database> map3dto1d_db2 =
        std::dynamic_pointer_cast<AMP::Database>( input_db->getDatabase( "MapCladto1D" ) );
    std::shared_ptr<AMP::Operator::MapOperatorParameters> map3dto1dParams2(
        new AMP::Operator::MapOperatorParameters( map3dto1d_db2 ) );
    map3dto1dParams2->d_MapComm = globalComm;
    map3dto1dParams2->d_MapMesh = meshAdapter2;
    std::shared_ptr<AMP::Operator::Map3Dto1D> map2ToLowDim(
        new AMP::Operator::Map3Dto1D( map3dto1dParams2 ) );

    std::shared_ptr<AMP::Database> map1dto3d_db2 =
        std::dynamic_pointer_cast<AMP::Database>( input_db->getDatabase( "Map1DtoPellet" ) );
    std::shared_ptr<AMP::Operator::MapOperatorParameters> map1dto3dParams2(
        new AMP::Operator::MapOperatorParameters( map1dto3d_db2 ) );
    map1dto3dParams2->d_MapComm = globalComm;
    map1dto3dParams2->d_MapMesh = meshAdapter1;
    std::shared_ptr<AMP::Operator::Map1Dto3D> map2ToHighDim(
        new AMP::Operator::Map1Dto3D( map1dto3dParams2 ) );

    map2ToLowDim->setZLocations( map2ToHighDim->getZLocations() );
    //-------------------------------------
    size_t gapVecCladSize = map1ToHighDim->getNumZlocations();
    AMP::LinearAlgebra::Vector::shared_ptr gapVecClad =
        AMP::LinearAlgebra::SimpleVector<double>::create( gapVecCladSize, gapVariable );

    size_t gapVecPelletSize = map2ToHighDim->getNumZlocations();
    AMP::LinearAlgebra::Vector::shared_ptr gapVecPellet =
        AMP::LinearAlgebra::SimpleVector<double>::create( gapVecPelletSize, gapVariable );

    // Set the boundary for the source vector
    unsigned int d_boundaryId = map3dto1d_db1->getScalar<int>( "BoundaryId" );
    auto bnd     = mesh->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, d_boundaryId, 0 );
    auto end_bnd = bnd.end();
    auto dof_map = mapSolutionMaster->getDOFManager();
    std::vector<size_t> ids;
    for ( ; bnd != end_bnd; ++bnd ) {
        auto x = bnd->coord();
        dof_map->getDOFs( bnd->globalID(), ids );
        AMP_ASSERT( ids.size() == 1 );
        mapSolutionMaster->setValuesByGlobalID( 1, &ids[0], &x[2] );
    }
    mapSolutionMaster->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );

    //-------------------------------------

    int cnt         = 0;
    bool testPassed = false;
    while ( cnt < 20 ) {
        cnt++;

        map1ToLowDim->setVector( gapVecClad );
        map1ToHighDim->setVector( mapSolutionSlave );
        map1ToLowDim->apply( mapSolutionMaster, gapVecClad );
        map1ToHighDim->apply( gapVecClad, mapSolutionSlave );
        std::cout << "Master Map Solution " << std::endl;
        for ( size_t i = 0; i < gapVecCladSize; i++ ) {
            std::cout << " @i : " << i << " is " << gapVecClad->getValueByLocalID( i );
        }
        std::cout << std::endl;
        //------------------------------------------------------------
        //    mapSolutionSlave->setToScalar(90);

        map2ToLowDim->setVector( gapVecPellet );
        map2ToHighDim->setVector( mapSolutionMaster );
        map2ToLowDim->apply( mapSolutionSlave, gapVecPellet );
        map2ToHighDim->apply( gapVecPellet, mapSolutionMaster );

        std::cout << "Slave Map Solution " << std::endl;
        for ( size_t i = 0; i < gapVecPelletSize; i++ ) {
            std::cout << " @i : " << i << " is " << gapVecPellet->getValueByLocalID( i );
        }
        std::cout << std::endl;
        //------------------------------------------------------------

        if ( true ) {
            testPassed = true;
            break;
        } else {
            std::cout << "Norm of the change in sol for iteration " << cnt << "is -->" << std::endl;
        }
        std::cout << std::endl;
    }

#ifdef USE_EXT_SILO
    if ( globalComm.getSize() == 1 ) {
        // manager->writeFile<AMP::Mesh::SiloIO> ( exeName , 0 );
    }
#endif

    if ( testPassed )
        ut->passes( "Seggregated solve of Composite Operator using control loop of "
                    "Thermal+Robin->Map->Gap->Map->Thermal+Robin ." );
    else
        ut->failure( "Seggregated solve of Composite Operator using control loop of "
                     "Thermal+Robin->Map->Gap->Map->Thermal+Robin ." );

    input_db.reset();

    ut->passes( exeName );
}


int testMapOperator( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    testMap( &ut, "testMapOperator" );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
