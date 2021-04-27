#include "AMP/operators/map/AsyncMapOperator.h"
#include "AMP/operators/map/AsyncMapOperatorParameters.h"
#include "AMP/utils/Database.h"


namespace AMP {
namespace Operator {


enum class MapDominance { Master, Slave };
enum class MapConstructionType { Synchronous, Asynchronous };
struct MapConstructionParam {
    int boundaryId;
    AMP_MPI comm;
    MapConstructionType construction;
    MapDominance dominance;
    std::shared_ptr<AMP::Database> database;
    std::string meshName;
    AMP::Mesh::Mesh::shared_ptr mesh;
    std::string mapType;
    size_t tagOffset;
};

extern size_t globalMapTagOffset; // We need a global unique tag offset for every map


template<typename MAP_TYPE>
std::shared_ptr<AsyncMapColumnOperator>
AsyncMapColumnOperator::build( AMP::Mesh::Mesh::shared_ptr manager,
                               std::shared_ptr<Database> database )
{

    auto newParams    = std::make_shared<AsyncMapColumnOperatorParameters>( database );
    auto newMapColumn = std::make_shared<AsyncMapColumnOperator>( newParams );

    // Check that the map type matches the maps we are building
    AMP_ASSERT( database->keyExists( "MapType" ) );
    AMP_ASSERT( MAP_TYPE::validMapType( database->getString( "MapType" ) ) );

    // Get the number of maps in the database
    AMP_INSIST( database->keyExists( "N_maps" ), "N_maps must exist in input database" );
    int N_maps = database->getScalar<int>( "N_maps" );

    // Create the databases for the individual maps
    auto map_databases = createDatabases( database );

    // Loop through the maps
    AMP_MPI managerComm = manager->getComm();
    for ( int i = 0; i < N_maps; i++ ) {
        // Get the names of the 2 meshes involved
        auto meshName1 = map_databases[i]->getString( "Mesh1" );
        auto meshName2 = map_databases[i]->getString( "Mesh2" );

        // Subset the multmesh for the 2 meshes
        auto mesh1 = manager->Subset( meshName1 );
        auto mesh2 = manager->Subset( meshName2 );
        int inComm = -1;
        if ( mesh1 || mesh2 )
            inComm = 1;

        // Increment the global tag offset by 2
        globalMapTagOffset = manager->getComm().maxReduce( globalMapTagOffset + 2 );

        // Create a comm spanning the meshes
        AMP_MPI mapComm = managerComm.split( inComm );
        if ( inComm == -1 )
            continue;

        // Create the map parameters
        std::shared_ptr<AsyncMapOperatorParameters> mapParams(
            new typename MAP_TYPE::Parameters( map_databases[i] ) );
        mapParams->d_MapComm             = mapComm;
        mapParams->d_Mesh1               = mesh1;
        mapParams->d_Mesh2               = mesh2;
        mapParams->d_BoundaryID1         = map_databases[i]->getScalar<int>( "Surface1" );
        mapParams->d_BoundaryID2         = map_databases[i]->getScalar<int>( "Surface2" );
        mapParams->d_commTag             = globalMapTagOffset;
        mapParams->callMakeConsistentSet = false;

        // Create the map
        std::shared_ptr<AsyncMapOperator> mapOp( new MAP_TYPE( mapParams ) );
        newMapColumn->append( mapOp );
    }

    return newMapColumn;
}
} // namespace Operator
} // namespace AMP
