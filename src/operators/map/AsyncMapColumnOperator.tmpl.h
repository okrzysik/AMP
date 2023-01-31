#include "AMP/operators/map/AsyncMapOperator.h"
#include "AMP/operators/map/AsyncMapOperatorParameters.h"
#include "AMP/utils/Database.h"


namespace AMP::Operator {


enum class MapDominance { Master, Slave };
enum class MapConstructionType { Synchronous, Asynchronous };
struct MapConstructionParam {
    int boundaryId;
    AMP_MPI comm;
    MapConstructionType construction;
    MapDominance dominance;
    std::shared_ptr<AMP::Database> database;
    std::string meshName;
    std::shared_ptr<AMP::Mesh::Mesh> mesh;
    std::string mapType;
    size_t tagOffset;
};

extern size_t globalMapTagOffset; // We need a global unique tag offset for every map


template<typename MAP_TYPE>
std::shared_ptr<AsyncMapColumnOperator>
AsyncMapColumnOperator::build( std::shared_ptr<AMP::Mesh::Mesh> manager,
                               std::shared_ptr<Database> database )
{

    auto newParams    = std::make_shared<AsyncMapColumnOperatorParameters>( database );
    auto newMapColumn = std::make_shared<AsyncMapColumnOperator>( newParams );

    // Check that the map type matches the maps we are building
    AMP_ASSERT( database->keyExists( "MapType" ) );
    AMP_ASSERT( MAP_TYPE::validMapType( database->getString( "MapType" ) ) );

    // Create the databases for the individual maps
    auto map_databases = createDatabases( database );

    // Loop through the maps
    AMP_MPI managerComm = manager->getComm();
    for ( auto db : map_databases ) {
        // Get the names of the 2 meshes involved
        auto meshName1 = db->getString( "Mesh1" );
        auto meshName2 = db->getString( "Mesh2" );

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
        auto mapParams                   = std::make_shared<typename MAP_TYPE::Parameters>( db );
        mapParams->d_MapComm             = mapComm;
        mapParams->d_Mesh1               = mesh1;
        mapParams->d_Mesh2               = mesh2;
        mapParams->d_BoundaryID1         = db->getScalar<int>( "Surface1" );
        mapParams->d_BoundaryID2         = db->getScalar<int>( "Surface2" );
        mapParams->d_commTag             = globalMapTagOffset;
        mapParams->callMakeConsistentSet = false;

        // Create the map
        auto mapOp = std::make_shared<MAP_TYPE>( mapParams );
        newMapColumn->append( mapOp );
    }

    return newMapColumn;
}

} // namespace AMP::Operator
