#include "operators/map/AsyncMapOperator.h"
#include "operators/map/AsyncMapOperatorParameters.h"
//#include "operators/map/Map3to1to3.h"
#include "utils/MemoryDatabase.h"
#include "utils/Database.h"

#include <iostream>

namespace AMP {
namespace Operator {


static std::vector<boost::shared_ptr<AMP::Database> >  createDatabases(boost::shared_ptr<AMP::Database>);

enum  MapDominance { Master , Slave };
enum  MapConstructionType { Synchronous , Asynchronous };
struct MapConstructionParam {
    int                                 boundaryId;
    AMP_MPI                             comm;
    MapConstructionType                 construction;
    MapDominance                        dominance;
    boost::shared_ptr<AMP::Database>    database;
    std::string                         meshName;
    AMP::Mesh::Mesh::shared_ptr         mesh;
    std::string                         mapType;
    size_t                              tagOffset;
};

extern size_t globalMapTagOffset;      // We need a global unique tag offset for every map


template <typename MAP_TYPE>
boost::shared_ptr<AsyncMapColumnOperator>  AsyncMapColumnOperator::build ( AMP::Mesh::Mesh::shared_ptr manager, boost::shared_ptr<Database> database )
{

    typedef boost::shared_ptr < AsyncMapOperator >             AsyncOp_ptr;
    typedef boost::shared_ptr < AsyncMapOperatorParameters >   AsyncOpParams_ptr;
    // Map3to1to3::spMap::type       sharedMap ( new std::multimap<double,double> );
    //sharedMap.isComputed = false;

    std::vector < AsyncOp_ptr >           asyncToDo;
    std::vector < AsyncOpParams_ptr >     asyncParams;

    boost::shared_ptr<AsyncMapColumnOperatorParameters>  newParams ( new AsyncMapColumnOperatorParameters ( database ) );
    boost::shared_ptr<AsyncMapColumnOperator>  newMapColumn ( new AsyncMapColumnOperator ( newParams ) );

    // Check that the map type matches the maps we are building
    AMP_ASSERT(database->keyExists("MapType"));
    AMP_ASSERT( MAP_TYPE::validMapType( database->getString("MapType") ) );

    // Get the number of maps in the database
    AMP_INSIST(database->keyExists("N_maps"),"N_maps must exist in input database");
    int N_maps = database->getInteger("N_maps");

    // Create the databases for the individual maps
    std::vector<boost::shared_ptr<AMP::Database> >  map_databases = createDatabases(database);

    // Loop through the maps
    AMP_MPI managerComm = manager->getComm();
    for ( int i=0; i<N_maps; i++ ) {

        // Get the names of the 2 meshes involved
        std::string meshName1 = map_databases[i]->getString("Mesh1");
        std::string meshName2 = map_databases[i]->getString("Mesh2");

        // Subset the multmesh for the 2 meshes
        boost::shared_ptr<AMP::Mesh::Mesh> mesh1 = manager->Subset(meshName1);
        boost::shared_ptr<AMP::Mesh::Mesh> mesh2 = manager->Subset(meshName2);
        int inComm = -1;
        if ( mesh1 || mesh2 )
            inComm = 1;

        // Increment the global tag offset by 2
        globalMapTagOffset = manager->getComm().maxReduce(globalMapTagOffset+2);

        // Create a comm spanning the meshes
        AMP_MPI mapComm = managerComm.split(inComm);
        if ( inComm==-1 )
            continue;

        // Create the map parameters
        boost::shared_ptr<AsyncMapOperatorParameters>  mapParams( new typename MAP_TYPE::Parameters( map_databases[i] ) );
        mapParams->d_MapComm = mapComm;
        mapParams->d_Mesh1 = mesh1;
        mapParams->d_Mesh2 = mesh2;
        mapParams->d_BoundaryID1 = map_databases[i]->getInteger("Surface1");
        mapParams->d_BoundaryID2 = map_databases[i]->getInteger("Surface2");

        // Create the map
        AsyncOp_ptr  mapOp ( new MAP_TYPE ( mapParams ) );
        //boost::shared_ptr<Map3to1to3>  curM313;
        //if ( curM313 = boost::dynamic_pointer_cast<Map3to1to3> ( mapOp ) )
        //    curM313->getMapData() = sharedMap;
        newMapColumn->append ( mapOp );

    }

    return newMapColumn;

}


/********************************************************
* Function to copy a key from database 1 to database 2  *
* If the key is an array of size N, it will only copy   *
* the ith value.                                        *
********************************************************/
static void copyKey(boost::shared_ptr<AMP::Database> &database1, 
    boost::shared_ptr<AMP::Database> &database2, std::string key, int N, int i ) 
{
    int size = database1->getArraySize(key);
    AMP::Database::DataType type = database1->getArrayType(key);
    switch (type) {
        case AMP::Database::AMP_INVALID: {
            // We don't know what this is
            AMP_ERROR("Invalid database object");
            } break;
        case AMP::Database::AMP_DATABASE: {
            // Copy the database
            AMP_ERROR("Not programmed for databases yet");
            } break;
        case AMP::Database::AMP_BOOL: {
            // Copy a bool
            std::vector<unsigned char> data = database1->getBoolArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            if ( N == size )
                database2->putBool(key,data[i]);
            else
                database2->putBoolArray(key,data);
            } break;
        case AMP::Database::AMP_CHAR: {
            // Copy a char
            std::vector<char> data = database1->getCharArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            if ( N == size ) {
                database2->putChar(key,data[i]);
            } else {
                // We need to try a search and replace
                database2->putCharArray(key,data);
            }
            } break;
        case AMP::Database::AMP_INT: {
            // Copy an int
            std::vector<int> data = database1->getIntegerArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            if ( N == size )
                database2->putInteger(key,data[i]);
            else
                database2->putIntegerArray(key,data);
            } break;
        case AMP::Database::AMP_COMPLEX: {
            // Copy a complex number
            std::vector< std::complex<double> > data = database1->getComplexArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            if ( N == size )
                database2->putComplex(key,data[i]);
            else
                database2->putComplexArray(key,data);
            } break;
        case AMP::Database::AMP_DOUBLE: {
            // Copy a double number
            std::vector<double> data = database1->getDoubleArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            if ( N == size )
                database2->putDouble(key,data[i]);
            else
                database2->putDoubleArray(key,data);
            } break;
        case AMP::Database::AMP_FLOAT: {
            // Copy a float
            std::vector<float> data = database1->getFloatArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            if ( N == size )
                database2->putFloat(key,data[i]);
            else
                database2->putFloatArray(key,data);
            } break;
        case AMP::Database::AMP_STRING: {
            // Copy a string
            std::vector<std::string> data = database1->getStringArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            if ( N == size )
                database2->putString(key,data[i]);
            else
                database2->putStringArray(key,data);
            } break;
        case AMP::Database::AMP_BOX: {
            // Copy a box
            AMP_ERROR("Not programmed for boxes yet");
            } break;
        default:
            AMP_ERROR("Unknown key type");
    }
}


/************************************************************
* Function to create the databases for the individual maps  *
************************************************************/
static std::vector<boost::shared_ptr<AMP::Database> >  createDatabases(boost::shared_ptr<AMP::Database> database1)
{
    AMP_INSIST(database1->keyExists("N_maps"),"N_maps must exist in input database");
    int N_maps = database1->getInteger("N_maps");
    // Create the basic databases for each mesh
    std::vector<boost::shared_ptr<AMP::Database> > meshDatabases;
    meshDatabases.reserve(N_maps);
    for (int i=0; i<N_maps; i++) {
        // Create a new database from the existing database
        boost::shared_ptr<AMP::Database> database2( new AMP::MemoryDatabase("MapDatabase") );
        std::vector<std::string> keys = database1->getAllKeys();
        for (size_t k=0; k<keys.size(); k++) {
            if ( keys[k].compare("N_maps")==0 ) {
                // These keys are used by the builder and should not be in the sub database
            } else {
                // We need to copy the key
                copyKey(database1,database2,keys[k],N_maps,i);
            }
        }
        meshDatabases.push_back(database2);
    }
    return meshDatabases;
}


}
}

