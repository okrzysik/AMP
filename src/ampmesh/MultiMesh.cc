#include "ampmesh/MultiMesh.h"
#include "ampmesh/MultiIterator.h"
#include "ampmesh/MeshElement.h"
#include "utils/Utilities.h"
#include "utils/Database.h"
#include "utils/MemoryDatabase.h"

#include <set>
#include <vector>
#include <iostream>
#include <string>

namespace AMP {
namespace Mesh {


// Function to check if a string contains another string as the prefix
static bool check_prefix(std::string prefix, std::string str) {
    if ( str.size() < prefix.size() )
        return false;
    if ( str.compare(0,prefix.size(),prefix)==0 )
        return true;
    return false;
}


// Some template functions for vector math
template <class T> 
static inline T vecmin(std::vector<T> x) {
    T min = x[0];
    for (size_t i=1; i<x.size(); i++)
        min = (x[i]<min) ? x[i] : min;
    return min;
}
template <class T> 
static inline T vecmax(std::vector<T> x) {
    T max = x[0];
    for (size_t i=1; i<x.size(); i++)
        max = (x[i]>max) ? x[i] : max;
    return max;
}


// Misc function declerations
static void copyKey(boost::shared_ptr<AMP::Database>&,boost::shared_ptr<AMP::Database>&,std::string,int,int);
static void replaceText(boost::shared_ptr<AMP::Database>&,std::string,std::string);


/********************************************************
* Constructors                                          *
********************************************************/
MultiMesh::MultiMesh( const MeshParameters::shared_ptr &params_in ):
    Mesh(params_in)
{
    // Create a database for each mesh within the multimesh
    AMP_ASSERT(d_db!=NULL);
    std::vector<boost::shared_ptr<AMP::Database> > meshDatabases = MultiMesh::createDatabases(d_db);
    // Create the input parameters and get the approximate number of elements for each mesh
    std::vector<size_t> meshSizes(meshDatabases.size(),0);
    std::vector< boost::shared_ptr<AMP::Mesh::MeshParameters> > meshParameters(meshDatabases.size());
    for (size_t i=0; i<meshDatabases.size(); i++) {
        boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(meshDatabases[i]));
        params->setComm(AMP::AMP_MPI(AMP_COMM_SELF));
        meshParameters[i] = params;
        meshSizes[i] = AMP::Mesh::Mesh::estimateMeshSize(params);
    }
    size_t totalMeshSize = 0;
    for (size_t i=0; i<meshSizes.size(); i++)
        totalMeshSize += meshSizes[i];
    // Determine the load balancing we want to use and create the communicator for the parameters for each mesh
    std::vector<double> weights(meshParameters.size(),0);
    for (size_t i=0; i<meshSizes.size(); i++)
        weights[i] = ((double)meshSizes[i])/((double)totalMeshSize);
    std::vector<AMP_MPI> comms = loadBalancer( weights, 1 );
    // Check that every mesh exist on some comm
    std::vector<int> onComm(meshParameters.size(),0);
    for (size_t i=0; i<meshParameters.size(); i++) {
        if ( !comms[i].isNull() )
            onComm[i] = 1;
    }
    d_comm.maxReduce(&onComm[0],onComm.size());
    for (size_t i=0; i<onComm.size(); i++)
        AMP_ASSERT(onComm[i]==1);
    // Create the meshes
    d_meshes = std::vector<AMP::Mesh::Mesh::shared_ptr>(0);
    for (size_t i=0; i<meshParameters.size(); i++) {
        if ( comms[i].isNull() )
            continue;
        AMP::Mesh::MeshParameters::shared_ptr  params = meshParameters[i];
        params->setComm(comms[i]);
        AMP::Mesh::Mesh::shared_ptr new_mesh = AMP::Mesh::Mesh::buildMesh(params);
        d_meshes.push_back(new_mesh);
    }
    // Get the physical dimension and the highest geometric type
    PhysicalDim = d_meshes[0]->getDim();
    GeomDim = d_meshes[0]->getGeomType();
    for (size_t i=1; i<d_meshes.size(); i++) {
        AMP_INSIST(PhysicalDim==d_meshes[i]->getDim(),"Physical dimension must match for all meshes in multimesh");
        if ( d_meshes[i]->getGeomType() > GeomDim )
            GeomDim = d_meshes[i]->getGeomType();
    }
    GeomDim = (GeomType) d_comm.maxReduce((int)GeomDim);
    // Compute the bounding box of the multimesh
    d_box = d_meshes[0]->getBoundingBox();
    for (size_t i=1; i<d_meshes.size(); i++) {
        std::vector<double> meshBox = d_meshes[i]->getBoundingBox();
        for (int j=0; j<PhysicalDim; j++) {
            if ( meshBox[2*j+0] < d_box[2*j+0] ) { d_box[2*j+0] = meshBox[2*j+0]; }
            if ( meshBox[2*j+1] > d_box[2*j+1] ) { d_box[2*j+1] = meshBox[2*j+1]; }
        }
    }
    // Displace the meshes
    std::vector<double> displacement(PhysicalDim,0.0);
    if ( d_db->keyExists("x_offset") )
        displacement[0] = d_db->getDouble("x_offset");
    if ( d_db->keyExists("y_offset") )
        displacement[1] = d_db->getDouble("y_offset");
    if ( d_db->keyExists("z_offset") )
        displacement[2] = d_db->getDouble("z_offset");
    bool test = false;
    for (size_t i=0; i<displacement.size(); i++) {
        if ( displacement[i] != 0.0 )
            test = true;
    }        
    if ( test )
        displaceMesh(displacement);
}


/********************************************************
* De-constructor                                        *
********************************************************/
MultiMesh::~MultiMesh()
{
}


/********************************************************
* Function to copy the mesh                             *
********************************************************/
Mesh MultiMesh::copy() const
{
    return MultiMesh(*this);
}


/********************************************************
* Function to estimate the mesh size                    *
********************************************************/
size_t MultiMesh::estimateMeshSize( const MeshParameters::shared_ptr &params )
{
    // Create a database for each mesh within the multimesh
    boost::shared_ptr<AMP::Database> database = params->getDatabase();
    AMP_ASSERT(database.get()!=NULL);
    std::vector<boost::shared_ptr<AMP::Database> > meshDatabases = MultiMesh::createDatabases(database);
    // Create the input parameters and get the approximate number of elements for each mesh
    size_t totalMeshSize = 0;
    std::vector< boost::shared_ptr<AMP::Mesh::MeshParameters> > meshParameters(meshDatabases.size());
    for (size_t i=0; i<meshDatabases.size(); i++) {
        boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(meshDatabases[i]));
        params->setComm(AMP::AMP_MPI(AMP_COMM_SELF));
        meshParameters.push_back(params);
        totalMeshSize += AMP::Mesh::Mesh::estimateMeshSize(params);
    }
    return totalMeshSize;
}


/********************************************************
* Function to create the databases for the meshes       *
* within the multimesh.                                 *
********************************************************/
std::vector<boost::shared_ptr<AMP::Database> >  MultiMesh::createDatabases(boost::shared_ptr<AMP::Database> database)
{
    // Find all of the meshes in the database
    AMP_ASSERT(database!=NULL);
    AMP_INSIST(database->keyExists("MeshDatabasePrefix"),"MeshDatabasePrefix must exist in input database");
    AMP_INSIST(database->keyExists("MeshArrayDatabasePrefix"),"MeshArrayDatabasePrefix must exist in input database");
    std::string MeshPrefix = database->getString("MeshDatabasePrefix");
    std::string MeshArrayPrefix = database->getString("MeshArrayDatabasePrefix");
    AMP_ASSERT(!check_prefix(MeshPrefix,MeshArrayPrefix));
    AMP_ASSERT(!check_prefix(MeshArrayPrefix,MeshPrefix));
    std::vector<std::string> keys = database->getAllKeys();
    std::vector<std::string> meshes, meshArrays;
    for (size_t i=0; i<keys.size(); i++) {
        if ( check_prefix(MeshPrefix,keys[i]) ) {
            meshes.push_back(keys[i]);
        } else if ( check_prefix(MeshArrayPrefix,keys[i]) ) {
            meshArrays.push_back(keys[i]);
        }
    }
    // Create the basic databases for each mesh
    std::vector<boost::shared_ptr<AMP::Database> > meshDatabases;
    for (size_t i=0; i<meshes.size(); i++) {
        // We are dealing with a single mesh object (it might be a multimesh), use the existing database
        boost::shared_ptr<AMP::Database> database2 = database->getDatabase(meshes[i]);
        meshDatabases.push_back(database2);
    }
    for (size_t i=0; i<meshArrays.size(); i++) {
        // We are dealing with an array of meshes, create a database for each
        boost::shared_ptr<AMP::Database> database1 = database->getDatabase(meshArrays[i]);
        int N = database1->getInteger("N");
        for (int j=0; j<N; j++) {
            // Create a new database from the existing database
            boost::shared_ptr<AMP::Database> database2( new AMP::MemoryDatabase(meshArrays[i]) );
            keys = database1->getAllKeys();
            for (size_t k=0; k<keys.size(); k++) {
                if ( keys[k].compare("N")==0 || keys[k].compare("iterator")==0 || keys[k].compare("indicies")==0 ) {
                    // These keys are used by the mesh-array and should not be copied
                } else {
                    // We need to copy the key
                    copyKey(database1,database2,keys[k],N,j);
                }
            }
            // Recursively replace all instances of iterator with the index
            if ( database1->keyExists("iterator") ) {
                std::string iterator = database1->getString("iterator");
                AMP_ASSERT(database1->keyExists("indicies"));
                AMP::Database::DataType dataType = database1->getArrayType("indicies");
                std::string index;
                if ( dataType==AMP::Database::AMP_INT ) {
                    std::vector<int> array = database1->getIntegerArray("indicies");
                    AMP_ASSERT((int)array.size()==N);
                    std::stringstream ss;
                    ss << array[j];
                    index = ss.str();
                } else if ( dataType==AMP::Database::AMP_STRING ) {
                    std::vector<std::string> array = database1->getStringArray("indicies");
                    AMP_ASSERT((int)array.size()==N);
                    index = array[j];
                } else {
                    AMP_ERROR("Unknown type for indicies");
                }
                replaceText(database2,iterator,index);
            }
            boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database2));
            meshDatabases.push_back(database2);
        }
    }
    return meshDatabases;
}


/********************************************************
* Return basic mesh info                                *
********************************************************/
size_t MultiMesh::numLocalElements( const GeomType type ) const
{
    // Should we cache this?
    size_t N = 0;
    for (size_t i=0; i<d_meshes.size(); i++)
        N += d_meshes[i]->numLocalElements(type);
    return N;
}
size_t MultiMesh::numGlobalElements( const GeomType type ) const
{
    // Should we cache this?
    size_t N = numLocalElements(type);
    return d_comm.sumReduce(N);
}
size_t MultiMesh::numGhostElements( const GeomType type, int gcw ) const
{
    // Should we cache this?
    size_t N = 0;
    for (size_t i=0; i<d_meshes.size(); i++)
        N += d_meshes[i]->numGhostElements(type,gcw);
    return N;
}


/********************************************************
* Return mesh iterators                                 *
********************************************************/
MeshIterator MultiMesh::getIterator( const GeomType type, const int gcw )
{
    std::vector<boost::shared_ptr<MeshIterator> > iterators(d_meshes.size());
    for (size_t i=0; i<d_meshes.size(); i++) {
        boost::shared_ptr<MeshIterator> iterator_ptr( new MeshIterator(d_meshes[i]->getIterator(type,gcw)) );
        iterators[i] = iterator_ptr;
    }
    return MultiIterator(iterators);
}
std::vector<int> MultiMesh::getIDSets( )
{
    // Get all local id sets
    std::set<int> ids_set;
    for (size_t i=0; i<d_meshes.size(); i++) {
        std::vector<int> mesh_idSet = d_meshes[i]->getIDSets();
        ids_set.insert(mesh_idSet.begin(),mesh_idSet.end());
    }
    std::vector<int> local_ids(ids_set.begin(),ids_set.end());
    // Perform a global communication to syncronize the id sets across all processors
    int N_id_local = (int) local_ids.size();
    std::vector<int> count(d_comm.getSize(),0);
    std::vector<int> disp(d_comm.getSize(),0);
    d_comm.allGather(N_id_local,&count[0]);
    for (int i=1; i<d_comm.getSize(); i++)
        disp[i] = disp[i-1] + count[i-1];
    int N_id_global = disp[d_comm.getSize()-1] + count[d_comm.getSize()-1];
    if ( N_id_global == 0 )
        return std::vector<int>();
    std::vector<int> global_id_list(N_id_global,0);
    int *ptr = NULL;
    if ( N_id_local > 0 )
        ptr = &local_ids[0];
    d_comm.allGather( ptr, N_id_local, &global_id_list[0], &count[0], &disp[0], true );
    // Get the unique set
    for (size_t i=0; i<global_id_list.size(); i++)
        ids_set.insert(global_id_list[i]);
    // Return the final vector of ids
    return std::vector<int>(ids_set.begin(),ids_set.end());
}
MeshIterator MultiMesh::getIDsetIterator( const GeomType type, const int id, const int gcw )
{
    std::vector<boost::shared_ptr<MeshIterator> > iterators;
    iterators.reserve(d_meshes.size());
    for (size_t i=0; i<d_meshes.size(); i++) {
        MeshIterator it = d_meshes[i]->getIDsetIterator(type,id,gcw);
        if ( it.size() > 0 )
            iterators.push_back( boost::shared_ptr<MeshIterator>( new MeshIterator(it) ) );
    }
    return MultiIterator(iterators);
}


/********************************************************
* Function to return the meshID composing the mesh      *
********************************************************/
std::vector<MeshID> MultiMesh::getAllMeshIDs() const
{
    std::set<MeshID> ids;
    ids.insert(d_meshID);
    for (size_t i=0; i<d_meshes.size(); i++) {
        std::vector<MeshID> mesh_ids = d_meshes[i]->getAllMeshIDs();
        for (size_t j=0; j<mesh_ids.size(); j++)
            ids.insert(mesh_ids[j]);
    }
    int send_cnt = ids.size();
    int recv_cnt = d_comm.sumReduce(send_cnt);
    MeshID *send_data = new MeshID[send_cnt];
    MeshID *recv_data = new MeshID[recv_cnt];
    std::set<MeshID>::iterator iterator = ids.begin();
    for (int i=0; i<send_cnt; i++) {
        send_data[i] = *iterator;
        iterator++;
    }
    d_comm.allGather( send_data, send_cnt, recv_data );
    for (int i=0; i<recv_cnt; i++)
        ids.insert(recv_data[i]);
    delete [] send_data;
    delete [] recv_data;
    return std::vector<MeshID>(ids.begin(),ids.end());
}
std::vector<MeshID> MultiMesh::getBaseMeshIDs() const
{
    std::set<MeshID> ids;
    for (size_t i=0; i<d_meshes.size(); i++) {
        std::vector<MeshID> mesh_ids = d_meshes[i]->getBaseMeshIDs();
        for (size_t j=0; j<mesh_ids.size(); j++)
            ids.insert(mesh_ids[j]);
    }
    int send_cnt = ids.size();
    int recv_cnt = d_comm.sumReduce(send_cnt);
    MeshID *send_data = new MeshID[send_cnt];
    MeshID *recv_data = new MeshID[recv_cnt];
    std::set<MeshID>::iterator iterator = ids.begin();
    for (int i=0; i<send_cnt; i++) {
        send_data[i] = *iterator;
        iterator++;
    }
    d_comm.allGather( send_data, send_cnt, recv_data );
    for (int i=0; i<recv_cnt; i++)
        ids.insert(recv_data[i]);
    delete [] send_data;
    delete [] recv_data;
    return std::vector<MeshID>(ids.begin(),ids.end());
}


/********************************************************
* Function to return the mesh with the given ID         *
********************************************************/
boost::shared_ptr<Mesh>  MultiMesh::Subset( MeshID meshID ) 
{
    if ( d_meshID==meshID )
        return shared_from_this();
    for (size_t i=0; i<d_meshes.size(); i++) {
        boost::shared_ptr<Mesh> mesh = d_meshes[i]->Subset(meshID);
        if ( mesh.get()!=NULL )
            return mesh;
    }
    return boost::shared_ptr<Mesh>();
}


/********************************************************
* Function to return the mesh with the given name       *
********************************************************/
boost::shared_ptr<Mesh>  MultiMesh::Subset( std::string name ) {
    if ( d_name==name )
        return shared_from_this();
    for (size_t i=0; i<d_meshes.size(); i++) {
        boost::shared_ptr<Mesh> mesh = d_meshes[i]->Subset(name);
        if ( mesh.get()!=NULL )
            return mesh;
    }
    return boost::shared_ptr<Mesh>();
}


/********************************************************
* Displace a mesh by a scalar ammount                   *
********************************************************/
void MultiMesh::displaceMesh( std::vector<double> x_in )
{
    // Check x
    AMP_INSIST((short int)x_in.size()==PhysicalDim,"Displacement vector size should match PhysicalDim");
    std::vector<double> x = x_in;
    d_comm.minReduce(&x[0],x.size());
    for (size_t i=0; i<x.size(); i++)
        AMP_INSIST(fabs(x[i]-x_in[i])<1e-12,"x does not match on all processors");
    // Displace the meshes
    for (size_t i=0; i<d_meshes.size(); i++)
        d_meshes[i]->displaceMesh(x);
    // Update the bounding box
    for (int i=0; i<PhysicalDim; i++) {
        d_box[2*i+0] += x[i];
        d_box[2*i+1] += x[i];
    }
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


/********************************************************
* Function to recursively replace all string values     *
* that match a substring.                               *
********************************************************/
static void replaceSubString( std::string& string, std::string search, std::string replace )
{
    while ( 1 ) {
        size_t pos = string.find(search);
        if ( pos == std::string::npos ) 
            break;
        string.replace( pos, search.size(), replace );
    }
}
static void replaceText( boost::shared_ptr<AMP::Database>& database, std::string search, std::string replace )
{
    std::vector<std::string> key = database->getAllKeys();
    for (size_t i=0; i<key.size(); i++) {
        AMP::Database::DataType type = database->getArrayType(key[i]);
        switch (type) {
            case AMP::Database::AMP_DATABASE: {
                // Search the database
                boost::shared_ptr<AMP::Database> database2 = database->getDatabase(key[i]);
                replaceText( database2, search, replace );
                } break;
            case AMP::Database::AMP_STRING: {
                // Search the string
                std::vector<std::string> data = database->getStringArray(key[i]);
                for (size_t j=0; j<data.size(); j++) {
                    replaceSubString(data[j],search,replace);
                }
                database->putStringArray(key[i],data);
                } break;
            // Cases that we don't need to do anything for
            case AMP::Database::AMP_CHAR:
            case AMP::Database::AMP_INVALID:
            case AMP::Database::AMP_BOOL:
            case AMP::Database::AMP_INT:
            case AMP::Database::AMP_COMPLEX:
            case AMP::Database::AMP_DOUBLE:
            case AMP::Database::AMP_FLOAT:
            case AMP::Database::AMP_BOX:
                break;
            // Unknown case
            default:
                AMP_ERROR("Unknown key type");
        }
    }
}


/********************************************************
* Function to create the sub-communicators              *
********************************************************/
std::vector<AMP_MPI> MultiMesh::loadBalancer( std::vector<double> &weights, int method)
{
    // Deal with the special cases directly
    if ( weights.size()<=1 || d_comm.getSize()==1 )
        return std::vector<AMP_MPI>(weights.size(),d_comm);
    // Choose the appropriate load balancer
    if ( method == 0 ) {
        // Everybody is on the same communicator
        return std::vector<AMP_MPI>(weights.size(),d_comm);
    } else if ( method == 1 ) {
        // We want to split the meshes onto independent processors
        std::vector<std::pair<double,int> >  ids(weights.size());
        for (size_t i=0; i<weights.size(); i++)
            ids[i] = std::pair<double,int>( weights[i], (int) i );
        std::vector<comm_groups> groups;
        if ( d_comm.getSize() >= (int)weights.size() )
            groups = independentGroups1( d_comm.getSize(), ids );
        else
            groups = independentGroups2( d_comm.getSize(), ids );
        // Split the comm into the appropriate groups
        std::vector<int> cum_N_proc(groups.size());
        cum_N_proc[0] = groups[0].N_procs;
        for (size_t i=1; i<groups.size(); i++)
            cum_N_proc[i] = cum_N_proc[i-1] + groups[i].N_procs;
        int myid = AMP::Utilities::findfirst(cum_N_proc,d_comm.getRank()+1);
        AMP_MPI myComm = d_comm.split(myid,d_comm.getRank());
        std::vector<AMP_MPI> newComms(weights.size(),AMP_MPI(AMP_COMM_NULL));
        for (size_t i=0; i<groups[myid].ids.size(); i++)
            newComms[groups[myid].ids[i]] = myComm;
        return newComms;
    } else if ( method == 2 ) {
        // We want to try to achieve a balanced approach to the load balance
        AMP_ERROR("This load balancer is not implimented yet");
    } else {
        AMP_ERROR("Unknown load balancer");
    }
    return std::vector<AMP_MPI>(0);
}


/********************************************************
* Functions to solve simple load balancing computations *
********************************************************/
std::vector<MultiMesh::comm_groups>  MultiMesh::independentGroups1(int N_procs, std::vector<std::pair<double,int> >  &ids_in)
{
    // This will distribute the groups onto the processors such that no groups share any processors
    // Note: this requires the number of processors to be >= the number of ids
    AMP_ASSERT(N_procs>=(int)ids_in.size());
    std::vector<std::pair<double,int> > ids = ids_in;
    AMP::Utilities::quicksort(ids);
    // Start with the smallest mesh and start assigning processors
    double weight_remaining = 0.0;
    for (size_t i=0; i<ids.size(); i++)
        weight_remaining += ids[i].first;
    int N_procs_remaining = N_procs;
    std::vector<comm_groups> groups(ids.size());
    for (size_t i=0; i<ids.size(); i++) {
        double N_proc_d = ids[i].first/weight_remaining * ((double) N_procs_remaining);
        int N_proc = (int) (N_proc_d+0.5);      // Round to the nearest processor count
        if ( N_proc < 1 ) { N_proc=1; }         // We require at least on processor
        if ( (N_procs_remaining-N_proc) < (int)(ids.size()-i-1) ) 
            N_proc--;       // We need at least the 1 processors for each of the remaining meshes
        groups[i].N_procs = N_proc;
        groups[i].ids = std::vector<int>(1,ids[i].second);
        N_procs_remaining -= N_proc;
        weight_remaining -= ids[i].first;
    }
    return groups;
}
std::vector<MultiMesh::comm_groups>  MultiMesh::independentGroups2(int N_procs, std::vector<std::pair<double,int> >  &ids_in)
{
    // This will distribute the groups onto the processors such that each group has exactly one processor.  
    // Note: this requires the number of processors to be <= the number of ids
    AMP_ASSERT(N_procs<=(int)ids_in.size());
    std::vector<std::pair<double,int> > ids = ids_in;
    double total_weight = 0.0;
    for (size_t i=0; i<ids.size(); i++)
        total_weight += ids[i].first;
    double weight_avg = total_weight / ((double) N_procs);
    std::vector<comm_groups> groups;
    // Remove any ids that require a processor by themselves
    std::vector<std::pair<double,int> >::iterator  iterator = ids.begin();
    while ( iterator != ids.end() ) {
        if ( iterator->first >= weight_avg ) {
            comm_groups tmp;
            tmp.N_procs = 1;
            tmp.ids = std::vector<int>(1,iterator->second);
            groups.push_back(tmp);
            iterator = ids.erase(iterator);
        } else {
            ++iterator;
        }
    }
    if ( groups.size() == 0 ) {
        // We did not remove any elements, lets begin grouping them together
        if ( ( ids[0].first + ids[ids.size()-1].first ) > weight_avg ) {
            // The first and last elements create a valid bin, add them together
            comm_groups tmp;
            tmp.N_procs = 1;
            tmp.ids = std::vector<int>(2);
            tmp.ids[0] = ids[0].second;
            tmp.ids[1] = ids[ids.size()-1].second;
            ids.erase(ids.begin());     // erase the first element
            ids.resize(ids.size()-1);   // erase the last element
            groups.push_back(tmp);      // add the new element
        } else {
            // Combine the first set of elements that sum to the desired value
            int n = 1;
            double weight_sum = ids[0].first;
            while ( weight_sum < weight_avg ) {
                weight_sum += ids[n].first;
                n++;
            }
            comm_groups tmp;
            tmp.N_procs = 1;
            tmp.ids = std::vector<int>(n);
            for (int i=0; i<n; i++)
                tmp.ids[i] = ids[i].second;
            ids.erase(ids.begin(),ids.begin()+n);
            groups.push_back(tmp);
        }
    }
    // Recursively add the remaining ids
    if ( ids.size() > 0 ) {
        std::vector<MultiMesh::comm_groups> groups2 = independentGroups2( N_procs-groups.size(), ids );
        groups.insert( groups.end(), groups2.begin(), groups2.end() );
    }
    return groups;
}


} // Mesh namespace
} // AMP namespace

