#include "ampmesh/MultiMesh.h"
#include "ampmesh/MultiMeshParameters.h"
#include "ampmesh/SubsetMesh.h"
#include "ampmesh/MultiIterator.h"
#include "ampmesh/MeshElement.h"
#include "ampmesh/loadBalance.h"

#include "utils/Utilities.h"
#include "utils/Database.h"
#include "utils/MemoryDatabase.h"

#ifdef USE_AMP_VECTORS
    #include "vectors/Vector.h"
    #include "vectors/MultiVector.h"
#endif

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
static void copyKey(AMP::Database::shared_ptr,std::vector<AMP::Database::shared_ptr>,
    const std::string&,bool,const std::string&,const std::vector<std::string>&);
static void replaceText(boost::shared_ptr<AMP::Database>&,const std::string&,const std::string&);
static void replaceSubString( std::string&, const std::string&, const std::string& );


/********************************************************
* Constructors                                          *
********************************************************/
MultiMesh::MultiMesh( const MeshParameters::shared_ptr &params_in ):
    Mesh(params_in)
{
    // Create an array of MeshParameters for each submesh
    AMP_ASSERT(d_db!=NULL);
    boost::shared_ptr<MultiMeshParameters> params = boost::dynamic_pointer_cast<MultiMeshParameters>(params_in);
    if ( params.get()==NULL ) {
        boost::shared_ptr<AMP::Database> database = params_in->getDatabase();
        // Create a database for each mesh within the multimesh
        std::vector<boost::shared_ptr<AMP::Database> > meshDatabases = MultiMesh::createDatabases(database);
        params = boost::shared_ptr<MultiMeshParameters>( new MultiMeshParameters(database) );
        params->params = std::vector<MeshParameters::shared_ptr>(meshDatabases.size());
        params->N_elements = std::vector<size_t>(meshDatabases.size());
        for (size_t i=0; i<meshDatabases.size(); i++) {
            params->params[i] = AMP::Mesh::MeshParameters::shared_ptr( new AMP::Mesh::MeshParameters(meshDatabases[i]) );
            params->N_elements[i] = Mesh::estimateMeshSize( params->params[i] );
        }
        params->method = 1;
    }
    // Determine the load balancing we want to use and create the communicator for the parameters for each mesh
    std::vector<AMP_MPI> comms = createComms( loadBalancer( d_comm.getSize(), params->params, params->N_elements, params->method ) );
    // Check that every mesh exist on some comm
    std::vector<int> onComm(comms.size(),0);
    for (size_t i=0; i<comms.size(); i++) {
        if ( !comms[i].isNull() )
            onComm[i] = 1;
    }
    d_comm.maxReduce(&onComm[0],(int)onComm.size());
    for (size_t i=0; i<onComm.size(); i++)
        AMP_ASSERT(onComm[i]==1);
    // Create the meshes
    d_meshes = std::vector<AMP::Mesh::Mesh::shared_ptr>(0);
    for (size_t i=0; i<comms.size(); i++) {
        if ( comms[i].isNull() )
            continue;
        params->params[i]->setComm(comms[i]);
        AMP::Mesh::Mesh::shared_ptr new_mesh = AMP::Mesh::Mesh::buildMesh(params->params[i]);
        d_meshes.push_back(new_mesh);
    }
    // Get the physical dimension and the highest geometric type
    PhysicalDim = d_meshes[0]->getDim();
    GeomDim = d_meshes[0]->getGeomType();
    d_max_gcw = 0;
    for (size_t i=1; i<d_meshes.size(); i++) {
        AMP_INSIST(PhysicalDim==d_meshes[i]->getDim(),"Physical dimension must match for all meshes in multimesh");
        if ( d_meshes[i]->getGeomType() > GeomDim )
            GeomDim = d_meshes[i]->getGeomType();
        if ( d_meshes[i]->getMaxGhostWidth() > d_max_gcw )
            d_max_gcw = d_meshes[i]->getMaxGhostWidth();
    }
    GeomDim = (GeomType) d_comm.maxReduce((int)GeomDim);
    d_max_gcw = d_comm.maxReduce(d_max_gcw);
    // Compute the bounding box of the multimesh
    d_box_local = d_meshes[0]->getBoundingBox();
    for (size_t i=1; i<d_meshes.size(); i++) {
        std::vector<double> meshBox = d_meshes[i]->getBoundingBox();
        for (int j=0; j<PhysicalDim; j++) {
            if ( meshBox[2*j+0] < d_box_local[2*j+0] ) { d_box_local[2*j+0] = meshBox[2*j+0]; }
            if ( meshBox[2*j+1] > d_box_local[2*j+1] ) { d_box_local[2*j+1] = meshBox[2*j+1]; }
        }
    }
    d_box = std::vector<double>(PhysicalDim*2);
    for (int i=0; i<PhysicalDim; i++) {
        d_box[2*i+0] = d_comm.minReduce( d_box_local[2*i+0] );
        d_box[2*i+1] = d_comm.maxReduce( d_box_local[2*i+1] );
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
MultiMesh::MultiMesh ( const AMP_MPI &comm, const std::vector<Mesh::shared_ptr> &meshes )
{
    d_comm = comm;
    this->setMeshID();
    // Get the list of non-null meshes
    d_meshes = std::vector<Mesh::shared_ptr>();
    for (size_t i=0; i<meshes.size(); i++) {
        if ( meshes[i].get()!=NULL )
            d_meshes.push_back( meshes[i] );
    }
    if ( d_comm.sumReduce(d_meshes.size())==0 ) {
        AMP_ERROR("Empty multimeshes have not been tested yet");
    }
    // Check the comm (note: the order for the comparison matters)
    for (size_t i=0; i<d_meshes.size(); i++) {
        AMP_ASSERT( d_meshes[i]->getComm() <= d_comm );
    }
    // Get the physical dimension and the highest geometric type
    PhysicalDim = d_meshes[0]->getDim();
    GeomDim = d_meshes[0]->getGeomType();
    d_max_gcw = 0;
    for (size_t i=1; i<d_meshes.size(); i++) {
        AMP_INSIST(PhysicalDim==d_meshes[i]->getDim(),"Physical dimension must match for all meshes in multimesh");
        if ( d_meshes[i]->getGeomType() > GeomDim )
            GeomDim = d_meshes[i]->getGeomType();
        if ( d_meshes[i]->getMaxGhostWidth() > d_max_gcw )
            d_max_gcw = d_meshes[i]->getMaxGhostWidth();
    }
    GeomDim = (GeomType) d_comm.maxReduce((int)GeomDim);
    d_max_gcw = d_comm.maxReduce(d_max_gcw);
    // Compute the bounding box of the multimesh
    d_box_local = d_meshes[0]->getBoundingBox();
    for (size_t i=1; i<d_meshes.size(); i++) {
        std::vector<double> meshBox = d_meshes[i]->getBoundingBox();
        for (int j=0; j<PhysicalDim; j++) {
            if ( meshBox[2*j+0] < d_box_local[2*j+0] ) { d_box_local[2*j+0] = meshBox[2*j+0]; }
            if ( meshBox[2*j+1] > d_box_local[2*j+1] ) { d_box_local[2*j+1] = meshBox[2*j+1]; }
        }
    }
    d_box = std::vector<double>(PhysicalDim*2);
    for (int i=0; i<PhysicalDim; i++) {
        d_box[2*i+0] = d_comm.minReduce( d_box_local[2*i+0] );
        d_box[2*i+1] = d_comm.maxReduce( d_box_local[2*i+1] );
    } 
}


/********************************************************
* Function to simulate the mesh building process        *
********************************************************/
AMP::Mesh::LoadBalance  MultiMesh::simulateBuildMesh( const MeshParameters::shared_ptr params, const std::vector<int> &comm_ranks ) 
{
    // Create the multimesh parameters
    boost::shared_ptr<MultiMeshParameters> multimeshParams = boost::dynamic_pointer_cast<MultiMeshParameters>( params );
    if ( multimeshParams.get()==NULL ) {
        boost::shared_ptr<AMP::Database> database = params->getDatabase();
        multimeshParams = boost::shared_ptr<MultiMeshParameters>( new MultiMeshParameters(database) );
        // Create a database for each mesh within the multimesh
        std::vector<boost::shared_ptr<AMP::Database> > meshDatabases = MultiMesh::createDatabases(database);
        multimeshParams->params = std::vector<MeshParameters::shared_ptr>(meshDatabases.size());
        multimeshParams->N_elements = std::vector<size_t>(meshDatabases.size());
        std::vector<int> rank1(1,0);
        for (size_t i=0; i<meshDatabases.size(); i++)  {
            AMP::Mesh::MeshParameters::shared_ptr meshParam( new AMP::Mesh::MeshParameters(meshDatabases[i]) );
            LoadBalance  mesh( meshParam, rank1 );
            multimeshParams->params[i] = mesh.getParams();
            multimeshParams->N_elements[i] = mesh.getSize();
        }
        multimeshParams->method = 1;
    }
    // Determine the load balancing we want to use and create the virtual communicators for each mesh
    std::vector<rank_list> groups = loadBalancer( comm_ranks.size(), multimeshParams->params, multimeshParams->N_elements, multimeshParams->method );
    AMP_ASSERT(groups.size()==multimeshParams->params.size());
    size_t N_proc_groups = 0;
    for (size_t i=0; i<groups.size(); i++) {
        AMP_ASSERT(groups[i].size()>0);
        for (size_t j=0; j<groups[i].size(); j++)
            groups[i][j] = comm_ranks[groups[i][j]];
        N_proc_groups += groups[i].size();
    }
    int decomp = (N_proc_groups==comm_ranks.size()) ? 1:0;
    // Create the simulated mesh structure
    std::vector<LoadBalance> submeshes(groups.size());
    for (size_t i=0; i<groups.size(); i++)
        submeshes[i] = LoadBalance( multimeshParams->params[i], groups[i], multimeshParams->N_elements[i] );
    return  LoadBalance( multimeshParams, comm_ranks, submeshes, decomp );
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
    // Create the multimesh parameters
    boost::shared_ptr<MultiMeshParameters> multimeshParams = boost::dynamic_pointer_cast<MultiMeshParameters>( params );
    if ( multimeshParams.get()==NULL ) {
        boost::shared_ptr<AMP::Database> database = params->getDatabase();
        multimeshParams = boost::shared_ptr<MultiMeshParameters>( new MultiMeshParameters(database) );
        // Create a database for each mesh within the multimesh
        std::vector<boost::shared_ptr<AMP::Database> > meshDatabases = MultiMesh::createDatabases(database);
        multimeshParams->params = std::vector<MeshParameters::shared_ptr>(meshDatabases.size());
        for (size_t i=0; i<meshDatabases.size(); i++) 
            multimeshParams->params[i] = AMP::Mesh::MeshParameters::shared_ptr( new AMP::Mesh::MeshParameters(meshDatabases[i]) );
    }
    // Get the approximate number of elements for each mesh
    size_t totalMeshSize = 0;
    for (size_t i=0; i<multimeshParams->params.size(); i++) {
        size_t localMeshSize = AMP::Mesh::Mesh::estimateMeshSize(multimeshParams->params[i]);
        AMP_ASSERT(localMeshSize>0);
        totalMeshSize += localMeshSize;
    }
    // Adjust the number of elements by a weight if desired
    if ( params->getDatabase()->keyExists("Weight") ) {
        double weight = params->getDatabase()->getDouble("Weight");
        totalMeshSize = (size_t) ceil(weight*((double)totalMeshSize));
    }
    return totalMeshSize;
}


/********************************************************
* Function to create the databases for the meshes       *
* within the multimesh.                                 *
********************************************************/
std::vector<boost::shared_ptr<AMP::Database> >  MultiMesh::createDatabases(boost::shared_ptr<AMP::Database> database)
{
    // We might have already created and stored the databases for each mesh
    if ( database->keyExists("submeshDatabases") ) {
        std::vector<std::string> databaseNames = database->getStringArray("submeshDatabases");
        AMP_ASSERT(databaseNames.size()>0);
        std::vector<boost::shared_ptr<AMP::Database> > meshDatabases(databaseNames.size());
        for (size_t i=0; i<databaseNames.size(); i++)
            meshDatabases[i] = database->getDatabase(databaseNames[i]);
        return meshDatabases;
    }
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
        // Get the iterator and indicies
        std::string iterator;
        std::vector<std::string> index(N);
        if ( database1->keyExists("iterator") ) {
            std::string iterator = database1->getString("iterator");
            AMP_ASSERT(database1->keyExists("indicies"));
            AMP::Database::DataType dataType = database1->getArrayType("indicies");
            std::vector<std::string> index(N);
            if ( dataType==AMP::Database::AMP_INT ) {
                std::vector<int> array = database1->getIntegerArray("indicies");
                AMP_ASSERT((int)array.size()==N);
                for (int j=0; j<N; j++) {
                    std::stringstream ss;
                    ss << array[j];
                    index[j] = ss.str();
                }
            } else if ( dataType==AMP::Database::AMP_STRING ) {
                index = database1->getStringArray("indicies");
            } else {
                AMP_ERROR("Unknown type for indicies");
            }
        }
        // Create the new databases 
        std::vector<AMP::Database::shared_ptr> databaseArray(N);
        for (int j=0; j<N; j++)
            databaseArray[j] = boost::shared_ptr<AMP::Database>( new AMP::MemoryDatabase(meshArrays[i]) );
        // Populate the databases with the proper keys
        keys = database1->getAllKeys();
        for (size_t k=0; k<keys.size(); k++) {
            if ( keys[k].compare("N")==0 || keys[k].compare("iterator")==0 || keys[k].compare("indicies")==0 ) {
                // These keys are used by the mesh-array and should not be copied
            } else if ( keys[k].compare("Size")==0 || keys[k].compare("Range")==0 ) {
                // These are special keys that should not be divided
                copyKey(database1,databaseArray,keys[k],false,std::string(),index);
            } else {
                // We need to copy the key (and possibly replace the iterator)
                copyKey(database1,databaseArray,keys[k],true,iterator,index);
            }
        }
        // Add the new databases to meshDatabases
        for (int j=0; j<N; j++) {
            meshDatabases.push_back(databaseArray[j]);
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
std::vector<Mesh::shared_ptr> MultiMesh::getMeshes() const
{
    return d_meshes;
}


/********************************************************
* Return mesh iterators                                 *
********************************************************/
MeshIterator MultiMesh::getIterator( const GeomType type, const int gcw ) const
{
    std::vector<boost::shared_ptr<MeshIterator> > iterators(d_meshes.size());
    for (size_t i=0; i<d_meshes.size(); i++) {
        boost::shared_ptr<MeshIterator> iterator_ptr( new MeshIterator(d_meshes[i]->getIterator(type,gcw)) );
        iterators[i] = iterator_ptr;
    }
    return MultiIterator(iterators);
}
MeshIterator MultiMesh::getSurfaceIterator ( const GeomType type, const int gcw ) const
{
    std::vector<boost::shared_ptr<MeshIterator> > iterators(d_meshes.size());
    for (size_t i=0; i<d_meshes.size(); i++) {
        boost::shared_ptr<MeshIterator> iterator_ptr( new MeshIterator(d_meshes[i]->getSurfaceIterator(type,gcw)) );
        iterators[i] = iterator_ptr;
    }
    return MultiIterator(iterators);
}
std::vector<int> MultiMesh::getBoundaryIDs( ) const
{
    // Get all local id sets
    std::set<int> ids_set;
    for (size_t i=0; i<d_meshes.size(); i++) {
        std::vector<int> mesh_idSet = d_meshes[i]->getBoundaryIDs();
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
MeshIterator MultiMesh::getBoundaryIDIterator( const GeomType type, const int id, const int gcw ) const
{
    std::vector<boost::shared_ptr<MeshIterator> > iterators;
    iterators.reserve(d_meshes.size());
    for (size_t i=0; i<d_meshes.size(); i++) {
        MeshIterator it = d_meshes[i]->getBoundaryIDIterator(type,id,gcw);
        if ( it.size() > 0 )
            iterators.push_back( boost::shared_ptr<MeshIterator>( new MeshIterator(it) ) );
    }
    return MultiIterator(iterators);
}
std::vector<int> MultiMesh::getBlockIDs( ) const
{
    // Get all local id sets
    std::set<int> ids_set;
    for (size_t i=0; i<d_meshes.size(); i++) {
        std::vector<int> mesh_idSet = d_meshes[i]->getBlockIDs();
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
MeshIterator MultiMesh::getBlockIDIterator( const GeomType type, const int id, const int gcw ) const
{
    std::vector<boost::shared_ptr<MeshIterator> > iterators;
    iterators.reserve(d_meshes.size());
    for (size_t i=0; i<d_meshes.size(); i++) {
        MeshIterator it = d_meshes[i]->getBlockIDIterator(type,id,gcw);
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
    std::vector<MeshID> tmp = this->getLocalMeshIDs();
    std::set<MeshID> ids(tmp.begin(),tmp.end());
    int send_cnt = (int) ids.size();
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
    std::vector<MeshID> tmp = this->getLocalBaseMeshIDs();
    std::set<MeshID> ids(tmp.begin(),tmp.end());
    int send_cnt = (int) ids.size();
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
std::vector<MeshID> MultiMesh::getLocalMeshIDs() const
{
    std::set<MeshID> ids;
    ids.insert(d_meshID);
    for (size_t i=0; i<d_meshes.size(); i++) {
        std::vector<MeshID> mesh_ids = d_meshes[i]->getLocalMeshIDs();
        for (size_t j=0; j<mesh_ids.size(); j++)
            ids.insert(mesh_ids[j]);
    }
    return std::vector<MeshID>(ids.begin(),ids.end());
}
std::vector<MeshID> MultiMesh::getLocalBaseMeshIDs() const
{
    std::set<MeshID> ids;
    for (size_t i=0; i<d_meshes.size(); i++) {
        std::vector<MeshID> mesh_ids = d_meshes[i]->getLocalBaseMeshIDs();
        for (size_t j=0; j<mesh_ids.size(); j++)
            ids.insert(mesh_ids[j]);
    }
    return std::vector<MeshID>(ids.begin(),ids.end());
}


/********************************************************
* Function to return the element given an ID            *
********************************************************/
MeshElement MultiMesh::getElement( const MeshElementID &elem_id ) const
{
    MeshID mesh_id = elem_id.meshID();
    for (size_t i=0; i<d_meshes.size(); i++) {
        std::vector<MeshID> ids = d_meshes[i]->getLocalBaseMeshIDs();
        bool mesh_found = false;
        for (size_t j=0; j<ids.size(); j++) {
            if ( ids[j]==mesh_id )
                mesh_found = true;
        }
        if ( mesh_found )
            return d_meshes[i]->getElement( elem_id );
    }
    AMP_ERROR("A mesh matching the element's mesh id was not found");
    return MeshElement();
}


/********************************************************
* Function to return the mesh with the given ID         *
********************************************************/
boost::shared_ptr<Mesh>  MultiMesh::Subset( MeshID meshID ) const
{
    if ( d_meshID==meshID )
        return boost::const_pointer_cast<Mesh>( shared_from_this() );
    for (size_t i=0; i<d_meshes.size(); i++) {
        boost::shared_ptr<Mesh> mesh = d_meshes[i]->Subset(meshID);
        if ( mesh.get()!=NULL )
            return mesh;
    }
    return boost::shared_ptr<Mesh>();
}


/********************************************************
* Function to subset a mesh using a mesh iterator       *
********************************************************/
boost::shared_ptr<Mesh> MultiMesh::Subset( const MeshIterator &iterator_in ) const
{
    if ( iterator_in.size()==0 )
        return boost::shared_ptr<Mesh>();
    GeomType type = null;
    AMP_ASSERT(iterator_in.size()>0);
    MeshIterator iterator = iterator_in.begin();
    for (size_t i=0; i<iterator.size(); i++) {
        if ( type==null )
            type = iterator->elementType();
        if ( type!= iterator->elementType() )
            AMP_ERROR("Subset mesh requires all of the elements to be the same type");
        ++iterator;
    }
    // Subset for the name in each submesh
    std::vector<Mesh::shared_ptr> subset;
    std::set<MeshID> subsetID;
    for (size_t i=0; i<d_meshes.size(); i++) {
        iterator = Mesh::getIterator( Intersection, iterator_in, d_meshes[i]->getIterator(type,d_meshes[i]->getMaxGhostWidth()) );
        if ( iterator.size() == 0 ) 
            continue;
        boost::shared_ptr<Mesh> mesh = d_meshes[i]->Subset( iterator );
        if ( mesh.get()!=NULL ) {
            subset.push_back( mesh );
            subsetID.insert( mesh->meshID() );
        }
    }
    // Count the number of globally unique sub-meshes
    d_comm.setGather( subsetID ); 
    if ( subsetID.size() <= 1 ) {
        if ( subset.size() == 0 ) {
            return boost::shared_ptr<Mesh>();
        } else {
            boost::shared_ptr<Mesh> subsetMultiMesh( new MultiMesh( subset[0]->getComm(), subset ) );
            subsetMultiMesh->setName( d_name+"_subset" );
            return subsetMultiMesh;
        }
    }
    // Create a new multi-mesh to contain the subset
    int color = subset.size()==0 ? -1:0;
    AMP::AMP_MPI new_comm = d_comm.split( color );
    if ( new_comm.isNull() )
        return boost::shared_ptr<Mesh>();
    boost::shared_ptr<Mesh> subsetMultiMesh( new MultiMesh( new_comm, subset ) );
    subsetMultiMesh->setName( d_name+"_subset" );
    return subsetMultiMesh;
}


/********************************************************
* Function to return the mesh with the given name       *
********************************************************/
boost::shared_ptr<Mesh>  MultiMesh::Subset( std::string name ) const
{
    if ( d_name==name )
        return boost::const_pointer_cast<Mesh>( shared_from_this() );
    // Subset for the name in each submesh
    std::vector<Mesh::shared_ptr> subset;
    std::set<MeshID> subsetID;
    for (size_t i=0; i<d_meshes.size(); i++) {
        Mesh::shared_ptr mesh = d_meshes[i]->Subset(name);
        if ( mesh.get()!=NULL ) {
            subset.push_back( mesh ); 
            subsetID.insert( mesh->meshID() ); 
        }
    }
    // Count the number of globally unique sub-meshes
    d_comm.setGather( subsetID ); 
    if ( subsetID.size() <= 1 ) {
        if ( subset.size() == 0 ) {
            return boost::shared_ptr<Mesh>();
        } else {
            boost::shared_ptr<Mesh> subsetMultiMesh( new MultiMesh( subset[0]->getComm(), subset ) );
            subsetMultiMesh->setName( d_name+"_subset" );
            return subsetMultiMesh;
        }
    }
    // Create a new multi-mesh to contain the subset
    int color = subset.size()==0 ? -1:0;
    AMP::AMP_MPI new_comm = d_comm.split( color );
    if ( new_comm.isNull() )
        return boost::shared_ptr<Mesh>();
    boost::shared_ptr<MultiMesh> subsetMultiMesh( new MultiMesh( new_comm, subset ) );
    subsetMultiMesh->setName( name );
    return subsetMultiMesh;
}


/********************************************************
* Displace a mesh                                       *
********************************************************/
void MultiMesh::displaceMesh( std::vector<double> x_in )
{
    // Check x
    AMP_INSIST((short int)x_in.size()==PhysicalDim,"Displacement vector size should match PhysicalDim");
    std::vector<double> x = x_in;
    d_comm.minReduce(&x[0],(int)x.size());
    for (size_t i=0; i<x.size(); i++)
        AMP_INSIST(fabs(x[i]-x_in[i])<1e-12,"x does not match on all processors");
    // Displace the meshes
    for (size_t i=0; i<d_meshes.size(); i++)
        d_meshes[i]->displaceMesh(x);
    // Update the bounding box
    for (int i=0; i<PhysicalDim; i++) {
        d_box[2*i+0] += x[i];
        d_box[2*i+1] += x[i];
        d_box_local[2*i+0] += x[i];
        d_box_local[2*i+1] += x[i];
    }
}
#ifdef USE_AMP_VECTORS
void MultiMesh::displaceMesh( const AMP::LinearAlgebra::Vector::const_shared_ptr x )
{
    // Displace the individual meshes
    for (size_t i=0; i<d_meshes.size(); i++)
        d_meshes[i]->displaceMesh(x);
    // Compute the bounding box of the multimesh
    d_box_local = d_meshes[0]->getBoundingBox();
    for (size_t i=1; i<d_meshes.size(); i++) {
        std::vector<double> meshBox = d_meshes[i]->getBoundingBox();
        for (int j=0; j<PhysicalDim; j++) {
            if ( meshBox[2*j+0] < d_box_local[2*j+0] ) { d_box_local[2*j+0] = meshBox[2*j+0]; }
            if ( meshBox[2*j+1] > d_box_local[2*j+1] ) { d_box_local[2*j+1] = meshBox[2*j+1]; }
        }
    }
    d_box = std::vector<double>(PhysicalDim*2);
    for (int i=0; i<PhysicalDim; i++) {
        d_box[2*i+0] = d_comm.minReduce( d_box_local[2*i+0] );
        d_box[2*i+1] = d_comm.maxReduce( d_box_local[2*i+1] );
    } 
}
#endif


/********************************************************
* Function to copy a key from database 1 to database 2  *
* If the key is an array of size N, it will only copy   *
* the ith value.                                        *
********************************************************/
static void copyKey( AMP::Database::shared_ptr database1, 
    std::vector<AMP::Database::shared_ptr> database2, const std::string& key, bool select,
    const std::string& iterator, const std::vector<std::string>& index ) 
{
    int size = database1->getArraySize(key);
    AMP::Database::DataType type = database1->getArrayType(key);
    int N = (int) database2.size();
    switch (type) {
        case AMP::Database::AMP_INVALID: {
            // We don't know what this is
            AMP_ERROR("Invalid database object");
            } break;
        case AMP::Database::AMP_DATABASE: {
            // Copy the database
            boost::shared_ptr<AMP::Database> subDatabase1 = database1->getDatabase(key);
            for (size_t i=0; i<database2.size(); i++) {
                std::vector<AMP::Database::shared_ptr> subDatabase2(1,database2[i]->putDatabase(key));
                std::vector<std::string> index2(1,index[i]);
                std::vector<std::string> subKeys = subDatabase1->getAllKeys();
                for (size_t j=0; j<subKeys.size(); j++)
                    copyKey( subDatabase1, subDatabase2, subKeys[j], false, iterator, index2 );
            }
            } break;
        case AMP::Database::AMP_BOOL: {
            // Copy a bool
            std::vector<unsigned char> data = database1->getBoolArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            for (size_t i=0; i<database2.size(); i++) {
                if ( N==size && select )
                    database2[i]->putBool(key,data[i]!=0);
                else
                    database2[i]->putBoolArray(key,data);
            }
            } break;
        case AMP::Database::AMP_CHAR: {
            // Copy a char
            std::vector<char> data = database1->getCharArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            for (size_t i=0; i<database2.size(); i++) {
                if ( N==size && select ) {
                    database2[i]->putChar(key,data[i]);
                } else {
                    // We need to try a search and replace
                    database2[i]->putCharArray(key,data);
                }
            }
            } break;
        case AMP::Database::AMP_INT: {
            // Copy an int
            std::vector<int> data = database1->getIntegerArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            for (size_t i=0; i<database2.size(); i++) {
                if ( N==size && select )
                    database2[i]->putInteger(key,data[i]);
                else
                    database2[i]->putIntegerArray(key,data);
            }
            } break;
        case AMP::Database::AMP_COMPLEX: {
            // Copy a complex number
            std::vector< std::complex<double> > data = database1->getComplexArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            for (size_t i=0; i<database2.size(); i++) {
                if ( N==size && select )
                    database2[i]->putComplex(key,data[i]);
                else
                    database2[i]->putComplexArray(key,data);
            }
            } break;
        case AMP::Database::AMP_DOUBLE: {
            // Copy a double number
            std::vector<double> data = database1->getDoubleArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            for (size_t i=0; i<database2.size(); i++) {
                if ( N==size && select )
                    database2[i]->putDouble(key,data[i]);
                else
                    database2[i]->putDoubleArray(key,data);
            }
            } break;
        case AMP::Database::AMP_FLOAT: {
            // Copy a float
            std::vector<float> data = database1->getFloatArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            for (size_t i=0; i<database2.size(); i++) {
                if ( N==size && select )
                    database2[i]->putFloat(key,data[i]);
                else
                    database2[i]->putFloatArray(key,data);
            }
            } break;
        case AMP::Database::AMP_STRING: {
            // Copy a string
            std::vector<std::string> data = database1->getStringArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            for (size_t i=0; i<database2.size(); i++) {
                std::vector<std::string> data2;
                if ( N==size && select ) {
                    data2.resize(1);
                    data2[0] = data[i];
                } else {
                    data2 = data;
                }
                if ( !iterator.empty() ) {
                    for (size_t j=0; j<data2.size(); j++) {
                        replaceSubString(data2[j],iterator,index[i]);
                    }
                }
                database2[i]->putStringArray(key,data2);
            }
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
static void replaceSubString( std::string& string, const std::string& search, const std::string& replace )
{
    while ( 1 ) {
        size_t pos = string.find(search);
        if ( pos == std::string::npos ) 
            break;
        string.replace( pos, search.size(), replace );
    }
}
static void replaceText( boost::shared_ptr<AMP::Database>& database, const std::string& search, const std::string& replace )
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
std::vector<AMP_MPI> MultiMesh::createComms( const std::vector<rank_list> &groups )
{
    int myRank = d_comm.getRank();
    std::vector<AMP_MPI> comms(groups.size());
    for (size_t i=0; i<groups.size(); i++) {
        int color = -1;
        for (size_t j=0; j<groups[i].size(); j++) {
            if ( groups[i][j]==myRank )
                color = 0;
        }
        comms[i] = d_comm.split(color,d_comm.getRank());
        if ( color!=-1 ) 
            AMP_ASSERT(comms[i].getSize()==(int)groups[i].size());
    }
    return comms;
}


/********************************************************
* Function to perform the load balance                  *
********************************************************/
std::vector<MultiMesh::rank_list> MultiMesh::loadBalancer( int N_procs, 
    const std::vector<MeshParameters::shared_ptr> &meshParams, const std::vector<size_t> &meshSizes, int method )
{
    // Deal with the special cases directly
    if ( meshParams.size()<=1 || N_procs==1 || method==0 ) {
        // Everybody is on the same communicator
        std::vector<int> allRanks(N_procs);
        for (int i=0; i<N_procs; i++)
            allRanks[i] = i;
        return std::vector<rank_list>(meshParams.size(),allRanks);
    } else if ( method == 1 ) {
        // We want to split the meshes onto independent processors
        std::vector<comm_groups> groups;
        if ( N_procs > (int)meshParams.size() ) {
            // We have more processors than meshes
            groups = independentGroups1( N_procs, meshParams, meshSizes );
        } else {
            // The number of meshes is <= the number of processors
            // Create the weights
            size_t totalMeshSize = 0;
            for (size_t i=0; i<meshSizes.size(); i++)
                totalMeshSize += meshSizes[i];
            std::vector<double> weights(meshParams.size(),0);
            for (size_t i=0; i<meshSizes.size(); i++)
                weights[i] = ((double)meshSizes[i])/((double)totalMeshSize);
            // Create the ids and recursively perform the load balance
            std::vector<std::pair<double,int> >  ids(weights.size());
            for (size_t i=0; i<weights.size(); i++)
                ids[i] = std::pair<double,int>( weights[i], (int) i );
            groups = independentGroups2( N_procs, ids );
        }
        // Split the comm into the appropriate groups
        std::vector<rank_list> commGroups(meshParams.size());
        int myStartRank=0;
        for (size_t i=0; i<groups.size(); i++) {
            for (int myRank=myStartRank; myRank<myStartRank+groups[i].N_procs; myRank++) {
                for (size_t j=0; j<groups[i].ids.size(); j++)
                    commGroups[groups[i].ids[j]].push_back( myRank );
            }
            myStartRank += groups[i].N_procs;
        }
        for (size_t i=0; i<commGroups.size(); i++)
            AMP_ASSERT(commGroups[i].size()>0);
        return commGroups;
    } else if ( method == 2 ) {
        // We want to try to achieve a more balanced approach to the load balance
        AMP_ERROR("This load balancer is not implimented yet");
    } else {
        AMP_ERROR("Unknown load balancer");
    }
    return std::vector<rank_list>(0);
}
void MultiMesh::addProcSimulation( const LoadBalance& mesh, std::vector<LoadBalance> &submeshes, int rank, char &decomp )
{
    boost::shared_ptr<MultiMeshParameters> multimeshParams = boost::dynamic_pointer_cast<MultiMeshParameters>( mesh.getParams() );
    AMP_ASSERT(multimeshParams.get());
    AMP_ASSERT(submeshes.size()==multimeshParams->params.size());
    if ( multimeshParams->method==1 ) {
        if ( mesh.getRanks().size()==submeshes.size() ) {
            // Special case where the domain decomposition changes
            std::vector<int> rank2(1,0);
            for (size_t i=0; i<submeshes.size(); i++) {
                rank2[0] = i;
                submeshes[i].changeRanks( rank2 );
            }
            decomp = 1;
        } else if ( mesh.getRanks().size()<submeshes.size() ) {
            // We need to create new group sets
            const std::vector<int>& ranks = mesh.getRanks();
            int N_procs = (int) ranks.size();
            // Create the weights
            std::vector<double> weights(submeshes.size(),0);
            for (size_t i=0; i<submeshes.size(); i++)
                weights[i] = ((double)submeshes[i].getSize())/((double)mesh.getSize());
            // Create the ids and recursively perform the load balance
            std::vector<std::pair<double,int> >  ids(weights.size());
            for (size_t i=0; i<weights.size(); i++)
                ids[i] = std::pair<double,int>( weights[i], (int) i );
            std::vector<comm_groups> groups = independentGroups2( N_procs, ids );
            // Create the comms for each mesh
            std::vector<rank_list> comms(submeshes.size(),rank_list(1,0));
            AMP_ASSERT((int)groups.size()==N_procs);
            for (size_t i=0; i<groups.size(); i++) {
                AMP_ASSERT(groups[i].N_procs==1);
                for (size_t j=0; j<groups[i].ids.size(); j++)
                    comms[groups[i].ids[j]][0] = ranks[i];
            }
            // Create the new submeshes
            for (size_t i=0; i<submeshes.size(); i++)
                submeshes[i].changeRanks( comms[i] );
        } else {
            // We need to add a processor to the mesh with the largest # of elements per processor
            size_t max = 0;
            int i_max = 0;
            std::vector<size_t> max_list(submeshes.size(),0);
            for (size_t i=0; i<submeshes.size(); i++) {
                max_list[i] = submeshes[i].max();
                if ( max_list[i] > max ) {
                    max = max_list[i];
                    i_max = i;
                }
            }
            submeshes[i_max].addProc( rank );
        }
    } else {
        AMP_ERROR("Not ready for other load balancers yet");
    }
}



/********************************************************
* Functions to solve simple load balancing computations *
********************************************************/
std::vector<MultiMesh::comm_groups>  MultiMesh::independentGroups1(
    int N_procs, const std::vector<MeshParameters::shared_ptr> &meshParameters, const std::vector<size_t> &size )
{
    // This will distribute the groups onto the processors such that no groups share any processors
    // Note: this requires the number of processors to be >= the number of ids
    AMP_ASSERT(N_procs>=(int)meshParameters.size());
    AMP_ASSERT(size.size()==meshParameters.size());
    // Handle the special case if the # of processors == the # of ids
    if ( N_procs == (int)meshParameters.size() ) {
        std::vector<comm_groups> groups(meshParameters.size());
        for (int i=0; i<N_procs; i++) {
            groups[i].N_procs = 1;
            groups[i].ids = std::vector<int>(1,i);
        }
        return groups;
    }
    // Start by using ~80% of the procs
    size_t N_total = 0;
    for (size_t i=0; i<meshParameters.size(); i++)
        N_total += size[i];
    std::vector<LoadBalance> load_balance(meshParameters.size());
    std::vector<size_t> max_size(meshParameters.size());
    std::vector<int> rank1(N_procs);
    int N_procs_remaining = N_procs;
    for (size_t i=0; i<meshParameters.size(); i++) {
        int N_proc_local = floor(0.8*((double)N_procs)*((double)size[i])/((double)N_total))-1;
        N_proc_local = std::min(N_proc_local,N_procs-(int)size.size());
        N_proc_local = std::max(N_proc_local,1);
        rank1.resize(N_proc_local);
        for (int j=0; j<N_proc_local; j++)
            rank1[j] = j;
        load_balance[i] = LoadBalance( meshParameters[i], rank1, size[i] );
        max_size[i] = load_balance[i].max();
        N_procs_remaining -= N_proc_local;
    }
    // While we have processors remaining, try to add them so that we minimize 
    // the maximum number of elements on any rank
    while ( N_procs_remaining > 0 ) {
        // Find the two meshes with the greatest number of elements
        size_t max1 = 0;
        size_t max2 = 0;
        size_t i1 = 0;
        size_t i2 = 0;
        for (size_t i=0; i<load_balance.size(); i++) {
            size_t local_max = max_size[i];
            if ( local_max > max1 ) {
                max2 = max1;
                i2 = i1;
                max1 = local_max;
                i1 = i;
            } else if ( local_max > max2 ) {
                max2 = local_max;
                i2 = i;
            }
        }
        if ( max1==max2 )
            max2--;
        // Add processors to the mesh pointed to by i1, until max1 is <= max2
        AMP_ASSERT( max1>max2 && max1>0 );
        while ( N_procs_remaining>0 && max1>max2 ) {
            N_procs_remaining--;
            load_balance[i1].addProc( (int) load_balance[i1].getRanks().size() );
            max_size[i1] = load_balance[i1].max();
            max1 = max_size[i1];
        }
    }
    // We have filled all processors, create the groups
    AMP_ASSERT(N_procs_remaining==0);
    std::vector<comm_groups> groups(meshParameters.size());
    for (size_t i=0; i<meshParameters.size(); i++) {
        groups[i].N_procs = (int) load_balance[i].getRanks().size();
        groups[i].ids = std::vector<int>(1,(int)i);
    }
    return groups;
}
std::vector<MultiMesh::comm_groups>  MultiMesh::independentGroups2(
    int N_procs, std::vector<std::pair<double,int> >  &ids_in )
{
    // This will distribute the groups onto the processors such that each group has exactly one processor.  
    // Note: this requires the number of processors to be <= the number of ids
    AMP_ASSERT(N_procs<=(int)ids_in.size());
    std::vector<comm_groups> groups;
    groups.reserve(N_procs);
    // Handle the special case if the # of processors == the # of ids
    if ( N_procs == (int)ids_in.size() ) {
        groups.resize(N_procs);
        for (int i=0; i<N_procs; i++) {
            groups[i].N_procs = 1;
            groups[i].ids = std::vector<int>(1,ids_in[i].second);
        }
        return groups;
    }
    // Sort by the weights
    std::vector<std::pair<double,int> > ids = ids_in;
    AMP::Utilities::quicksort( ids );
    comm_groups tmp_group;
    tmp_group.N_procs = 1;
    tmp_group.ids.reserve(8);
    if ( ids[0].first==ids[ids.size()-1].first ) {
        // Special case where all meshes are the same size  
        groups.resize(N_procs);
        for (size_t i=0; i<(size_t)N_procs; i++) {
            size_t i0 = (i*ids.size())/N_procs;
            size_t i1 = ((i+1)*ids.size())/N_procs;
            AMP_ASSERT(i1>i0&&i1<=ids.size());
            groups[i].N_procs = 1;
            groups[i].ids = std::vector<int>(i1-i0);
            for (size_t j=0; j<i1-i0; j++)
                groups[i].ids[j] = (int) (i0+j);
        }
        return groups;
    }
    // Remove any ids that require a processor by themselves
    for (size_t i=ids.size()-1; i>=0; i--) {
        double total_weight = 0.0;
        for (size_t j=0; j<ids.size(); j++)
            total_weight += ids[i].first;
        double weight_avg = total_weight / ((double) N_procs-groups.size());
        if ( ids[i].first < weight_avg )
            break;
        tmp_group.ids.resize(1);
        tmp_group.ids[0] = ids[i].second;
        groups.push_back(tmp_group);
        ids.resize(i);
    }
    while ( !ids.empty() ) {
        // Group meshes until we reach or exceed the average number of elements remaining
        if ( (ids.size()+groups.size()) == (size_t) N_procs ) {
            tmp_group.ids.resize(1);
            for (size_t i=0; i<ids.size(); i++) {
                tmp_group.ids = std::vector<int>(1,ids_in[i].second);
                groups.push_back( tmp_group );
            }
            ids.resize(0);
            break;
        }
        double total_weight = 0.0;
        for (size_t j=0; j<ids.size(); j++)
            total_weight += ids[j].first;
        double weight_avg = total_weight / ((double) N_procs-groups.size());
        tmp_group.ids.resize(0);
        double weight_sum = 0.0;
        while ( weight_sum < weight_avg ) {
            size_t i = ids.size()-1;
            if ( weight_sum+ids[i].first < weight_avg ) {
                weight_sum += ids[i].first;
                tmp_group.ids.push_back( ids[i].second );
                ids.resize(i);
            } else {
                for (size_t j=0; j<=i; j++) {
                    if ( weight_sum+ids[j].first >= weight_avg ) {
                        weight_sum += ids[j].first;
                        tmp_group.ids.push_back( ids[j].second );
                        for (size_t k=j; k<i; k++)
                            ids[k] = ids[k+1];
                        ids.resize(i);
                        break;
                    }
                }
            }
        }
        groups.push_back( tmp_group );
    }
    AMP_ASSERT((int)groups.size()==N_procs);
    AMP_ASSERT(ids.size()==0);
    return groups;
}


} // Mesh namespace
} // AMP namespace

