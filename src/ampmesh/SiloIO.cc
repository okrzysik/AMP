#include "SiloIO.h"
#include "ampmesh/MultiMesh.h"
#include "utils/ProfilerApp.h"

namespace AMP { 
namespace Mesh {


/************************************************************
* Constructor                                               *
************************************************************/
SiloIO::SiloIO( )
{
    d_comm = AMP_MPI(AMP_COMM_WORLD);
    d_dim = -1;
    decomposition = 0;
}


/************************************************************
* Some basic functions                                      *
************************************************************/
std::string SiloIO::getExtension() 
{ 
    return "silo"; 
}
void SiloIO::setDecomposition( int d )
{
    AMP_INSIST(d>=0&&d<=1,"decomposition must be 0 or 1");
    decomposition = d;
}


#ifdef USE_EXT_SILO

// Some internal functions
static void createSiloDirectory( DBfile *FileHandle, std::string path );


/************************************************************
* Function to read a silo file                              *
************************************************************/
void SiloIO::readFile( const std::string &fname )
{ 
    AMP_ERROR("readFile is not implimented yet");
}


/************************************************************
* Function to write a silo file                             *
* Note: it appears that only one prcoessor may write to a   *
* file at a time, and that once a processor closes the file *
* it cannot reopen it (or at least doing this on the        *
* processor that created the file creates problems).        *
************************************************************/
void SiloIO::writeFile( const std::string &fname_in, size_t iteration_count )
{ 
    PROFILE_START("writeFile");
    // Create the file name
    std::stringstream tmp;
    tmp << fname_in << "_" << iteration_count << "." << getExtension();
    std::string fname = tmp.str();
    // Check that the dimension is matched across all processors
    PROFILE_START("sync dim");
    int dim2 = d_comm.maxReduce(d_dim);
    if ( d_dim == -1 )
        d_dim = dim2;
    else
        AMP_ASSERT(dim2==d_dim);
    d_comm.barrier();
    PROFILE_STOP("sync dim");
    // Syncronize all vectors
    #ifdef USE_AMP_VECTORS
    PROFILE_START("makeConsistent");
    for (size_t i=0; i<d_vectors.size(); ++i) {
        AMP::LinearAlgebra::Vector::UpdateState localState = d_vectors[i]->getUpdateStatus();
        if ( localState==AMP::LinearAlgebra::Vector::ADDING ) 
            d_vectors[i]->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_ADD );
        else
            d_vectors[i]->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    }
    PROFILE_STOP("makeConsistent");
    #endif
    // Write the data for each base mesh
    if ( decomposition==0 ) {
        // Write all mesh data to the main file
        for (int i=0; i<d_comm.getSize(); ++i) {
            if ( d_comm.getRank()==i ) {
                // Open the file
                DBfile *FileHandle;
                if ( d_comm.getRank()==0 ) {
                    FileHandle = DBCreate( fname.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_HDF5 );
                } else {
                    FileHandle = DBOpen ( fname.c_str(), DB_HDF5, DB_APPEND );
                }
                // Write the base meshes
                std::map<AMP::Mesh::MeshID,siloBaseMeshData>::iterator iterator;
                for (iterator=d_baseMeshes.begin(); iterator!=d_baseMeshes.end(); ++iterator) {
                    siloBaseMeshData &data = iterator->second;
                    data.file = fname.c_str();
                    AMP_ASSERT(data.id==iterator->first);
                    writeMesh( FileHandle, iterator->second );
                }
                // Close the file
                DBClose ( FileHandle );
            }
            d_comm.barrier();
        }
    } else if ( decomposition==1 ) {
        // Every rank will write a seperate file
        std::stringstream tmp2;
        tmp2 << fname_in << "_" << iteration_count << "." << d_comm.getRank()+1 << "." << getExtension();
        std::string fname_rank = tmp2.str();
        DBfile *FileHandle = DBCreate( fname_rank.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_HDF5 );
        // Write the base meshes
        std::map<AMP::Mesh::MeshID,siloBaseMeshData>::iterator iterator;
        for (iterator=d_baseMeshes.begin(); iterator!=d_baseMeshes.end(); ++iterator) {
            siloBaseMeshData &data = iterator->second;
            data.file = fname_rank.c_str();
            AMP_ASSERT(data.id==iterator->first);
            writeMesh( FileHandle, iterator->second );
        }
        // Close the file
        DBClose ( FileHandle );
    } else {
        AMP_ERROR("Unknown file decomposition");
    }
    // Write the summary results (multimeshes, multivariables, etc.)
    if ( decomposition!=0 ) {
        if ( d_comm.getRank()==0 ) {
            DBfile *FileHandle = DBCreate( fname.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_HDF5 );
            DBClose ( FileHandle );
        }
        d_comm.barrier();
    }
    writeSummary( fname );
    PROFILE_STOP("writeFile");
}


/************************************************************
* Function to register a mesh with silo                     *
************************************************************/
void SiloIO::registerMesh( AMP::Mesh::Mesh::shared_ptr mesh, int level, std::string path )
{ 
    if ( d_dim == -1 )
        d_dim = mesh->getDim();
    else
        AMP_INSIST(d_dim==mesh->getDim(),"All meshes must have the same number of physical dimensions");
    AMP_INSIST(level>=0&&level<=3,"Invalid value for level");
    boost::shared_ptr<AMP::Mesh::MultiMesh> multimesh = boost::dynamic_pointer_cast<AMP::Mesh::MultiMesh>( mesh );
    if ( multimesh.get()==NULL ) {
        // We are dealing with a single mesh
        siloBaseMeshData data;
        data.id = mesh->meshID();
        data.mesh = mesh;
        data.rank = mesh->getComm().getRank()+1;
        std::stringstream stream;
        stream << data.rank;
        std::string rank = stream.str();
        data.ownerRank = d_comm.getRank(); // Everybody owns the mesh because every rank is an independent mesh
        data.meshName = "rank_" + rank;
        data.path = path + mesh->getName() + "_/";
        if ( d_baseMeshes.find(mesh->meshID()) == d_baseMeshes.end() )
            d_baseMeshes.insert( std::pair<AMP::Mesh::MeshID,siloBaseMeshData>(mesh->meshID(),data) );
        // Create and register a multimesh for the current mesh
        if ( level>0 ) {
            siloMultiMeshData data2;
            data2.id = mesh->meshID();
            data2.mesh = mesh;
            data2.name = path + mesh->getName();
            data.ownerRank = mesh->getComm().bcast(d_comm.getRank(),0);
            d_multiMeshes.insert( std::pair<AMP::Mesh::MeshID,siloMultiMeshData>(data2.id,data2) );
        }
        // Create and register a multimesh for the rank
        if ( level==3 ) {
            // Create a unique id for each rank
            AMP::Mesh::uint64 tmp_id = mesh->meshID().getData();
            AMP::Mesh::uint64 root2 = d_comm.getRank()+1;
            tmp_id = (root2<<48) + tmp_id;
            siloMultiMeshData data2;
            data2.id = AMP::Mesh::MeshID(tmp_id);
            data2.mesh = mesh;
            data2.name = path + mesh->getName() + "_/rank_" + rank;
            data.ownerRank = d_comm.getRank();
            d_multiMeshes.insert( std::pair<AMP::Mesh::MeshID,siloMultiMeshData>(data2.id,data2) );
        }
    } else {
        // We are dealining with a multimesh, register the current mesh and sub meshes
        int level2 = level;
        if ( level == 1 ) { level2 = 0; }
        std::string new_path = path + mesh->getName() + "_/";
        std::vector<AMP::Mesh::Mesh::shared_ptr> submeshes = multimesh->getMeshes();
        for (size_t i=0; i<submeshes.size(); ++i)
            registerMesh( submeshes[i], level2, new_path );
        if ( level > 0 ) {
            siloMultiMeshData data;
            data.id = mesh->meshID();
            data.mesh = mesh;
            data.name = path + mesh->getName();
            data.ownerRank = mesh->getComm().bcast(d_comm.getRank(),0);
            d_multiMeshes.insert( std::pair<AMP::Mesh::MeshID,siloMultiMeshData>(mesh->meshID(),data) );
        }
    }
}

/************************************************************
* Function to get the mesh ids to use for registering       *
************************************************************/
std::vector<AMP::Mesh::MeshID> SiloIO::getMeshIDs( AMP::Mesh::Mesh::shared_ptr mesh )
{
    std::vector<AMP::Mesh::MeshID> ids;
    boost::shared_ptr<AMP::Mesh::MultiMesh> multimesh = boost::dynamic_pointer_cast<AMP::Mesh::MultiMesh>( mesh );
    if ( multimesh.get()==NULL ) {
        // We are dealing with a single mesh
        ids = std::vector<AMP::Mesh::MeshID>(1,mesh->meshID());
    } else {
        // We are dealining with a multimesh
        std::vector<AMP::Mesh::Mesh::shared_ptr> meshes = multimesh->getMeshes();
        for (size_t i=0; i<meshes.size(); ++i) {
            std::vector<AMP::Mesh::MeshID> ids2 = getMeshIDs( meshes[i] );
            for (size_t j=0; j<ids2.size(); ++j)
                ids.push_back( ids2[j] );
        }
    }
    return ids;
}


/************************************************************
* Function to register a vector with silo                   *
************************************************************/
#ifdef USE_AMP_VECTORS
void SiloIO::registerVector( AMP::LinearAlgebra::Vector::shared_ptr vec, 
    AMP::Mesh::Mesh::shared_ptr mesh, AMP::Mesh::GeomType type, const std::string &name_in )
{ 
    // Make sure the mesh has been registered
    registerMesh( mesh );
    // Perform some error checking
    AMP::Discretization::DOFManager::shared_ptr DOFs = vec->getDOFManager();
    AMP::Mesh::MeshIterator iterator1 = mesh->getIterator( type, 0 );
    AMP::Mesh::MeshIterator iterator2 = DOFs->getIterator( );
    AMP::Mesh::MeshIterator iterator3 = AMP::Mesh::Mesh::getIterator( AMP::Mesh::Intersection, iterator1, iterator2 );
    if ( iterator1.size() != iterator3.size() )
        AMP_ERROR("vector does not cover the entire mesh for the given entity type");
    std::vector<size_t> dofs;
    DOFs->getDOFs( iterator1->globalID(), dofs );
    int DOFsPerPoint = dofs.size();
    if ( type==AMP::Mesh::Vertex )
        iterator1 = mesh->getIterator( type, 1 );
    for (size_t i=0; i<iterator1.size(); ++i) {
        DOFs->getDOFs( iterator1->globalID(), dofs );
        AMP_ASSERT((int)dofs.size()==DOFsPerPoint);
        ++iterator1;
    }
    // Register the vector with the appropriate base meshes
    std::string name = vec->type();
    if ( !name.empty() )
        name = name_in;
    std::vector<AMP::Mesh::MeshID> ids = getMeshIDs( mesh );
    for (size_t i=0; i<ids.size(); ++i) {
        std::map<AMP::Mesh::MeshID,siloBaseMeshData>::iterator it = d_baseMeshes.find(ids[i]);
        if ( it == d_baseMeshes.end() )
            continue;
        it->second.varName.push_back( name );
        it->second.varType.push_back( type );
        it->second.varSize.push_back( DOFsPerPoint );
        it->second.vec.push_back( vec );
    }
    // Register the vector with the appropriate multi-meshes
    std::map<AMP::Mesh::MeshID,siloMultiMeshData>::iterator it = d_multiMeshes.find(mesh->meshID());
    AMP_ASSERT(it!=d_multiMeshes.end());
    it->second.varName.push_back( name );
    // Add the vector to the list of vectors so we can perform makeConsistent
    d_vectors.push_back(vec);
    // Add the variable name to the list of variables
    d_varNames.insert(name);
}
#endif


/************************************************************
* Function to write a mesh                                  *
************************************************************/
void SiloIO::writeMesh( DBfile *FileHandle, const siloBaseMeshData &data )
{ 
    PROFILE_START("writeMesh");
    AMP::Mesh::Mesh::shared_ptr mesh = data.mesh;
    // Get the zone (element) lists
    AMP::Mesh::MeshIterator elem_iterator = mesh->getIterator(mesh->getGeomType(),0);
    AMP_ASSERT(elem_iterator.size()>0);
    std::vector<AMP::Mesh::MeshElement> nodes = elem_iterator->getElements(AMP::Mesh::Vertex);
    int shapesize = nodes.size();
    int shapetype;
    if ( shapesize == 8 )
        shapetype = DB_ZONETYPE_HEX;
    else if ( shapesize == 4 )
        shapetype = DB_ZONETYPE_QUAD;
    else
        AMP_ERROR("Unknown element type");
    int shapecnt = elem_iterator.size();
    // Get the node list (unique integer for each node) and coordinates
    AMP::Mesh::MeshIterator node_iterator = mesh->getIterator(AMP::Mesh::Vertex,1);
    std::vector<AMP::Mesh::MeshElementID> nodelist_ids(node_iterator.size());
    for (size_t i=0; i<node_iterator.size(); ++i) {
        nodelist_ids[i] = node_iterator->globalID();
        ++node_iterator;
    }
    AMP::Utilities::quicksort(nodelist_ids);
    double *coord[3];
    for (int i=0; i<d_dim; ++i)
        coord[i] = new double[node_iterator.size()];
    node_iterator = mesh->getIterator(AMP::Mesh::Vertex,1);
    for (size_t i=0; i<node_iterator.size(); ++i) {
        size_t index = AMP::Utilities::findfirst( nodelist_ids, node_iterator->globalID() );
        AMP_ASSERT(nodelist_ids[index]==node_iterator->globalID());
        std::vector<double> elem_coord = node_iterator->coord();
        for (int j=0; j<d_dim; ++j)
            coord[j][index] = elem_coord[j];
        ++node_iterator;
    }
    elem_iterator = mesh->getIterator(mesh->getGeomType(),0);
    std::vector<int> nodelist;
    nodelist.reserve(shapesize*elem_iterator.size());
    for (size_t i=0; i<elem_iterator.size(); ++i) {
        nodes = elem_iterator->getElements(AMP::Mesh::Vertex);
        AMP_INSIST((int)nodes.size()==shapesize,"Mixed element types is currently not supported");
        for (size_t j=0; j<nodes.size(); ++j) {
            size_t index = AMP::Utilities::findfirst( nodelist_ids, nodes[j].globalID() );
            AMP_ASSERT(nodelist_ids[index]==nodes[j].globalID());
            nodelist.push_back( (int) index );
        }
        ++elem_iterator;
    }
    // Create the directory for the mesh
    std::string tmp_path = data.path;
    while ( tmp_path.size() > 0 ) {
        if ( tmp_path[0]=='/' ) {
            tmp_path.erase(0,1);
            continue;
        }
        size_t pos = tmp_path.find_first_of('/');
        if ( pos==std::string::npos ) { pos = tmp_path.size(); }
        std::string subdir = tmp_path.substr(0,pos);
        DBtoc *toc = DBGetToc( FileHandle );
        bool subdir_found = false;
        for (int i=0; i<toc->ndir; ++i) {
            if ( subdir.compare(toc->dir_names[i])==0 )
                subdir_found = true;
        }
        if ( !subdir_found )
            DBMkDir( FileHandle, subdir.c_str() );
        DBSetDir( FileHandle, subdir.c_str() );
        tmp_path.erase(0,pos);
    }
    DBSetDir( FileHandle, "/" );
    DBSetDir( FileHandle, data.path.c_str() );
    // Write the elements (connectivity)
    std::stringstream  stream;
    stream << data.rank;
    std::string rank = stream.str();
    std::string meshName = data.meshName;
    std::string zoneName = "zone_" + rank;
    AMP::Mesh::MeshIterator  element_iterator = mesh->getIterator(mesh->getGeomType(),0);
    int num_elems = (int) element_iterator.size();
    DBPutZonelist2( FileHandle, zoneName.c_str(), num_elems, d_dim, 
        &nodelist[0], nodelist.size(), 0, 0, 0, &shapetype, &shapesize, &shapecnt, 1, 0 );
    // Write the mesh
    DBPutUcdmesh( FileHandle, meshName.c_str(), d_dim,
        NULL, coord, node_iterator.size(), nodelist.size(),
        zoneName.c_str(), 0, DB_DOUBLE, 0 );
    for (int i=0; i<d_dim; ++i)
        delete [] coord[i];
    // Write the variables
    #ifdef USE_AMP_VECTORS
        for (size_t i=0; i<data.varName.size(); ++i) {
            AMP::Discretization::DOFManager::shared_ptr DOFs = data.vec[i]->getDOFManager();
            int nvar = 0;
            int centering = 0;
            double **var = new double*[data.varSize[i]];
            for (int j=0; j<data.varSize[i]; ++j)
                var[j] = NULL;
            const char *varnames[] = {"1","2","3"};
            if ( data.varType[i]==AMP::Mesh::Vertex ) {
                // We are saving node-centered data
                centering = DB_NODECENT;
                nvar = (int) nodelist_ids.size();
                for (int j=0; j<data.varSize[i]; ++j)
                    var[j] = new double[nvar];
                std::vector<size_t> dofs(data.varSize[i]);
                std::vector<double> vals(data.varSize[i]);
                for (int j=0; j<nvar; ++j) {
                    DOFs->getDOFs( nodelist_ids[j], dofs );
                    AMP_ASSERT((int)dofs.size()==data.varSize[i]);
                    data.vec[i]->getValuesByGlobalID ( data.varSize[i], &dofs[0], &vals[0] );
                    for (int k=0; k<data.varSize[i]; ++k)
                        var[k][j] = vals[k];
                }
            } else if ( data.varType[i]==mesh->getGeomType() ) {
                // We are saving cell-centered data
                centering = DB_ZONECENT;
                nvar = (int) num_elems;
                for (int j=0; j<data.varSize[i]; ++j)
                    var[j] = new double[nvar];
                std::vector<size_t> dofs(data.varSize[i]);
                std::vector<double> vals(data.varSize[i]);
                AMP::Mesh::MeshIterator it = element_iterator.begin();
                for (int j=0; j<nvar; ++j) {
                    DOFs->getDOFs( it->globalID(), dofs );
                    data.vec[i]->getValuesByGlobalID ( data.varSize[i], &dofs[0], &vals[0] );
                    for (int k=0; k<data.varSize[i]; ++k)
                        var[k][j] = vals[k];
                    ++it;
                }
            } else {
                // We are storing edge or face data
                AMP_ERROR("The silo writer currently only supports Vertex and Cell data");
            }
            std::string varNameRank = data.varName[i]+"P"+rank;
            if ( data.varSize[i]==1 || data.varSize[i]==d_dim || data.varSize[i]==d_dim*d_dim ) {
                // We are writing a scalar, vector, or tensor variable
                DBPutUcdvar( FileHandle, varNameRank.c_str(), meshName.c_str(), 
                    data.varSize[i], (char**) varnames, var, nvar, NULL, 0, DB_DOUBLE, centering, NULL);
            } else {
                // Write each component
                for (int j=0; j<data.varSize[i]; ++j) {
                    std::stringstream  stream;
                    stream << varNameRank << "_" << j;
                    DBPutUcdvar( FileHandle, stream.str().c_str(), meshName.c_str(), 
                        1, (char**) varnames, &var[j], nvar, NULL, 0, DB_DOUBLE, centering, NULL);
                }
            }
            for (int j=0; j<data.varSize[i]; ++j) {
                if ( var[j] != NULL )
                    delete [] var[j];
            }
            delete [] var;
        }
    #endif
    // Change the directory back to root
    DBSetDir( FileHandle, "/" );
    PROFILE_STOP("writeMesh");
}


/************************************************************
* Function to syncronize the multimesh data                 *
* If root==-1, the data will be synced across all procs     *
************************************************************/
void SiloIO::syncMultiMeshData( std::map<AMP::Mesh::MeshID,siloMultiMeshData> &data, int root ) const
{
    PROFILE_START("syncMultiMeshData");
    // Convert the data to vectors
    std::vector<AMP::Mesh::MeshID> ids;
    std::vector<siloMultiMeshData> meshdata;
    std::map<AMP::Mesh::MeshID,siloMultiMeshData>::iterator iterator;
    int myRank = d_comm.getRank();
    for (iterator=data.begin(); iterator!=data.end(); ++iterator) {
        // Only send the base meshes that I own
        siloMultiMeshData tmp = iterator->second;
        tmp.meshes.resize(0);
        for (size_t i=0; i<iterator->second.meshes.size(); ++i) {
            if ( iterator->second.meshes[i].ownerRank==myRank )
                tmp.meshes.push_back( iterator->second.meshes[i] );
        }
        // Only the owner rank will send the variable list
        if ( tmp.ownerRank != myRank )
            tmp.varName.resize(0);
        // Only send the multimesh if there are base meshes that need to be sent or I own the mesh
        if ( !tmp.meshes.empty() || tmp.ownerRank==myRank ) {
            ids.push_back( iterator->first );
            meshdata.push_back( iterator->second );
        }
    }
    // Create buffers to store the data
    size_t send_size = 0;
    for (size_t i=0; i<meshdata.size(); ++i) {
        AMP_ASSERT(ids[i]==meshdata[i].id);
        send_size += meshdata[i].size();
    }
    char *send_buf = new char[send_size];
    char *ptr = send_buf;
    for (size_t i=0; i<meshdata.size(); ++i) {
        meshdata[i].pack(ptr);
        ptr = &ptr[meshdata[i].size()];
    }
    // Send the data and unpack the buffer to a vector
    size_t tot_num = d_comm.sumReduce(meshdata.size());
    if ( root==-1 ) {
        // Everybody gets a copy
        size_t tot_size = d_comm.sumReduce(send_size);
        char *recv_buf = new char[tot_size];
        meshdata.resize(tot_num);
        d_comm.allGather( send_buf, send_size, recv_buf );
        ptr = recv_buf;
        for (size_t i=0; i<tot_num; ++i) {
            meshdata[i] = siloMultiMeshData::unpack(ptr);
            ptr = &ptr[meshdata[i].size()];
        }
        delete [] recv_buf;
    } else {
        AMP_ASSERT(root>=0&&root<d_comm.getSize());
        // Only the root gets a copy
        // Note: the root already has his own data
        size_t max_size = d_comm.maxReduce(send_size);
        std::vector<int> recv_num(d_comm.getSize());
        d_comm.allGather((int)meshdata.size(),&recv_num[0]);
        if ( root == d_comm.getRank() ) {
            // Recieve all data
            meshdata.resize(0);
            meshdata.reserve(tot_num);
            char *recv_buf = new char[max_size];
            for (int i=0; i<d_comm.getSize(); ++i) {
                if ( i==root )
                    continue;
                int recv_size = d_comm.probe( i, 24987);
                AMP_ASSERT(recv_size<=(int)max_size);
                d_comm.recv( recv_buf, recv_size, i, false, 24987 );
                char *ptr = recv_buf;
                for (int j=0; j<recv_num[i]; ++j) {
                    siloMultiMeshData tmp = siloMultiMeshData::unpack(ptr);
                    ptr = &ptr[tmp.size()];
                    meshdata.push_back( tmp );
                }
            }
            delete [] recv_buf;
        } else {
            // Send my data
            d_comm.send( send_buf, send_size, root, 24987 );
        }
    }
    delete [] send_buf;
    // Add the meshes from other processors (keeping the existing meshes)
    for (size_t i=0; i<meshdata.size(); ++i) {
        iterator = data.find( meshdata[i].id );
        if ( iterator==data.end() ) {
            // Add the multimesh
            data.insert( std::pair<AMP::Mesh::MeshID,siloMultiMeshData>(meshdata[i].id,meshdata[i]) );
        } else {
            // Add the submeshes
            for (size_t j=0; j<meshdata[i].meshes.size(); ++j) {
                bool found = false;
                for (size_t k=0; k<iterator->second.meshes.size(); ++k) {
                    if ( meshdata[i].meshes[j].id==iterator->second.meshes[k].id &&
                         meshdata[i].meshes[j].meshName==iterator->second.meshes[k].meshName &&
                         meshdata[i].meshes[j].path==iterator->second.meshes[k].path &&
                         meshdata[i].meshes[j].path==iterator->second.meshes[k].file )
                        found = true;
                }
                if ( !found )
                    iterator->second.meshes.push_back( meshdata[i].meshes[j] );
            }
            // Add the variables if we don't have them yet
            if ( meshdata[i].varName.size() > 0 ) {
                if ( !iterator->second.varName.empty() )
                    AMP_ASSERT(iterator->second.varName.size()==meshdata[i].varName.size());
                iterator->second.varName = meshdata[i].varName;
            }
        }
    }
    PROFILE_STOP("syncMultiMeshData");
}


/************************************************************
* Function to syncronize a variable list                    *
* If root==-1, the data will be synced across all procs     *
************************************************************/
void SiloIO::syncVariableList( std::set<std::string> &data_set, int root ) const
{
    PROFILE_START("syncVariableList");
    std::vector<std::string> data(data_set.begin(),data_set.end());
    size_t N_local = data.size();
    size_t N_global = d_comm.sumReduce(N_local);
    size_t *size_local = new size_t[N_local];
    for (size_t i=0; i<N_local; ++i)
        size_local[i] = data[i].size();
    size_t *size_global = new size_t[N_global];
    d_comm.allGather( size_local, N_local, size_global );
    size_t tot_size_local = 0;
    for (size_t i=0; i<N_local; ++i)
        tot_size_local += size_local[i];
    size_t tot_size_global = 0;
    for (size_t i=0; i<N_global; ++i)
        tot_size_global += size_global[i];
    char *send_buf = new char[tot_size_local];
    char *recv_buf = new char[tot_size_global];
    size_t k=0;
    for (size_t i=0; i<N_local; ++i) {
        data[i].copy( &send_buf[k], data[i].size(), 0 );
        k += size_local[i];
    }
    if ( root==-1 ) {
        // Everybody gets a copy
        d_comm.allGather( send_buf, tot_size_local, recv_buf );
        k = 0;
        for (size_t i=0; i<N_global; ++i) {
            std::string tmp( &recv_buf[k], size_global[i] );
            data_set.insert(tmp);
            k += size_global[i];
        }
    } else {
        // Only the root gets a copy
        // Note: the root already has his own data
        AMP_ASSERT(root>=0&&root<d_comm.getSize());
        std::vector<int> recv_num(d_comm.getSize());
        d_comm.allGather((int)N_local,&recv_num[0]);
        if ( root == d_comm.getRank() ) {
            // Recieve all data
            int index = 0;
            for (int i=0; i<d_comm.getSize(); ++i) {
                if ( i==root ) {
                    index += recv_num[i];
                    continue;
                }
                int recv_size = d_comm.probe( i, 24987);
                d_comm.recv( recv_buf, recv_size, i, false, 24987 );
                k = 0;
                for (int j=0; j<recv_num[i]; ++j) {
                    std::string tmp( &recv_buf[k], size_global[index] );
                    data_set.insert(tmp);
                    k += size_global[index];
                    index++;
                }
                AMP_ASSERT((int)k==recv_size);
            }
        } else {
            // Send my data
            d_comm.send( send_buf, tot_size_local, root, 24987 );
        }
    }
    delete [] send_buf;
    delete [] recv_buf;
    delete [] size_local;
    delete [] size_global;
    PROFILE_STOP("syncVariableList");
}


/************************************************************
* Function to write the summary data                        *
************************************************************/
void SiloIO::writeSummary( std::string filename )
{
    PROFILE_START("writeSummary");
    // Add the siloBaseMeshData to the multimeshes
    std::map<AMP::Mesh::MeshID,siloMultiMeshData> multiMeshes = d_multiMeshes;
    std::map<AMP::Mesh::MeshID,siloMultiMeshData>::iterator iterator;
    for (iterator=multiMeshes.begin(); iterator!=multiMeshes.end(); ++iterator) {
        //AMP::Mesh::MeshID id = iterator->first;
        AMP::Mesh::Mesh::shared_ptr mesh = iterator->second.mesh;
        std::vector<AMP::Mesh::MeshID> base_ids = getMeshIDs( mesh );
        for (size_t i=0; i<base_ids.size(); ++i) {
            std::map<AMP::Mesh::MeshID,siloBaseMeshData>::iterator it = d_baseMeshes.find(base_ids[i]);
            if ( it != d_baseMeshes.end() ) {
                siloBaseMeshData data = it->second;
                AMP_ASSERT(it->first==data.id);
                iterator->second.meshes.push_back(data);
            }
        }
    }
    // Add the whole mesh
    /*if ( multiMeshes.size()==0 ) {
        siloMultiMeshData wholemesh;
        wholemesh.id = AMP::Mesh::MeshID((unsigned int)-1,0);
        wholemesh.name = "whole_mesh";
        std::map<AMP::Mesh::MeshID,siloBaseMeshData>::iterator iterator2;
        for (iterator2=d_baseMeshes.begin(); iterator2!=d_baseMeshes.end(); ++iterator2) {
            siloBaseMeshData data = iterator2->second;
            AMP_ASSERT(iterator2->first==data.id);
            wholemesh.meshes.push_back(data);
        }
        wholemesh.owner_rank = 0;
        multimeshes.insert( std::pair<AMP::Mesh::MeshID,siloMultiMeshData>(wholemesh.id,wholemesh) );
    }*/
    // Gather the results
    // Note: we only need to guarantee that rank 0 has all the data
    syncMultiMeshData( multiMeshes, 0 );
    syncVariableList( d_varNames, 0 );
    // Write the multimeshes
    if ( d_comm.getRank()==0 ) {
        DBfile  *FileHandle;
        FileHandle = DBOpen ( filename.c_str(), DB_HDF5, DB_APPEND );
        std::map<AMP::Mesh::MeshID,siloMultiMeshData>::iterator it;
        // Create the subdirectories
        PROFILE_START("create directories");
        std::set<std::string> subdirs;
        for (it=multiMeshes.begin(); it!=multiMeshes.end(); ++it) {
            siloMultiMeshData data = it->second;
            size_t pos = data.name.find_last_of("/");
            if ( pos!=std::string::npos )
                subdirs.insert( data.name.substr(0,pos) );
        }
        for (std::set<std::string>::iterator it2=subdirs.begin(); it2!=subdirs.end(); ++it2)
            createSiloDirectory( FileHandle, *it2 );
        PROFILE_STOP("create directories");
        // Create the multimeshes
        PROFILE_START("write multimeshes");
        for (it=multiMeshes.begin(); it!=multiMeshes.end(); ++it) {
            siloMultiMeshData data = it->second;
            std::vector<std::string> meshNames(data.meshes.size());
            for (size_t i=0; i<data.meshes.size(); ++i)
                meshNames[i] = data.meshes[i].file+":"+data.meshes[i].path+"/"+data.meshes[i].meshName;
            char **meshnames = new char*[data.meshes.size()];
            int *meshtypes = new int[data.meshes.size()];
            for (size_t i=0; i<data.meshes.size(); ++i) {
                meshnames[i] = (char*) meshNames[i].c_str();
                meshtypes[i] = DB_UCDMESH;
            }
            std::string tree_name = data.name+"_tree";
            DBoptlist *optList = DBMakeOptlist(10);
            DBAddOption( optList, DBOPT_MRGTREE_NAME, (char*)tree_name.c_str() );
            DBPutMultimesh( FileHandle, data.name.c_str(), meshNames.size(), meshnames, meshtypes, NULL );
            DBFreeOptlist( optList );
            delete [] meshnames;
            delete [] meshtypes;
        }
        PROFILE_STOP("write multimeshes");
        // Generate the multi-variables
        PROFILE_START("write multivariables");
        for (it=multiMeshes.begin(); it!=multiMeshes.end(); ++it) {
            siloMultiMeshData data = it->second;
            //std::cout << data.name << std::endl;
            for (size_t i=0; i<data.varName.size(); ++i) {
                std::string varName = data.varName[i];
                std::vector<std::string> varNames(data.meshes.size());
                char **varnames = new char*[data.meshes.size()];
                int *vartypes = new int[data.meshes.size()];
                for (size_t i=0; i<data.meshes.size(); ++i) {
                    std::stringstream  stream;
                    stream << data.meshes[i].rank;
                    varNames[i] = data.meshes[i].file+":"+data.meshes[i].path+"/"+varName+"P"+stream.str();
                    varnames[i] = (char*) varNames[i].c_str();
                    vartypes[i] = DB_UCDVAR;
                    //std::cout << varNames[i] << std::endl;
                }
                int varSize = 0;
                for (size_t i=0; i<data.meshes[0].varName.size(); ++i) {
                    if ( data.meshes[0].varName[i]==varName ) {
                        varSize = data.meshes[0].varSize[i];
                        break;
                    }
                }
                std::string multiMeshName = data.name;
                std::string visitVarName = multiMeshName+"_"+varName;
                DBoptlist *opts = NULL;
                //DBoptlist *opts = DBMakeOptlist(1);
                //DBAddOption( opts, DBOPT_MMESH_NAME, (char*) multiMeshName.c_str() );
                if ( varSize==1 || varSize==d_dim || varSize==d_dim*d_dim ) {
                    // We are writing a scalar, vector, or tensor variable
                    DBPutMultivar( FileHandle, visitVarName.c_str(), varNames.size(), varnames, vartypes, opts );
                } else {
                    // Write each component
                    for (int j=0; j<varSize; ++j) {
                        std::stringstream  stream;
                        stream << "_" << j;
                        std::string postfix = stream.str();
                        std::vector<std::string> varNames2(data.meshes.size());
                        for (size_t i=0; i<data.meshes.size(); ++i) {
                            varNames2[i] = varNames[i] + postfix;
                            varnames[i] = (char*) varNames2[i].c_str();
                        }
                        DBPutMultivar( FileHandle, (visitVarName+postfix).c_str(), varNames.size(), varnames, vartypes, opts );
                    }
                }
                //DBFreeOptlist( opts );
                delete [] varnames;
                delete [] vartypes;
            }
        }
        PROFILE_STOP("write multivariables");
        DBClose ( FileHandle );
    }
    PROFILE_STOP("writeSummary");
}


/************************************************************
* Functions for siloBaseMeshData                            *
************************************************************/
size_t SiloIO::siloBaseMeshData::size()
{
    size_t N_bytes = sizeof(AMP::Mesh::MeshID);     // Store the mesh id
    N_bytes += sizeof(int);                         // Store the processor rank
    N_bytes += sizeof(int);                         // Store the owner rank
    N_bytes += sizeof(int)+meshName.size();         // Store the mesh name
    N_bytes += sizeof(int)+path.size();             // Store the mesh path
    N_bytes += sizeof(int)+file.size();             // Store the mesh file
    N_bytes += sizeof(int);                         // Store the number of variables
    for (size_t i=0; i<varName.size(); ++i) {
        N_bytes += sizeof(int)+varName[i].size();   // Store the variable name
        N_bytes += sizeof(AMP::Mesh::GeomType);     // Store the variable type
        N_bytes += sizeof(int);                     // Store the number of unknowns per point
    }
    return N_bytes;
}
void SiloIO::siloBaseMeshData::pack( char* ptr )
{
    size_t pos = 0;
    // Store the mesh id
    *((AMP::Mesh::MeshID*) &ptr[pos]) = id;
    pos += sizeof(AMP::Mesh::MeshID);
    // Store the mesh rank
    *((int*) &ptr[pos]) = rank;
    pos += sizeof(int);
    // Store the owner rank
    *((int*) &ptr[pos]) = ownerRank;
    pos += sizeof(int);
    // Store the mesh name
    *((int*) &ptr[pos]) = (int) meshName.size();
    pos += sizeof(int);
    meshName.copy( (char*) &ptr[pos], meshName.size(), 0 );
    pos += meshName.size();
    // Store the mesh path
    *((int*) &ptr[pos]) = (int) path.size();
    pos += sizeof(int);
    path.copy( (char*) &ptr[pos], path.size(), 0 );
    pos += path.size();
    // Store the mesh file
    *((int*) &ptr[pos]) = (int) file.size();
    pos += sizeof(int);
    file.copy( (char*) &ptr[pos], file.size(), 0 );
    pos += file.size();
    // Store the number of variables
    *((int*) &ptr[pos]) = (int) varName.size();
    pos += sizeof(int);
    for (size_t i=0; i<varName.size(); ++i) {
        // Store the variable name
        *((int*) &ptr[pos]) = (int) varName[i].size();
        pos += sizeof(int);
        varName[i].copy( (char*) &ptr[pos], varName[i].size(), 0 );
        pos += varName[i].size();
        // Store the variable type
        *((AMP::Mesh::GeomType*) &ptr[pos]) = varType[i];
        pos += sizeof(AMP::Mesh::GeomType);
        // Store the number of unknowns per point
        *((int*) &ptr[pos]) = (int) varSize[i];
        pos += sizeof(int);
    }
    AMP_ASSERT(pos==size());
}
SiloIO::siloBaseMeshData SiloIO::siloBaseMeshData::unpack( char* ptr )
{
    siloBaseMeshData data;
    size_t pos = 0;
    // Store the mesh id
    data.id = *((AMP::Mesh::MeshID*) &ptr[pos]);
    pos += sizeof(AMP::Mesh::MeshID);
    // Store the mesh rank
    data.rank = *((int*) &ptr[pos]);
    pos += sizeof(int);
    // Store the owner rank
    data.ownerRank = *((int*) &ptr[pos]);
    pos += sizeof(int);
    // Store the mesh name
    int size = *((int*) &ptr[pos]);
    pos += sizeof(int);
    data.meshName = std::string( (char*) &ptr[pos], size );
    pos += size;
    // Store the mesh path
    size = *((int*) &ptr[pos]);
    pos += sizeof(int);
    data.path = std::string( (char*) &ptr[pos], size );
    pos += size;
    // Store the mesh file
    size = *((int*) &ptr[pos]);
    pos += sizeof(int);
    data.file = std::string( (char*) &ptr[pos], size );
    pos += size;
    // Store the variables
    int N_var = *((int*) &ptr[pos]);
    pos += sizeof(int);
    data.varName.resize(N_var);
    data.varType.resize(N_var);
    data.varSize.resize(N_var);
    #ifdef USE_AMP_VECTORS
        data.vec.resize(N_var);  // Set the vec to NULL
    #endif
    for (int i=0; i<N_var; ++i) {
        // Store the variable name
        size = *((int*) &ptr[pos]);
        pos += sizeof(int);
        data.varName[i] = std::string( (char*) &ptr[pos], size );
        pos += size;
        // Store the variable type
        data.varType[i] = *((AMP::Mesh::GeomType*) &ptr[pos]);
        pos += sizeof(AMP::Mesh::GeomType);
        // Store the number of unknowns per point
        data.varSize[i] = *((int*) &ptr[pos]);
        pos += sizeof(int);
    }
    AMP_ASSERT(pos==data.size());
    return data;
}


/************************************************************
* Functions for siloMultiMeshData                           *
************************************************************/
size_t SiloIO::siloMultiMeshData::size()
{
    size_t N_bytes = sizeof(AMP::Mesh::MeshID);     // Store the mesh id
    N_bytes += sizeof(int);                         // Store the owner rank
    N_bytes += sizeof(int)+name.size();             // Store the mesh name
    N_bytes += sizeof(int);                         // Store the number of sub meshes
    for (size_t i=0; i<meshes.size(); ++i)
        N_bytes += meshes[i].size();                // Store the sub meshes
    N_bytes += sizeof(int);                         // Store the number of variables
    for (size_t i=0; i<varName.size(); ++i) {
        N_bytes += sizeof(int);                     // Store the length of the variable name
        N_bytes += varName[i].size();               // Store the variable name
    }
    return N_bytes;
}
void SiloIO::siloMultiMeshData::pack( char* ptr )
{
    size_t pos = 0;
    // Store the mesh id
    *((AMP::Mesh::MeshID*) &ptr[pos]) = id;
    pos += sizeof(AMP::Mesh::MeshID);
    // Store the owner rank
    *((int*) &ptr[pos]) = ownerRank;
    pos += sizeof(int);
    // Store name
    *((int*) &ptr[pos]) = (int) name.size();
    pos += sizeof(int);
    name.copy( (char*) &ptr[pos], name.size(), 0 );
    pos += name.size();
    // Store the base meshes
    *((int*) &ptr[pos]) = (int) meshes.size();
    pos += sizeof(int);
    for (size_t i=0; i<meshes.size(); ++i) {
        meshes[i].pack( &ptr[pos] );
        pos += meshes[i].size();
    }
    // Store the variables
    *((int*) &ptr[pos]) = (int) varName.size();
    pos += sizeof(int);
    for (size_t i=0; i<varName.size(); ++i) {
        *((int*) &ptr[pos]) = (int) varName[i].size();
        pos += sizeof(int);
        varName[i].copy( (char*) &ptr[pos], varName[i].size(), 0 );
        pos += varName[i].size();
    }
    AMP_ASSERT(pos==size());
}
SiloIO::siloMultiMeshData SiloIO::siloMultiMeshData::unpack( char* ptr )
{
    siloMultiMeshData data;
    size_t pos = 0;
    // Store the mesh id
    data.id = *((AMP::Mesh::MeshID*) &ptr[pos]);
    pos += sizeof(AMP::Mesh::MeshID);
    // Store the owner rank
    data.ownerRank = *((int*) &ptr[pos]);
    pos += sizeof(int);
    // Store name
    int size = *((int*) &ptr[pos]);
    pos += sizeof(int);
    data.name = std::string( (char*) &ptr[pos], size );
    pos += size;
    // Store the base meshes
    int N_meshes = *((int*) &ptr[pos]);
    pos += sizeof(int);
    data.meshes.resize(N_meshes);
    for (int i=0; i<N_meshes; ++i) {
        data.meshes[i] = siloBaseMeshData::unpack( &ptr[pos] );
        pos += data.meshes[i].size();
    }
    // Store the variables
    int N_var = *((int*) &ptr[pos]);
    pos += sizeof(int);
    data.varName = std::vector<std::string>(N_var);
    for (size_t i=0; i<data.varName.size(); ++i) {
        size = *((int*) &ptr[pos]);
        pos += sizeof(int);
        data.varName[i] = std::string( (char*) &ptr[pos], size );
        pos += size;
    }
    AMP_ASSERT(pos==data.size());
    return data;
}


/************************************************************
* Some utilit functions                                     *
************************************************************/
void createSiloDirectory( DBfile *FileHandle, std::string path )
{
    // Create a subdirectory tree from the current working path if it does not exist
    char current_dir[256];
    DBGetDir( FileHandle, current_dir );
    // Get the list of directories that may need to be created
    std::vector<std::string> subdirs;
    std::string path2 = path + "/";
    while ( !path2.empty() ) {
        size_t pos = path2.find("/");
        if ( pos > 0 ) {
            subdirs.push_back( path2.substr(0,pos) );
        }
        path2.erase(0,pos+1);
    }
    // Create the directories as necessary
    for (size_t i=0; i<subdirs.size(); ++i) {
        DBtoc *toc = DBGetToc( FileHandle );
        bool exists = false;
        for (int j=0; j<toc->ndir; ++j) {
            if ( subdirs[i].compare(toc->dir_names[j])==0 )
                exists = true;
        }
        if ( !exists )
            DBMkDir ( FileHandle, subdirs[i].c_str() );
        DBSetDir( FileHandle, subdirs[i].c_str() );
    }
    // Return back to the original working directory
    DBSetDir( FileHandle, current_dir );
}



#else
void SiloIO::readFile( const std::string& ) { AMP_ERROR("SILO not configured"); }
void SiloIO::writeFile( const std::string&, size_t ) { AMP_ERROR("SILO not configured"); }
void SiloIO::registerMesh( AMP::Mesh::Mesh::shared_ptr, int, std::string ) { AMP_ERROR("SILO not configured"); }
#ifdef USE_AMP_VECTORS
void SiloIO::registerVector( AMP::LinearAlgebra::Vector::shared_ptr vec, AMP::Mesh::Mesh::shared_ptr,
    AMP::Mesh::GeomType, const std::string& ) { AMP_ERROR("SILO not configured"); }
#endif

#endif


} // Mesh namespace
} // AMP namespace


