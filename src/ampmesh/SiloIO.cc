#include "SiloIO.h"
#include "ampmesh/MultiMesh.h"

namespace AMP { 
namespace Mesh {


/************************************************************
* Constructor                                               *
************************************************************/
SiloIO::SiloIO( )
{
    d_comm = AMP_MPI(AMP_COMM_WORLD);
    dim = -1;
}


#ifdef USE_SILO

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
    std::string fname = fname_in + "." + getExtension();
    for (int i=0; i<d_comm.getSize(); i++) {
        if ( d_comm.getRank()==i ) {
            // Open the file
            DBfile *FileHandle;
            if ( d_comm.getRank()==0 ) {
                FileHandle = DBCreate( fname.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_HDF5 );
                // Create the necessary directories
                DBMkDir( FileHandle, "All_Meshes" );
            } else {
                FileHandle = DBOpen ( fname.c_str(), DB_HDF5, DB_APPEND );
            }
            // Write the base meshes
            std::map<AMP::Mesh::MeshID,siloBaseMeshData>::iterator iterator;
            for (iterator=d_baseMeshes.begin(); iterator!=d_baseMeshes.end(); iterator++) {
                siloBaseMeshData &data = iterator->second;
                data.file = fname;
                AMP_ASSERT(data.id==iterator->first);
                writeMesh( FileHandle, iterator->second );
            }
            // Close the file
            DBClose ( FileHandle );
        }
        d_comm.barrier();
    }
    // Write the summary results (multimeshes, multivariables, etc.)
    writeSummary( fname );
}


/************************************************************
* Function to register a mesh with silo                     *
************************************************************/
void SiloIO::registerMesh( AMP::Mesh::Mesh::shared_ptr mesh )
{ 
    if ( dim == -1 )
        dim = mesh->getDim();
    else
        AMP_INSIST(dim==mesh->getDim(),"All meshes must have the same number of physical dimensions");
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
        data.meshName = "rank_" + rank;
        data.path = "All_Meshes/"+mesh->getName();
        d_baseMeshes.insert( std::pair<AMP::Mesh::MeshID,siloBaseMeshData>(mesh->meshID(),data) );
        siloMultiMeshData data2;
        data2.id = mesh->meshID();
        data2.mesh = mesh;
        data2.name = mesh->getName();
        d_multiMeshes.insert( std::pair<AMP::Mesh::MeshID,siloMultiMeshData>(mesh->meshID(),data2) );
    } else {
        // We are dealining with a multimesh, register the current mesh and all base meshes
        std::vector<AMP::Mesh::MeshID> ids = mesh->getBaseMeshIDs();
        for (size_t i=0; i<ids.size(); i++) {
            AMP::Mesh::Mesh::shared_ptr mesh2 = mesh->Subset(ids[i]);
            if ( mesh2.get() != NULL )
                registerMesh( mesh2 );
        }
        siloMultiMeshData data;
        data.id = mesh->meshID();
        data.mesh = mesh;
        data.name = mesh->getName();
        d_multiMeshes.insert( std::pair<AMP::Mesh::MeshID,siloMultiMeshData>(mesh->meshID(),data) );
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
        for (size_t i=0; i<meshes.size(); i++) {
            std::vector<AMP::Mesh::MeshID> ids2 = getMeshIDs( meshes[i] );
            for (size_t j=0; j<ids2.size(); j++)
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
    for (size_t i=0; i<iterator1.size(); i++) {
        DOFs->getDOFs( iterator1->globalID(), dofs );
        AMP_ASSERT((int)dofs.size()==DOFsPerPoint);
        ++iterator1;
    }
    // Register the vector with the appropriate base meshes
    std::string name = vec->type();
    if ( !name.empty() )
        name = name_in;
    std::vector<AMP::Mesh::MeshID> ids = getMeshIDs( mesh );
    for (size_t i=0; i<ids.size(); i++) {
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
    // Add the variable name to the list of meshes
    d_varNames.insert(name);
}
#endif


/************************************************************
* Function to write a mesh                                  *
************************************************************/
void SiloIO::writeMesh( DBfile *FileHandle, const siloBaseMeshData &data )
{ 
    AMP::Mesh::Mesh::shared_ptr mesh = data.mesh;
    int dim = mesh->getDim();
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
    for (size_t i=0; i<node_iterator.size(); i++) {
        nodelist_ids[i] = node_iterator->globalID();
        ++node_iterator;
    }
    AMP::Utilities::quicksort(nodelist_ids);
    double *coord[3];
    for (int i=0; i<dim; i++)
        coord[i] = new double[node_iterator.size()];
    node_iterator = mesh->getIterator(AMP::Mesh::Vertex,1);
    for (size_t i=0; i<node_iterator.size(); i++) {
        size_t index = AMP::Utilities::findfirst( nodelist_ids, node_iterator->globalID() );
        AMP_ASSERT(nodelist_ids[index]==node_iterator->globalID());
        std::vector<double> elem_coord = node_iterator->coord();
        for (int j=0; j<dim; j++)
            coord[j][index] = elem_coord[j];
        ++node_iterator;
    }
    elem_iterator = mesh->getIterator(mesh->getGeomType(),0);
    std::vector<int> nodelist;
    nodelist.reserve(shapesize*elem_iterator.size());
    for (size_t i=0; i<elem_iterator.size(); i++) {
        nodes = elem_iterator->getElements(AMP::Mesh::Vertex);
        AMP_INSIST((int)nodes.size()==shapesize,"Mixed element types is currently not supported");
        for (size_t j=0; j<nodes.size(); j++) {
            size_t index = AMP::Utilities::findfirst( nodelist_ids, nodes[j].globalID() );
            AMP_ASSERT(nodelist_ids[index]==nodes[j].globalID());
            nodelist.push_back( (int) index );
        }
        ++elem_iterator;
    }
    // Write the elements (connectivity)
    if ( mesh->getComm().getRank()==0 ) {
        DBSetDir( FileHandle, "All_Meshes" );
        DBMkDir( FileHandle, mesh->getName().c_str() );
        DBSetDir( FileHandle, "/" );
    }
    DBSetDir( FileHandle, data.path.c_str() );
    std::stringstream  stream;
    stream << data.rank;
    std::string rank = stream.str();
    std::string meshName = data.meshName;
    std::string zoneName = "zone_" + rank;
    AMP::Mesh::MeshIterator  element_iterator = mesh->getIterator(mesh->getGeomType(),0);
    int num_elems = (int) element_iterator.size();
    DBPutZonelist2( FileHandle, zoneName.c_str(), num_elems, dim, 
        &nodelist[0], nodelist.size(), 0, 0, 0, &shapetype, &shapesize, &shapecnt, 1, 0 );
    // Write the mesh
    DBPutUcdmesh( FileHandle, meshName.c_str(), mesh->getDim(),
        NULL, coord, node_iterator.size(), nodelist.size(),
        zoneName.c_str(), 0, DB_DOUBLE, 0 );
    for (int i=0; i<dim; i++)
        delete [] coord[i];
    // Write the variables
    #ifdef USE_AMP_VECTORS
        for (size_t i=0; i<data.varName.size(); i++) {
            AMP::Discretization::DOFManager::shared_ptr DOFs = data.vec[i]->getDOFManager();
            int nvar = 0;
            int centering;
            double **var = new double*[data.varSize[i]];
            const char *varnames[] = {"1","2","3"};
            if ( data.varType[i]==AMP::Mesh::Vertex ) {
                // We are saving node-centered data
                centering = DB_NODECENT;
                nvar = (int) nodelist_ids.size();
                for (int j=0; j<data.varSize[i]; j++)
                    var[j] = new double[nvar];
                std::vector<size_t> dofs(data.varSize[i]);
                std::vector<double> vals(data.varSize[i]);
                for (int j=0; j<nvar; j++) {
                    DOFs->getDOFs( nodelist_ids[j], dofs );
                    AMP_ASSERT((int)dofs.size()==data.varSize[i]);
                    data.vec[i]->getValuesByGlobalID ( data.varSize[i], &dofs[0], &vals[0] );
                    for (int k=0; k<data.varSize[i]; k++)
                        var[k][j] = vals[k];
                }
            } else if ( data.varType[i]==mesh->getGeomType() ) {
                // We are saving cell-centered data
                centering = DB_ZONECENT;
                nvar = (int) num_elems;
                for (int j=0; j<data.varSize[i]; j++)
                    var[j] = new double[nvar];
                std::vector<size_t> dofs(data.varSize[i]);
                std::vector<double> vals(data.varSize[i]);
                AMP::Mesh::MeshIterator it = element_iterator.begin();
                for (int j=0; j<nvar; j++) {
                    DOFs->getDOFs( it->globalID(), dofs );
                    data.vec[i]->getValuesByGlobalID ( data.varSize[i], &dofs[0], &vals[0] );
                    for (int k=0; k<data.varSize[i]; k++)
                        var[k][j] = vals[k];
                    ++it;
                }
            } else {
                // We are storing edge or face data
                AMP_ERROR("The silo writer currently only supports Vertex and Cell data");
            }
            std::string varNameRank = data.varName[i]+"P"+rank;
            if ( data.varSize[i]==1 || data.varSize[i]==dim || data.varSize[i]==dim*dim ) {
                // We are writing a scalar, vector, or tensor variable
                DBPutUcdvar( FileHandle, varNameRank.c_str(), meshName.c_str(), 
                    data.varSize[i], (char**) varnames, var, nvar, NULL, 0, DB_DOUBLE, centering, NULL);
            } else {
                // Write each component
                for (int j=0; j<data.varSize[i]; j++) {
                    std::stringstream  stream;
                    stream << varNameRank << "_" << j;
                    DBPutUcdvar( FileHandle, stream.str().c_str(), meshName.c_str(), 
                        1, (char**) varnames, &var[j], nvar, NULL, 0, DB_DOUBLE, centering, NULL);
                }
            }
            for (int j=0; j<data.varSize[i]; j++)
                delete [] var[j];
            delete [] var;
        }
    #endif
    // Change the directory back to root
    DBSetDir( FileHandle, "/" );
}


/************************************************************
* Function to syncronize the multimesh data                 *
************************************************************/
void SiloIO::syncMultiMeshData( std::map<AMP::Mesh::MeshID,siloMultiMeshData> &data ) const
{
    // Convert the data to vectors
    std::vector<AMP::Mesh::MeshID> ids;
    std::vector<siloMultiMeshData> meshdata;
    std::map<AMP::Mesh::MeshID,siloMultiMeshData>::iterator iterator;
    for (iterator=data.begin(); iterator!=data.end(); iterator++) {
        ids.push_back( iterator->first );
        meshdata.push_back( iterator->second );
    }
    // Create buffers to store the data
    size_t send_size = 0;
    for (size_t i=0; i<meshdata.size(); i++) {
        AMP_ASSERT(ids[i]==meshdata[i].id);
        send_size += meshdata[i].size();
    }
    char *send_buf = new char[send_size];
    char *ptr = send_buf;
    for (size_t i=0; i<meshdata.size(); i++) {
        meshdata[i].pack(ptr);
        ptr = &ptr[meshdata[i].size()];
    }
    // Send the data
    size_t tot_num = d_comm.sumReduce(meshdata.size());
    size_t tot_size = d_comm.sumReduce(send_size);
    char *recv_buf = new char[tot_size];
    d_comm.allGather( send_buf, send_size, recv_buf );
    // Unpack the buffer to a vector
    meshdata = std::vector<siloMultiMeshData>(tot_num);
    ptr = recv_buf;
    for (size_t i=0; i<tot_num; i++) {
        meshdata[i] = siloMultiMeshData::unpack(ptr);
        ptr = &ptr[meshdata[i].size()];
    }
    delete [] send_buf;
    delete [] recv_buf;
    // Combine the results
    data = std::map<AMP::Mesh::MeshID,siloMultiMeshData>();
    while ( meshdata.size() > 0 ) {
        siloMultiMeshData current = meshdata[0];
        std::vector<siloMultiMeshData>::iterator it = meshdata.begin();
        it = meshdata.erase(it);
        while ( it != meshdata.end() ) {
            siloMultiMeshData tmp = *it;
            if ( tmp.id==current.id ) {
                for (size_t i=0; i<tmp.meshes.size(); i++)
                    current.meshes.push_back(tmp.meshes[i]);
                it = meshdata.erase(it);
            } else {
                it++;
            }
        }
        data.insert( std::pair<AMP::Mesh::MeshID,siloMultiMeshData>(current.id,current) );
    }
}


/************************************************************
* Function to syncronize a variable list                    *
************************************************************/
void SiloIO::syncVariableList( std::set<std::string> &data_set ) const
{
    std::vector<std::string> data(data_set.begin(),data_set.end());
    size_t N_local = data.size();
    size_t N_global = d_comm.sumReduce(N_local);
    size_t *size_local = new size_t[N_local];
    for (size_t i=0; i<N_local; i++)
        size_local[i] = data[i].size();
    size_t *size_global = new size_t[N_global];
    d_comm.allGather( size_local, N_local, size_global );
    size_t tot_size_local = 0;
    for (size_t i=0; i<N_local; i++)
        tot_size_local += size_local[i];
    size_t tot_size_global = 0;
    for (size_t i=0; i<N_global; i++)
        tot_size_global += size_global[i];
    char *send_buf = new char[tot_size_local];
    char *recv_buf = new char[tot_size_global];
    size_t k=0;
    for (size_t i=0; i<N_local; i++) {
        data[i].copy( &send_buf[k], data[i].size(), 0 );
        k += size_local[i];
    }
    d_comm.allGather( send_buf, tot_size_local, recv_buf );
    k = 0;
    for (size_t i=0; i<N_global; i++) {
        std::string tmp( &recv_buf[k], size_global[i] );
        data_set.insert(tmp);
        k += size_global[i];
    }
    delete [] send_buf;
    delete [] recv_buf;
    delete [] size_local;
    delete [] size_global;
}


/************************************************************
* Function to write the summary data                        *
************************************************************/
void SiloIO::writeSummary( std::string filename )
{
    int dim2 = d_comm.bcast(dim,0);
    if ( dim == -1 )
        dim = dim2;
    else
        AMP_ASSERT(dim2==dim);
    // Add the siloBaseMeshData to the multimeshes
    std::map<AMP::Mesh::MeshID,siloMultiMeshData> multiMeshes = d_multiMeshes;
    std::map<AMP::Mesh::MeshID,siloMultiMeshData>::iterator iterator;
    for (iterator=multiMeshes.begin(); iterator!=multiMeshes.end(); iterator++) {
        //AMP::Mesh::MeshID id = iterator->first;
        AMP::Mesh::Mesh::shared_ptr mesh = iterator->second.mesh;
        std::vector<AMP::Mesh::MeshID> base_ids = getMeshIDs( mesh );
        for (size_t i=0; i<base_ids.size(); i++) {
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
        for (iterator2=d_baseMeshes.begin(); iterator2!=d_baseMeshes.end(); iterator2++) {
            siloBaseMeshData data = iterator2->second;
            AMP_ASSERT(iterator2->first==data.id);
            wholemesh.meshes.push_back(data);
        }
        multimeshes.insert( std::pair<AMP::Mesh::MeshID,siloMultiMeshData>(wholemesh.id,wholemesh) );
    }*/
    // Gather the results
    syncMultiMeshData( multiMeshes );
    syncVariableList( d_varNames );
    // Write the multimeshes
    if ( d_comm.getRank()==0 ) {
        DBfile  *FileHandle;
        FileHandle = DBOpen ( filename.c_str(), DB_HDF5, DB_APPEND );
        DBMkDir ( FileHandle, "test" );
        //DBSetDir( FileHandle, "test" );
        std::map<AMP::Mesh::MeshID,siloMultiMeshData>::iterator it;
        for (it=multiMeshes.begin(); it!=multiMeshes.end(); it++) {
            // Create the multimesh            
            siloMultiMeshData data = it->second;
            std::vector<std::string> meshNames(data.meshes.size());
            for (size_t i=0; i<data.meshes.size(); i++)
                meshNames[i] = data.meshes[i].file+":"+data.meshes[i].path+"/"+data.meshes[i].meshName;
            char **meshnames = new char*[data.meshes.size()];
            int *meshtypes = new int[data.meshes.size()];
            for (size_t i=0; i<data.meshes.size(); i++) {
                meshnames[i] = (char*) meshNames[i].c_str();
                meshtypes[i] = DB_UCDMESH;
            }
            //DBPutMultimesh( FileHandle, "all", meshNames.size(), meshnames, meshtypes, NULL );
            std::string tree_name = data.name+"_tree";
            DBoptlist *optList = DBMakeOptlist(10);
            DBAddOption( optList, DBOPT_MRGTREE_NAME, (char*)tree_name.c_str() );
            DBPutMultimesh( FileHandle, data.name.c_str(), meshNames.size(), meshnames, meshtypes, NULL );
            DBFreeOptlist( optList );
            delete [] meshnames;
            delete [] meshtypes;
        }
        //DBSetDir( FileHandle, "/" );
        // Generate the multi-variables
        for (it=multiMeshes.begin(); it!=multiMeshes.end(); it++) {
            siloMultiMeshData data = it->second;
            //std::cout << data.name << std::endl;
            for (size_t i=0; i<data.varName.size(); i++) {
                std::string varName = data.varName[i];
                std::vector<std::string> varNames(data.meshes.size());
                char **varnames = new char*[data.meshes.size()];
                int *vartypes = new int[data.meshes.size()];
                for (size_t i=0; i<data.meshes.size(); i++) {
                    std::stringstream  stream;
                    stream << data.meshes[i].rank;
                    varNames[i] = data.meshes[i].file+":"+data.meshes[i].path+"/"+varName+"P"+stream.str();
                    varnames[i] = (char*) varNames[i].c_str();
                    vartypes[i] = DB_UCDVAR;
                    //std::cout << varNames[i] << std::endl;
                }
                int varSize = 0;
                for (size_t i=0; i<data.meshes[0].varName.size(); i++) {
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
                if ( varSize==1 || varSize==dim || varSize==dim*dim ) {
                    // We are writing a scalar, vector, or tensor variable
                    DBPutMultivar( FileHandle, visitVarName.c_str(), varNames.size(), varnames, vartypes, opts );
                } else {
                    // Write each component
                    for (int j=0; j<varSize; j++) {
                        std::stringstream  stream;
                        stream << "_" << j;
                        std::string postfix = stream.str();
                        std::vector<std::string> varNames2(data.meshes.size());
                        for (size_t i=0; i<data.meshes.size(); i++) {
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
        DBClose ( FileHandle );
    }
}


/************************************************************
* Functions for siloBaseMeshData                            *
************************************************************/
size_t SiloIO::siloBaseMeshData::size()
{
    size_t N_bytes = sizeof(AMP::Mesh::MeshID);     // Store the mesh id
    N_bytes += sizeof(int);                         // Store the processor rank
    N_bytes += sizeof(int)+meshName.size();         // Store the mesh name
    N_bytes += sizeof(int)+path.size();             // Store the mesh path
    N_bytes += sizeof(int)+file.size();             // Store the mesh file
    N_bytes += sizeof(int);                         // Store the number of variables
    for (size_t i=0; i<varName.size(); i++) {
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
    for (size_t i=0; i<varName.size(); i++) {
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
    for (int i=0; i<N_var; i++) {
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
    N_bytes += sizeof(int)+name.size();             // Store the mesh name
    N_bytes += sizeof(int);                         // Store the number of sub meshes
    for (size_t i=0; i<meshes.size(); i++)
        N_bytes += meshes[i].size();                // Store the sub meshes
    N_bytes += sizeof(int);                         // Store the number of variables
    for (size_t i=0; i<varName.size(); i++) {
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
    // Store name
    *((int*) &ptr[pos]) = (int) name.size();
    pos += sizeof(int);
    name.copy( (char*) &ptr[pos], name.size(), 0 );
    pos += name.size();
    // Store the base meshes
    *((int*) &ptr[pos]) = (int) meshes.size();
    pos += sizeof(int);
    for (size_t i=0; i<meshes.size(); i++) {
        meshes[i].pack( &ptr[pos] );
        pos += meshes[i].size();
    }
    // Store the variables
    *((int*) &ptr[pos]) = (int) varName.size();
    pos += sizeof(int);
    for (size_t i=0; i<varName.size(); i++) {
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
    // Store name
    int size = *((int*) &ptr[pos]);
    pos += sizeof(int);
    data.name = std::string( (char*) &ptr[pos], size );
    pos += size;
    // Store the base meshes
    int N_meshes = *((int*) &ptr[pos]);
    pos += sizeof(int);
    data.meshes.resize(N_meshes);
    for (int i=0; i<N_meshes; i++) {
        data.meshes[i] = siloBaseMeshData::unpack( &ptr[pos] );
        pos += data.meshes[i].size();
    }
    // Store the variables
    int N_var = *((int*) &ptr[pos]);
    pos += sizeof(int);
    data.varName = std::vector<std::string>(N_var);
    for (size_t i=0; i<data.varName.size(); i++) {
        size = *((int*) &ptr[pos]);
        pos += sizeof(int);
        data.varName[i] = std::string( (char*) &ptr[pos], size );
        pos += size;
    }
    AMP_ASSERT(pos==data.size());
    return data;
}


#else
void SiloIO::readFile( const std::string& ) { AMP_ERROR("SILO not configured"); }
void SiloIO::writeFile( const std::string&, size_t ) { AMP_ERROR("SILO not configured"); }
void SiloIO::registerMesh( AMP::Mesh::Mesh::shared_ptr ) { AMP_ERROR("SILO not configured"); }
#ifdef USE_AMP_VECTORS
void SiloIO::registerVector( AMP::LinearAlgebra::Vector::shared_ptr vec, AMP::Mesh::Mesh::shared_ptr,
    AMP::Mesh::GeomType, const std::string& ) { AMP_ERROR("SILO not configured"); }
#endif

#endif


} // Mesh namespace
} // AMP namespace


