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
    // Create the master file
    DBfile  *FileHandle;
    std::string fname = fname_in + "." + getExtension();
    std::map<AMP::Mesh::MeshID,std::string> pathnames;
    for (int i=0; i<d_comm.getSize(); i++) {
        if ( d_comm.getRank()==i ) {
            // Open the file
            DBfile *FileHandle;
            if ( d_comm.getRank()==0 ) {
                FileHandle = DBCreate( fname.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_HDF5 );
                DBMkDir( FileHandle, "All_Meshes" );
            } else {
                FileHandle = DBOpen ( fname.c_str(), DB_HDF5, DB_APPEND );
            }
            // Write the base meshes
            std::map<AMP::Mesh::MeshID,AMP::Mesh::Mesh::shared_ptr>::iterator iterator;
            for (iterator=d_baseMeshes.begin(); iterator!=d_baseMeshes.end(); iterator++) {
                std::string dbpath = writeMesh( FileHandle, iterator->second );
                pathnames.insert( std::pair<AMP::Mesh::MeshID,std::string>( iterator->first, fname+":"+dbpath ) );
            }
            // Close the file
            DBClose ( FileHandle );
        }
        d_comm.barrier();
    }
    // Create the multimeshes
    std::map<AMP::Mesh::MeshID,siloMultiMeshData> multimeshes;
    std::map<AMP::Mesh::MeshID,AMP::Mesh::Mesh::shared_ptr>::iterator iterator;
    for (iterator=d_multiMeshes.begin(); iterator!=d_multiMeshes.end(); iterator++) {
        AMP::Mesh::MeshID id = iterator->first;
        AMP::Mesh::Mesh::shared_ptr mesh = iterator->second;
        std::vector<AMP::Mesh::MeshID> base_ids = mesh->getBaseMeshIDs();
        siloMultiMeshData multimesh;
        multimesh.id = id;
        multimesh.name = mesh->getName();
        for (size_t i=0; i<base_ids.size(); i++) {
            std::map<AMP::Mesh::MeshID,std::string>::iterator it = pathnames.find(base_ids[i]);
            if ( it != pathnames.end() )
                multimesh.paths.push_back(it->second);
        }
        if ( multimesh.paths.size()>0 ) 
        multimeshes.insert( std::pair<AMP::Mesh::MeshID,siloMultiMeshData>(id,multimesh) );
    }
    // Add the whole mesh
    siloMultiMeshData wholemesh;
    wholemesh.id = AMP::Mesh::MeshID(-1,0);
    wholemesh.name = "whole_mesh";
    std::map<AMP::Mesh::MeshID,std::string>::iterator iterator2;
    for (iterator2=pathnames.begin(); iterator2!=pathnames.end(); iterator2++)
        wholemesh.paths.push_back(iterator2->second);
    multimeshes.insert( std::pair<AMP::Mesh::MeshID,siloMultiMeshData>(wholemesh.id,wholemesh) );
    // Gather the results
    syncMultiMeshData( multimeshes );
    // Write the multimeshes
    if ( d_comm.getRank()==0 ) {
        FileHandle = DBOpen ( fname.c_str(), DB_HDF5, DB_APPEND );
        std::map<AMP::Mesh::MeshID,siloMultiMeshData>::iterator it;
        for (it=multimeshes.begin(); it!=multimeshes.end(); it++) {
            siloMultiMeshData data = it->second;
            char **meshnames = new char*[data.paths.size()];
            int *meshtypes = new int[data.paths.size()];
            for (size_t i=0; i<data.paths.size(); i++) {
                meshnames[i] = (char*) data.paths[i].c_str();
                meshtypes[i] = DB_UCDMESH;
            }
            DBPutMultimesh( FileHandle, data.name.c_str(), data.paths.size(), meshnames, meshtypes, NULL );
            delete [] meshnames;
            delete [] meshtypes;
        }
        DBClose ( FileHandle );
    }
}


/************************************************************
* Function to register a mesh with silo                     *
************************************************************/
void SiloIO::registerMesh( AMP::Mesh::Mesh::shared_ptr mesh )
{ 
    std::vector<AMP::Mesh::MeshID> ids = mesh->getAllMeshIDs();
    if ( ids.size()==1 ) {
        // We are dealing with a single mesh
        d_baseMeshes.insert( std::pair<AMP::Mesh::MeshID,AMP::Mesh::Mesh::shared_ptr>(mesh->meshID(),mesh) );
    } else {
        // We are dealining with a multimesh, register the current meshes and all submeshes
        for (size_t i=0; i<ids.size(); i++) {
            if ( ids[i] == mesh->meshID() ) {
                d_multiMeshes.insert( std::pair<AMP::Mesh::MeshID,AMP::Mesh::Mesh::shared_ptr>(mesh->meshID(),mesh) );
            } else {
                AMP::Mesh::Mesh::shared_ptr mesh2 = mesh->Subset(ids[i]);
                if ( mesh2.get() != NULL )
                    registerMesh( mesh2 );
            }
        }
    }
}


/************************************************************
* Function to register a vector with silo                   *
************************************************************/
#ifdef USE_AMP_VECTORS
void SiloIO::registerVector( AMP::LinearAlgebra::Vector::shared_ptr vec, const std::string &s )
{ 
    AMP_ERROR("registerVector is not implimented yet");
}
#endif


/************************************************************
* Function to write a mesh                                  *
************************************************************/
std::string SiloIO::writeMesh( DBfile *FileHandle, AMP::Mesh::Mesh::shared_ptr mesh )
{ 
    int dim = mesh->getDim();
    // Get the zone (element) lists
    AMP::Mesh::MeshIterator elem_iterator = mesh->getIterator(mesh->getGeomType(),0);
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
    DBSetDir( FileHandle, "All_Meshes" );
    if ( mesh->getComm().getRank()==0 )
        DBMkDir( FileHandle, mesh->getName().c_str() );
    DBSetDir( FileHandle, mesh->getName().c_str() );
    std::stringstream  stream;
    stream << mesh->getName() << "_" << mesh->getComm().getRank()+1;
    std::string meshName = stream.str();
    std::string zoneName = meshName + "_zone";
    AMP::Mesh::MeshIterator  element_iterator = mesh->getIterator(mesh->getGeomType(),0);
    int num_elems = (int) element_iterator.size();
    DBPutZonelist2( FileHandle, zoneName.c_str(), num_elems, dim, 
        &nodelist[0], nodelist.size(), 0, 0, 0, &shapetype, &shapesize, &shapecnt, 1, 0 );
    // Write the mesh
    DBPutUcdmesh( FileHandle, meshName.c_str(), mesh->getDim(),
        NULL, coord, node_iterator.size(), nodelist.size(),
        zoneName.c_str(), 0, DB_DOUBLE, 0 );
    DBSetDir( FileHandle, "/" );
    // Delete temporary memory
    for (int i=0; i<dim; i++)
        delete [] coord[i];
    return "All_Meshes/"+mesh->getName()+"/"+meshName;
}


/************************************************************
* Function to syncronize the multimesh data                 *
************************************************************/
void SiloIO::syncMultiMeshData( std::map<AMP::Mesh::MeshID,siloMultiMeshData> &data )
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
    // Combine the results
    data = std::map<AMP::Mesh::MeshID,siloMultiMeshData>();
    while ( meshdata.size() > 0 ) {
        siloMultiMeshData current = meshdata[0];
        std::vector<siloMultiMeshData>::iterator it = meshdata.begin();
        it = meshdata.erase(it);
        while ( it != meshdata.end() ) {
            siloMultiMeshData tmp = *it;
            if ( tmp.id==current.id ) {
                for (size_t i=0; i<tmp.paths.size(); i++)
                    current.paths.push_back(tmp.paths[i]);
                it = meshdata.erase(it);
            } else {
                it++;
            }
        }
        data.insert( std::pair<AMP::Mesh::MeshID,siloMultiMeshData>(current.id,current) );
    }
}


/************************************************************
* Functions for siloMultiMeshData                           *
************************************************************/
size_t SiloIO::siloMultiMeshData::size()
{
    size_t N_bytes = sizeof(AMP::Mesh::MeshID);
    N_bytes += sizeof(int);
    N_bytes += name.size();
    N_bytes += sizeof(int);
    N_bytes += sizeof(int)*paths.size();
    for (size_t i=0; i<paths.size(); i++)
        N_bytes += paths[i].size();
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
    // Store paths
    *((int*) &ptr[pos]) = (int) paths.size();
    pos += sizeof(int);
    for (size_t i=0; i<paths.size(); i++) {
        *((int*) &ptr[pos]) = (int) paths[i].size();
        pos += sizeof(int);
        paths[i].copy( (char*) &ptr[pos], paths[i].size(), 0 );
        pos += paths[i].size();
    }
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
    // Store paths
    size = *((int*) &ptr[pos]);
    pos += sizeof(int);
    data.paths.resize(size);
    for (size_t i=0; i<data.paths.size(); i++) {
        size = *((int*) &ptr[pos]);
        pos += sizeof(int);
        data.paths[i] = std::string( &ptr[pos], size );
        pos += size;
    }
    return data;
}


#else
void SiloIO::readFile( const std::string &fname ) { AMP_ERROR("SILO not configured"); }
void SiloIO::writeFile( const std::string &fname, size_t iteration_count ) { AMP_ERROR("SILO not configured"); }
void SiloIO::registerMesh( AMP::Mesh::Mesh::shared_ptr mesh, const std::string &s ) { AMP_ERROR("SILO not configured"); }
#ifdef USE_AMP_VECTORS
void SiloIO::registerVector( AMP::LinearAlgebra::Vector::shared_ptr vec, const std::string &s ) { AMP_ERROR("SILO not configured"); }
#endif

#endif


} // Mesh namespace
} // AMP namespace


