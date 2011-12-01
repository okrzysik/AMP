#include "ampmesh/Mesh.h"
#include "utils/Utilities.h"

#include "ampmesh/MultiMesh.h"
#ifdef USE_LIBMESH
#include "ampmesh/libmesh/libMesh.h"
#endif

namespace AMP {
namespace Mesh {

static unsigned int nextLocalMeshID = 1; 

/********************************************************
* Constructors                                          *
********************************************************/
Mesh::Mesh( const MeshParameters::shared_ptr &params_in )
{
    // Set the base properties
    AMP_ASSERT(sizeof(MeshElementID)==16);
    params = params_in;
    GeomDim = null;
    comm = params->comm;
    d_db = params->d_db;
    AMP_ASSERT(comm!=AMP_MPI(AMP_COMM_NULL));
    setMeshID();
    d_name = "NULL";
}
Mesh::Mesh( const Mesh::shared_ptr &old_mesh )
{
    AMP_ERROR("Not Implimented Yet");
}


/********************************************************
* De-constructor                                        *
********************************************************/
Mesh::~Mesh()
{
}


/********************************************************
* Assignment operator                                   *
********************************************************/
Mesh Mesh::operator=(const Mesh& rhs)
{
    return rhs.copy();
}
Mesh Mesh::copy() const
{
    return Mesh(*this);
}


/********************************************************
* Create a mesh from the input database                 *
********************************************************/
boost::shared_ptr<AMP::Mesh::Mesh> Mesh::buildMesh( const MeshParameters::shared_ptr &params )
{
    boost::shared_ptr<AMP::Database> database = params->d_db;
    AMP_ASSERT(database!=NULL);
    AMP_INSIST(database->keyExists("MeshType"),"MeshType must exist in input database");
    AMP_INSIST(database->keyExists("MeshName"),"MeshName must exist in input database");
    std::string MeshType = database->getString("MeshType");
    std::string MeshName = database->getString("MeshName");
    boost::shared_ptr<AMP::Mesh::Mesh> mesh;
    if ( MeshType == std::string("Multimesh") ) {
        // The mesh is a multimesh
        mesh = boost::shared_ptr<AMP::Mesh::MultiMesh>(new AMP::Mesh::MultiMesh(params) );
    } else if ( MeshType == std::string("libMesh") ) {
        // The mesh is a libmesh mesh
        #ifdef USE_LIBMESH
            mesh = boost::shared_ptr<AMP::Mesh::libMesh>(new AMP::Mesh::libMesh(params) );
        #else
            AMP_ERROR("AMP was compiled without support for libMesh");
        #endif
    } else {
        // Unknown mesh type
        AMP_ERROR( std::string("Unknown mesh type (") + MeshType + std::string(")") );
    }
    mesh->setName(MeshName);
    return mesh;
}


/********************************************************
* Estimate the mesh size                                *
********************************************************/
size_t Mesh::estimateMeshSize( const MeshParameters::shared_ptr &params )
{
    boost::shared_ptr<AMP::Database> database = params->d_db;
    AMP_ASSERT(database!=NULL);
    size_t meshSize = 0;
    // This is being called through the base class, call the appropriate function
    AMP_INSIST(database->keyExists("MeshType"),"MeshType must exist in input database");
    std::string MeshType = database->getString("MeshType");
    boost::shared_ptr<AMP::Mesh::Mesh> mesh;
    if ( MeshType == std::string("Multimesh") ) {
        // The mesh is a multimesh
        meshSize = AMP::Mesh::MultiMesh::estimateMeshSize(params);
    } else if ( MeshType == std::string("libMesh") ) {
        // The mesh is a libmesh mesh
        #ifdef USE_LIBMESH
            meshSize = AMP::Mesh::libMesh::estimateMeshSize(params);
        #else
            AMP_ERROR("AMP was compiled without support for libMesh");
        #endif
    } else if ( database->keyExists("NumberOfElements") ) {
        int NumberOfElements = database->getInteger("NumberOfElements");
        meshSize = NumberOfElements;
    } else {
        // Unknown mesh type
        AMP_ERROR( "Unknown mesh type and NumberOfElements does not exist in database" );
    }
    return meshSize;
}


/********************************************************
* Function to set the mesh ID                           *
* This function will create a unique ID for every mesh. *
* To accomplish this goal, the ID will consist of the   *
* rank of the root processor (from the global comm),    *
* and the number of meshes created by that processor.   *
********************************************************/
void Mesh::setMeshID( )
{
    if ( comm.getRank()==0 ) {
        AMP_MPI globalComm = AMP_MPI(AMP_COMM_WORLD);
        // I am the root processor, create the mesh ID
        if ( sizeof(size_t) == 8 ) {
            // 64-bit integer, use the first 32-bits for the processor
            // and the second 32-bits for the local id
            d_meshID = globalComm.getRank();
            d_meshID <<= 32;
            d_meshID += (size_t) nextLocalMeshID;
            nextLocalMeshID++;
        } else {
            // 32-bit integer, use the first 16-bits for the processor
            // and the second 16-bits for the local id
            if ( globalComm.getSize() > 0xFFFF )
                AMP_ERROR("Too many processors for meshID");
            d_meshID = globalComm.getRank();
            d_meshID <<= 16;
            d_meshID += (size_t) nextLocalMeshID;
            nextLocalMeshID++;
            if ( nextLocalMeshID > 0xFFFF )
                AMP_ERROR("Ran out of unique mesh IDs");
        }
    }
    // Broadcast the meshID to all processors
    d_meshID = comm.bcast(d_meshID,0);
}


/********************************************************
* Function to return the mesh with the given ID         *
********************************************************/
boost::shared_ptr<Mesh>  Mesh::Subset( size_t meshID ) {
    if ( d_meshID==meshID ) 
        return shared_from_this();
    else
        return boost::shared_ptr<Mesh>();
}


/********************************************************
* Function to return the mesh with the given name       *
********************************************************/
boost::shared_ptr<Mesh>  Mesh::Subset( std::string name ) {
    if ( d_name==name ) 
        return shared_from_this();
    else
        return boost::shared_ptr<Mesh>();
}


/********************************************************
* Functions that aren't implimented for the base class  *
********************************************************/
boost::shared_ptr<Mesh> Mesh::Subset( MeshIterator::shared_ptr & )
{
    AMP_ERROR("Not Implimented Yet");
    return boost::shared_ptr<Mesh>();
}
boost::shared_ptr<Mesh> Mesh::Subset( Mesh & )
{
    AMP_ERROR("Not Implimented Yet");
    return boost::shared_ptr<Mesh>();
}
MeshIterator Mesh::getIterator( const GeomType, const int )
{
    AMP_ERROR("Not Implimented Yet");
    return MeshIterator();
}
MeshIterator Mesh::getSurfaceIterator( const GeomType, const int ) 
{
    AMP_ERROR("Not Implimented Yet");
    return MeshIterator();
}
std::vector<int> Mesh::getIDSets ( )
{
    AMP_ERROR("Not Implimented Yet");
    return std::vector<int>();
}
MeshIterator Mesh::getIDsetIterator ( const GeomType, const int, const int )
{
    AMP_ERROR("Not Implimented Yet");
    return MeshIterator();
}
MeshIterator Mesh::getIterator( SetOP &, MeshIterator::shared_ptr &, MeshIterator::shared_ptr &)
{
    AMP_ERROR("Not Implimented Yet");
    return MeshIterator();
}
size_t Mesh::numLocalElements( const GeomType type ) const
{
    AMP_ERROR("Not Implimented Yet");
    return 0;
}
size_t Mesh::numGlobalElements( const GeomType type ) const
{
    AMP_ERROR("Not Implimented Yet");
    return 0;
}
size_t Mesh::numGhostElements( const GeomType type, int gcw ) const
{
    AMP_ERROR("Not Implimented Yet");
    return 0;
}


} // Mesh namespace
} // AMP namespace

