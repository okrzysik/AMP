#include "ampmesh/libmesh/libMesh.h"
#include "ampmesh/libmesh/libMeshIterator.h"
#include "utils/MemoryDatabase.h"
#include "utils/AMPManager.h"


// LibMesh include
#include "mesh.h"
#include "mesh_data.h"

namespace AMP {
namespace Mesh {


/********************************************************
* Constructors                                          *
********************************************************/
libMesh::libMesh( const MeshParameters::shared_ptr &params_in ):
    Mesh(params_in)
{
    // Check for valid inputs
    AMP_INSIST(params.get(),"Params must not be null");
    AMP_INSIST(comm!=AMP_MPI(AMP_COMM_NULL),"Communicator must be set");
    // Intialize libMesh (this needs to be moved oout of AMPManager)
    AMPManager::initializeLibmesh(comm);
    // Load the mesh
    if ( d_db.get() ) {
        // Database exists
        if ( d_db->keyExists("Filename") ) {
            // Read an existing mesh
            AMP_INSIST(d_db->keyExists("dim"),"Variable 'dim' must be set in the database");
            PhysicalDim = d_db->getInteger("dim");
            GeomDim = (GeomType) PhysicalDim;
            AMP_INSIST(PhysicalDim>0&&PhysicalDim<10,"Invalid dimension");
            // Create the libMesh objects
            d_libMesh = boost::shared_ptr< ::Mesh>( new ::Mesh(PhysicalDim) );
            d_libMeshData = boost::shared_ptr< ::MeshData>( new ::MeshData(*d_libMesh) );
            // Use libMesh to read the data
            d_libMesh->read(d_db->getString("Filename"));
        } else {
            AMP_ERROR("Unable to construct mesh with given parameters");
        }
    } else {
        AMP_ERROR("Error: params must contain a database object");
    }    
}


/********************************************************
* De-constructor                                        *
********************************************************/
libMesh::~libMesh()
{
}


/********************************************************
* Function to copy the mesh                             *
********************************************************/
Mesh libMesh::copy() const
{
    return libMesh(*this);
}


/********************************************************
* Return the number of elements                         *
********************************************************/
size_t libMesh::numLocalElements( GeomType &type )
{
    size_t N=0;
    if ( type==PhysicalDim ) {
        // This is a libMesh element
        N = d_libMesh->n_local_elem();
    } else if ( type==Vertex ) {
        // This is a libMesh node
        N = d_libMesh->n_local_nodes();
    } else {
        AMP_ERROR("Not implimented for this type");
    }
    return N;
}
size_t libMesh::numTotalElements( GeomType &type )
{
    size_t N=0;
    if ( type==PhysicalDim ) {
        // This is a libMesh element
        N = d_libMesh->parallel_n_elem();
    } else if ( type==Vertex ) {
        // This is a libMesh node
        N = d_libMesh->parallel_n_nodes();
    } else {
        AMP_ERROR("Not implimented for this type");
    }
    return N;
}


/********************************************************
* Return an iterator over the given geometric type      *
********************************************************/
MeshIterator libMesh::getIterator( GeomType &type, int gcw )
{
    libMeshIterator iterator;
    if ( type==PhysicalDim ) {
        // This is a libMesh element
        if ( gcw==0 ) {
            ::Mesh::element_iterator begin = d_libMesh->local_elements_begin();
            ::Mesh::element_iterator end   = d_libMesh->local_elements_end();
            iterator = libMeshIterator( 1, d_libMesh, gcw, (void*) &(begin), (void*) &(end), (void*) &(begin) );
        } else if ( gcw==0 ) {
            AMP_ERROR("Unfinished");
        } else {
            AMP_ERROR("Unsupported ghost cell width");
        }
    } else if ( type==Vertex ) {
        // This is a libMesh node
        if ( gcw==0 ) {
            ::Mesh::node_iterator begin = d_libMesh->local_nodes_begin();
            ::Mesh::node_iterator end   = d_libMesh->local_nodes_end();
            iterator = libMeshIterator( 0, d_libMesh, gcw, (void*) &(begin), (void*) &(end), (void*) &(begin) );
        } else if ( gcw==0 ) {
            AMP_ERROR("Unfinished");
        } else {
            AMP_ERROR("Unsupported ghost cell width");
        }
    } else {
        AMP_ERROR("Not implimented for this type");
    }
    MeshIterator test = iterator;
    return test;
}


} // Mesh namespace
} // AMP namespace

