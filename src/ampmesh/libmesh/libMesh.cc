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
    // Intialize libMesh (this needs to be moved out of AMPManager)
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
            // Construct the neighbor information
            d_libMesh->find_neighbors();
        } else {
            AMP_ERROR("Unable to construct mesh with given parameters");
        }
    } else {
        AMP_ERROR("Error: params must contain a database object");
    }
    // Count the elements 
    n_local = std::vector<size_t>(PhysicalDim+1,0);
    n_global = std::vector<size_t>(PhysicalDim+1,0);
    n_ghost = std::vector<size_t>(PhysicalDim+1,0);
    for (int i=0; i<=(int)GeomDim; i++) {
        if ( i == (int) Vertex ) {
            // We are counting the nodes
            n_local[i] = d_libMesh->n_local_nodes();
            n_global[i] = d_libMesh->parallel_n_nodes();
            ::Mesh::node_iterator pos = d_libMesh->nodes_begin();
            ::Mesh::node_iterator end = d_libMesh->nodes_end();
            int N = 0;
            while ( pos != end ) {
                N++;
                ++pos;
            }
            n_ghost[i] = N - n_local[i];
        } else if ( i == (int) GeomDim ) {
            // We are counting the elements
            n_local[i] = d_libMesh->n_local_elem();
            n_global[i] = d_libMesh->parallel_n_elem();
            ::Mesh::element_iterator pos = d_libMesh->elements_begin();
            ::Mesh::element_iterator end = d_libMesh->elements_end();
            int N = 0;
            while ( pos != end ) {
                N++;
                ++pos;
            }
            n_ghost[i] = N - n_local[i];
        } else {
            // We are counting an intermediate element (not finished)
            n_local[i] = static_cast<size_t>(-1);
            n_global[i] = static_cast<size_t>(-1);
            n_ghost[i] = static_cast<size_t>(-1);
        }
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
size_t libMesh::numLocalElements( const GeomType type ) const
{
    if ( n_local[type] == static_cast<size_t>(-1) )
        AMP_ERROR("numLocalElements is not implimented for this type");
    return n_local[type];
}
size_t libMesh::numGlobalElements( const GeomType type ) const
{
    if ( n_global[type] == static_cast<size_t>(-1) )
        AMP_ERROR("numLocalElements is not implimented for this type");
    return n_global[type];
}
size_t libMesh::numGhostElements( const GeomType type, int gcw ) const
{
    if ( gcw == 0 )
        return 0;
    if ( gcw > 1 )
        AMP_ERROR("Libmesh only supports a ghost cell width of 1");
    if ( n_ghost[type] == static_cast<size_t>(-1) )
        AMP_ERROR("numLocalElements is not implimented for this type");
    return n_ghost[type];
}


/********************************************************
* Return an iterator over the given geometric type      *
********************************************************/
MeshIterator libMesh::getIterator( const GeomType type, const int gcw )
{
    libMeshIterator iterator;
    if ( type==PhysicalDim ) {
        // This is a libMesh element
        if ( gcw==0 ) {
            ::Mesh::element_iterator begin = d_libMesh->local_elements_begin();
            ::Mesh::element_iterator end   = d_libMesh->local_elements_end();
            iterator = libMeshIterator( 1, this, gcw, (void*) &(begin), (void*) &(end), (void*) &(begin) );
        } else if ( gcw==1 ) {
            ::Mesh::element_iterator begin = d_libMesh->elements_begin();
            ::Mesh::element_iterator end   = d_libMesh->elements_end();
            iterator = libMeshIterator( 1, this, gcw, (void*) &(begin), (void*) &(end), (void*) &(begin) );
        } else {
            AMP_ERROR("Unsupported ghost cell width");
        }
    } else if ( type==Vertex ) {
        // This is a libMesh node
        if ( gcw==0 ) {
            ::Mesh::node_iterator begin = d_libMesh->local_nodes_begin();
            ::Mesh::node_iterator end   = d_libMesh->local_nodes_end();
            iterator = libMeshIterator( 0, this, gcw, (void*) &(begin), (void*) &(end), (void*) &(begin) );
        } else if ( gcw==1 ) {
            ::Mesh::node_iterator begin = d_libMesh->nodes_begin();
            ::Mesh::node_iterator end   = d_libMesh->nodes_end();
            iterator = libMeshIterator( 0, this, gcw, (void*) &(begin), (void*) &(end), (void*) &(begin) );
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

