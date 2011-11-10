#include "ampmesh/Mesh.h"
#include "utils/Utilities.h"

namespace AMP {
namespace Mesh {

static unsigned int nextLocalMeshID = 1; 

/********************************************************
* Constructors                                          *
********************************************************/
Mesh::Mesh( const MeshParameters::shared_ptr &params_in )
{
    params = params_in;
    GeomDim = null;
    comm = params->comm;
    d_db = params->d_db;
    AMP_ASSERT(comm!=AMP_MPI(AMP_COMM_NULL));
    setMeshID();
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
* Functions that aren't implimented for teh base class  *
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
MeshIterator Mesh::getSurfaceIterator( const GeomType ) 
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

