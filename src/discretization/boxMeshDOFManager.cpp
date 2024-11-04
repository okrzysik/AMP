#include "AMP/discretization/boxMeshDOFManager.h"
#include "AMP/mesh/structured/BoxMesh.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Utilities.h"
#include <set>
#include <vector>


namespace AMP::Discretization {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
boxMeshDOFManager::boxMeshDOFManager( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                      AMP::Mesh::GeomType type,
                                      int gcw,
                                      int DOFsPerObject )
    : simpleDOFManager(
          mesh, mesh->getIterator( type, 0 ), mesh->getIterator( type, gcw ), type, DOFsPerObject ),
      d_boxMesh( std::dynamic_pointer_cast<AMP::Mesh::BoxMesh>( mesh ) )
{
    AMP_INSIST( d_boxMesh, "Mesh must be a BoxMesh" );
}


/****************************************************************
 * Get the array size for the local variables                    *
 ****************************************************************/
ArraySize boxMeshDOFManager::getArraySize() const
{
    if ( d_DOFsPerElement == 0 || d_local_id.empty() )
        return {};
    auto box   = d_boxMesh->getLocalBlock( d_comm.getRank() );
    auto range = d_boxMesh->getIteratorRange( box, d_type, 0 );
    if ( range.size() != 1u )
        return {};
    auto first = range[0].first;
    auto last  = range[0].second;
    int ndim   = static_cast<int>( d_boxMesh->getGeomType() );
    size_t Nx  = last.index( 0 ) - first.index( 0 ) + 1;
    size_t Ny  = last.index( 1 ) - first.index( 1 ) + 1;
    size_t Nz  = last.index( 2 ) - first.index( 2 ) + 1;
    ArraySize size;
    if ( d_DOFsPerElement == 1 )
        size = ArraySize( { Nx, Ny, Nz }, ndim );
    else
        size = ArraySize( { (size_t) d_DOFsPerElement, Nx, Ny, Nz }, ndim + 1 );
    AMP_ASSERT( size.length() == d_DOFsPerElement * d_local_id.size() );
    return size;
}


} // namespace AMP::Discretization
