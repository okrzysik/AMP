#include "discretization/structuredFaceDOFManager.h"
#include "ampmesh/MultiIterator.h"
#include "ampmesh/MultiMesh.h"
#include "ampmesh/StructuredMeshHelper.h"
#include "discretization/MultiDOF_Manager.h"
#include "utils/Utilities.h"
#include <algorithm>
#include <set>
#include <vector>


namespace AMP {
namespace Discretization {


/****************************************************************
* Constructors                                                  *
****************************************************************/
DOFManager::shared_ptr structuredFaceDOFManager::create( AMP::shared_ptr<AMP::Mesh::Mesh> mesh,
                                                         int DOFsPerFace[3],
                                                         int gcw )
{
    if ( mesh.get() == NULL )
        return DOFManager::shared_ptr();
    if ( mesh->getGeomType() != AMP::Mesh::Volume || mesh->getDim() != 3 )
        AMP_ERROR( "The mesh must be a volume/3d mesh for structuredFaceDOFManager" );
    AMP::shared_ptr<structuredFaceDOFManager> manager( new structuredFaceDOFManager() );
    manager->d_comm = mesh->getComm();
    manager->d_mesh = mesh;
    manager->d_gcw  = gcw;
    for ( int i                   = 0; i < 3; i++ )
        manager->d_DOFsPerFace[i] = DOFsPerFace[i];
    manager->initialize();
    return manager;
}


/****************************************************************
* Deconstructor                                                 *
****************************************************************/
structuredFaceDOFManager::~structuredFaceDOFManager() {}


/****************************************************************
* Initialize the data                                           *
****************************************************************/
void structuredFaceDOFManager::initialize()
{
    // Create a sorted list of the local and remote types
    for ( int d = 0; d < 3; d++ ) {
        d_local_ids[d].resize( 0 );
        d_remote_ids[d].resize( 0 );
        if ( d_DOFsPerFace[d] == 0 )
            continue;
        const AMP::Mesh::MeshIterator localIterator =
            AMP::Mesh::StructuredMeshHelper::getFaceIterator( d_mesh, 0, d );
        const AMP::Mesh::MeshIterator ghostIterator =
            AMP::Mesh::StructuredMeshHelper::getFaceIterator( d_mesh, d_gcw, d );
        d_local_ids[d].resize( localIterator.size() );
        d_local_dofs[d].resize( localIterator.size() );
        AMP::Mesh::MeshIterator it = localIterator.begin();
        for ( size_t i = 0; i < localIterator.size(); ++i, ++it ) {
            d_local_ids[d][i]  = it->globalID();
            d_local_dofs[d][i] = i;
            AMP_ASSERT( d_local_ids[d][i].is_local() );
        }
        d_remote_ids[d].reserve( ghostIterator.size() - localIterator.size() );
        it = ghostIterator.begin();
        for ( size_t i = 0; i < ghostIterator.size(); ++i, ++it ) {
            AMP::Mesh::MeshElementID id = it->globalID();
            if ( !id.is_local() )
                d_remote_ids[d].push_back( id );
        }
        AMP_ASSERT( d_local_ids[d].size() == localIterator.size() );
        AMP_ASSERT( d_remote_ids[d].size() == ghostIterator.size() - localIterator.size() );
        // Sort the ids (keeping track of the original order)
        AMP::Utilities::quicksort( d_local_ids[d], d_local_dofs[d] );
        AMP::Utilities::quicksort( d_remote_ids[d] );
    }
    // Get the number of local elements per processor and the global number of DOFs
    size_t N_local = 0;
    for ( int d = 0; d < 3; d++ )
        N_local += d_local_ids[d].size() * d_DOFsPerFace[d];
    d_comm.sumScan<size_t>( &N_local, &d_end );
    d_begin  = d_end - N_local;
    d_global = d_comm.bcast( d_end, d_comm.getSize() - 1 );
    // Correct the local dof indicies
    size_t offset = d_begin;
    for ( int d = 0; d < 3; d++ ) {
        for ( size_t i         = 0; i < d_local_ids[d].size(); ++i )
            d_local_dofs[d][i] = offset + d_local_dofs[d][i] * d_DOFsPerFace[d];
        offset += d_local_ids[d].size() * d_DOFsPerFace[d];
    }
    // Determine the remote DOFs
    // Note: this must be done after d_local_id is set, d_begin and d_global are set, and remote_ids
    // must be sorted.
    for ( int d          = 0; d < 3; d++ )
        d_remote_dofs[d] = getRemoteDOF( d_remote_ids[d] );
}


/****************************************************************
* Get the DOFs for the element                                  *
****************************************************************/
void structuredFaceDOFManager::getDOFs( const AMP::Mesh::MeshElementID &id,
                                        std::vector<size_t> &dofs ) const
{
    dofs.resize( 0 );
    if ( id.type() != AMP::Mesh::Face )
        return;
    // Search for the dof locally
    for ( int d = 0; d < 3; d++ ) {
        if ( !d_local_ids[d].empty() ) {
            size_t index = AMP::Utilities::findfirst( d_local_ids[d], id );
            if ( index == d_local_ids[d].size() ) {
                index--;
            }
            if ( id == d_local_ids[d][index] ) {
                // The id was found
                dofs.resize( d_DOFsPerFace[d] );
                for ( int j = 0; j < d_DOFsPerFace[d]; j++ )
                    dofs[j] = d_local_dofs[d][index] + j;
                return;
            }
        }
    }
    // Search for the dof in the remote list
    for ( int d = 0; d < 3; d++ ) {
        if ( !d_remote_ids[d].empty() ) {
            size_t index = AMP::Utilities::findfirst( d_remote_ids[d], id );
            if ( index == d_remote_ids[d].size() ) {
                index--;
            }
            if ( id == d_remote_ids[d][index] ) {
                // The id was found
                dofs.resize( d_DOFsPerFace[d] );
                for ( int j = 0; j < d_DOFsPerFace[d]; j++ )
                    dofs[j] = d_remote_dofs[d][index] + j;
                return;
            }
        }
    }
}


/****************************************************************
* Get an entry over the mesh elements associated with the DOFs  *
****************************************************************/
AMP::Mesh::MeshIterator structuredFaceDOFManager::getIterator() const
{
    std::vector<AMP::Mesh::MeshIterator::shared_ptr> faces( 3 );
    for ( int i = 0; i < 3; i++ ) {
        if ( d_DOFsPerFace[i] > 0 )
            faces[i].reset( new AMP::Mesh::MeshIterator(
                AMP::Mesh::StructuredMeshHelper::getFaceIterator( d_mesh, 0, i ) ) );
    }
    return AMP::Mesh::MultiIterator( faces, 0 );
}


/****************************************************************
* Return the remote DOFs                                        *
****************************************************************/
std::vector<size_t> structuredFaceDOFManager::getRemoteDOFs() const
{
    // Create the list of remote DOFs
    size_t N_remote = 0;
    for ( size_t d = 0; d < 3; d++ )
        N_remote += d_remote_ids[d].size() * d_DOFsPerFace[d];
    std::vector<size_t> remote_DOFs( N_remote, (size_t) -1 );
    int i = 0;
    for ( size_t d = 0; d < 3; d++ ) {
        for ( size_t j = 0; j < d_remote_ids[d].size(); j++ ) {
            for ( int k = 0; k < d_DOFsPerFace[d]; k++ ) {
                remote_DOFs[i] = d_remote_dofs[d][j] + k;
                ++i;
            }
        }
    }
    AMP::Utilities::quicksort( remote_DOFs );
    return remote_DOFs;
}


/****************************************************************
* Return the row DOFs                                           *
****************************************************************/
std::vector<size_t> structuredFaceDOFManager::getRowDOFs( const AMP::Mesh::MeshElement &obj ) const
{
    if ( obj.elementType() != AMP::Mesh::Face )
        return std::vector<size_t>();
    // Get a list of all element ids that are part of the row
    // Only faces that share an element are part of the row
    std::vector<AMP::Mesh::MeshElement> parents =
        d_mesh->getElementParents( obj, AMP::Mesh::Volume );
    AMP_ASSERT( parents.size() == 1 || parents.size() == 2 );
    // Temporarily add neighbor elements
    size_t p_size = parents.size();
    for ( size_t i = 0; i < p_size; i++ ) {
        std::vector<AMP::Mesh::MeshElement::shared_ptr> neighbors = parents[i].getNeighbors();
        for ( size_t j = 0; j < neighbors.size(); j++ ) {
            if ( neighbors[j] != NULL )
                parents.push_back( *neighbors[j] );
        }
    }
    AMP::Utilities::unique( parents );
    // Get the face ids of interest
    std::vector<AMP::Mesh::MeshElementID> ids;
    ids.reserve( 6 * parents.size() );
    for ( size_t i = 0; i < parents.size(); i++ ) {
        std::vector<AMP::Mesh::MeshElement> children = parents[i].getElements( AMP::Mesh::Face );
        AMP_ASSERT( children.size() == 6 );
        for ( size_t j = 0; j < children.size(); j++ )
            ids.push_back( children[j].globalID() );
    }
    AMP::Utilities::unique( ids );
    // AMP_ASSERT(ids.size()==6||ids.size()==11);
    // Get all dofs for each element id
    int maxDOFsPerFace = 0;
    for ( int i = 0; i < 3; i++ )
        maxDOFsPerFace = std::max( maxDOFsPerFace, d_DOFsPerFace[i] );
    std::vector<size_t> dofs;
    dofs.reserve( ids.size() * maxDOFsPerFace );
    std::vector<size_t> dofs2( maxDOFsPerFace );
    for ( size_t i = 0; i < ids.size(); i++ ) {
        getDOFs( ids[i], dofs2 );
        for ( size_t j = 0; j < dofs2.size(); j++ )
            dofs.push_back( dofs2[j] );
    }
    // Sort the row dofs
    AMP::Utilities::quicksort( dofs );
    return dofs;
}


/****************************************************************
* Find the remote DOF given a set of mesh element IDs           *
* Note: for this function to work correctly, the remote ids     *
* must be sorted, and d_local_id must be set                    *
****************************************************************/
std::vector<size_t> structuredFaceDOFManager::getRemoteDOF(
    const std::vector<AMP::Mesh::MeshElementID> &remote_ids ) const
{
    if ( d_comm.sumReduce<size_t>( remote_ids.size() ) == 0 )
        return std::vector<size_t>(); // There are no remote DOFs
    // Get the set of mesh ids (must match on all processors)
    std::vector<AMP::Mesh::MeshID> meshIDs = d_mesh->getBaseMeshIDs();
    // Get the rank that will own each MeshElement on the current communicator
    std::vector<int> owner_rank( remote_ids.size(), -1 );
    for ( size_t it = 0; it < meshIDs.size(); it++ ) {
        // Get the mesh with the given meshID
        AMP::Mesh::MeshID meshID            = meshIDs[it];
        AMP::Mesh::Mesh::shared_ptr submesh = d_mesh->Subset( meshID );
        // Create a map from the rank of the submesh to the current mesh
        int rank_submesh = -1;
        int root_submesh = d_comm.getSize();
        int subcommSize  = -1;
        int myRank       = -1;
        if ( submesh.get() != NULL ) {
            AMP_MPI subcomm = submesh->getComm();
            rank_submesh    = subcomm.getRank();
            root_submesh    = d_comm.getRank();
            subcommSize     = subcomm.getSize();
            myRank          = subcomm.getRank();
        }
        root_submesh = d_comm.minReduce( root_submesh );
        if ( root_submesh == d_comm.getSize() )
            AMP_ERROR( "Not processors on the current comm exist on the submesh comm" );
        subcommSize = d_comm.bcast( subcommSize, root_submesh );
        std::vector<int> subrank( d_comm.getSize() );
        d_comm.allGather( rank_submesh, &subrank[0] );
        std::vector<int> rank_map( subcommSize, -1 );
        for ( size_t i = 0; i < subrank.size(); i++ ) {
            if ( subrank[i] != -1 )
                rank_map[subrank[i]] = i;
        }
        // Get the rank of the proccessor that will own each meshElement
        for ( size_t i = 0; i < remote_ids.size(); i++ ) {
            if ( remote_ids[i].meshID() == meshID ) {
                int subowner_rank = remote_ids[i].owner_rank();
                AMP_ASSERT( subowner_rank != myRank );
                AMP_ASSERT( rank_map[subowner_rank] != -1 );
                owner_rank[i] = rank_map[subowner_rank];
            }
        }
    }
    // Check that each element has a vaild owner rank
    int commSize = d_comm.getSize();
    for ( size_t i = 0; i < remote_ids.size(); i++ )
        AMP_ASSERT( owner_rank[i] >= 0 && owner_rank[i] < commSize );
    // Resort the remote ids according the the owner rank
    std::vector<AMP::Mesh::MeshElementID> remote_ids2 = remote_ids;
    AMP::Utilities::quicksort( owner_rank, remote_ids2 );
    // Determine the send count and displacements for each processor
    std::vector<int> send_cnt( d_comm.getSize(), 0 );
    std::vector<int> send_disp( d_comm.getSize(), 0 );
    int rank     = 0;
    size_t start = 0;
    size_t index = 0;
    while ( index < owner_rank.size() ) {
        if ( owner_rank[index] < rank ) {
            AMP_ERROR( "This should not occur" );
        } else if ( owner_rank[index] == rank ) {
            // Move to the next element
            index++;
        } else {
            // Store the number of elements with the given rank, and move to the next rank
            send_disp[rank] = start;
            send_cnt[rank]  = index - start;
            start           = index;
            rank++;
        }
    }
    send_disp[rank] = start;
    send_cnt[rank]  = index - start;
    // Perform an allToAll to send the remote ids for DOF identification
    std::vector<int> recv_cnt( d_comm.getSize() );
    d_comm.allToAll<int>( 1, &send_cnt[0], &recv_cnt[0] );
    std::vector<int> recv_disp( d_comm.getSize() );
    size_t tot_size = recv_cnt[0];
    recv_disp[0]    = 0;
    for ( int i = 1; i < d_comm.getSize(); i++ ) {
        tot_size += recv_cnt[i];
        recv_disp[i] = recv_disp[i - 1] + recv_cnt[i - 1];
    }
    std::vector<AMP::Mesh::MeshElementID> recv_id( tot_size + 1 );
    AMP::Mesh::MeshElementID *send_buffer = NULL;
    if ( !remote_ids2.empty() )
        send_buffer = &remote_ids2[0];
    size_t N        = d_comm.allToAll<AMP::Mesh::MeshElementID>(
        send_buffer, &send_cnt[0], &send_disp[0], &recv_id[0], &recv_cnt[0], &recv_disp[0], true );
    AMP_INSIST( N == tot_size, "Unexpected recieve size" );
    recv_id.resize( tot_size );
    // Determine the DOF for each recieved id
    std::vector<size_t> recieved_DOF( tot_size ), dofs;
    for ( size_t i = 0; i < tot_size; i++ ) {
        getDOFs( recv_id[i], dofs );
        AMP_ASSERT( dofs.size() > 0 );
        for ( size_t j = 1; j < dofs.size(); j++ )
            AMP_ASSERT( dofs[j] == dofs[j - 1] + 1 );
        recieved_DOF[i] = dofs[0];
    }
    // Send the DOFs back to the original processor
    std::vector<size_t> remote_dof;
    remote_dof.resize( remote_ids2.size() + 1, static_cast<size_t>( -1 ) );
    size_t *send_buffer_DOFs = NULL;
    if ( tot_size > 0 )
        send_buffer_DOFs = &recieved_DOF[0];
    N                    = d_comm.allToAll<size_t>( send_buffer_DOFs,
                                 &recv_cnt[0],
                                 &recv_disp[0],
                                 &remote_dof[0],
                                 &send_cnt[0],
                                 &send_disp[0],
                                 true );
    AMP_INSIST( N == remote_ids2.size(), "Unexpected recieve size" );
    remote_dof.resize( remote_ids2.size() );
    // Sort the dofs back to the original order for the remote_ids
    AMP::Utilities::quicksort( remote_ids2, remote_dof );
    for ( size_t i = 0; i < remote_ids.size(); i++ )
        AMP_ASSERT( remote_ids[i] == remote_ids2[i] );
    return remote_dof;
}
}
}
