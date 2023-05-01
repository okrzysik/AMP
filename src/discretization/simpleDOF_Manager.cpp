#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/discretization/MultiDOF_Manager.h"
#include "AMP/discretization/boxMeshDOFManager.h"
#include "AMP/mesh/MultiMesh.h"
#include "AMP/mesh/structured/BoxMesh.h"
#include "AMP/utils/AMP_MPI.I"
#include "AMP/utils/Utilities.h"

#include <set>
#include <vector>


namespace AMP::Discretization {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
std::shared_ptr<DOFManager> simpleDOFManager::create( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                                      AMP::Mesh::GeomType type,
                                                      int gcw,
                                                      int DOFsPerObject,
                                                      bool split )
{
    if ( !mesh )
        return std::shared_ptr<DOFManager>();
    if ( split && std::dynamic_pointer_cast<AMP::Mesh::MultiMesh>( mesh ) ) {
        // We want to split the DOFs by the mesh
        auto meshIDs = mesh->getLocalBaseMeshIDs();
        std::vector<std::shared_ptr<DOFManager>> managers;
        for ( auto &meshID : meshIDs ) {
            auto subMesh = mesh->Subset( meshID );
            if ( subMesh )
                managers.push_back( create( subMesh, type, gcw, DOFsPerObject, false ) );
        }
        auto rtn = std::make_shared<multiDOFManager>( mesh->getComm(), managers );
        return rtn;
    }
    // Check if the mesh is a BoxMesh
    if ( std::dynamic_pointer_cast<AMP::Mesh::BoxMesh>( mesh ) ) {
        return std::make_shared<boxMeshDOFManager>( mesh, type, gcw, DOFsPerObject );
    }
    // Create a simpleDOFManager
    auto it1 = mesh->getIterator( type, 0 );
    auto it2 = mesh->getIterator( type, gcw );
    return std::make_shared<simpleDOFManager>( mesh, it1, it2, type, DOFsPerObject );
}
std::shared_ptr<DOFManager> simpleDOFManager::create( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                                      const AMP::Mesh::MeshIterator &it1,
                                                      const AMP::Mesh::MeshIterator &it2,
                                                      int DOFsPerObject )
{
    // Check the iterators
    auto intersection = AMP::Mesh::Mesh::getIterator( AMP::Mesh::SetOP::Intersection, it1, it2 );
    AMP_INSIST( intersection.size() == it2.size(), "it1 must include it2" );
    auto tmp = it2.begin();
    for ( size_t i = 0; i < tmp.size(); i++ ) {
        auto id = tmp->globalID();
        AMP_INSIST( id.is_local(), "it2 may not contain any ghost elements" );
        ++tmp;
    }
    tmp       = it1.begin();
    auto type = tmp->globalID().type();
    for ( size_t i = 0; i < tmp.size(); i++ ) {
        auto id = tmp->globalID();
        AMP_INSIST( id.type() == type, "All elements in the iterator must be the same type" );
        ++tmp;
    }
    // Create the simpleDOFManager
    return std::make_shared<simpleDOFManager>( mesh, it2, it1, type, DOFsPerObject );
}
std::shared_ptr<DOFManager> simpleDOFManager::create( const AMP::Mesh::MeshIterator &it,
                                                      int DOFsPerObject )
{
    // Check the iterator
    auto tmp  = it.begin();
    auto type = tmp->globalID().type();
    for ( size_t i = 0; i < tmp.size(); i++ ) {
        auto id = tmp->globalID();
        AMP_INSIST( id.type() == type, "All elements in the iterator must be the same type" );
        ++tmp;
    }
    // Create the simpleDOFManager
    return std::make_shared<simpleDOFManager>( nullptr, it, it, type, DOFsPerObject );
}
simpleDOFManager::simpleDOFManager( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                    const AMP::Mesh::MeshIterator &it1,
                                    const AMP::Mesh::MeshIterator &it2,
                                    AMP::Mesh::GeomType type,
                                    int DOFsPerObject )
    : d_type( type ),
      d_DOFsPerElement( DOFsPerObject ),
      d_mesh( mesh ),
      d_localIterator( it1 ),
      d_ghostIterator( it2 )
{
    AMP_ASSERT( d_ghostIterator.size() >= d_localIterator.size() );
    d_comm = AMP_MPI( AMP_COMM_SELF );
    if ( mesh )
        d_comm = mesh->getComm();
    initialize();
}


/****************************************************************
 * Deconstructor                                                 *
 ****************************************************************/
simpleDOFManager::~simpleDOFManager() = default;


/****************************************************************
 * Initialize the data                                           *
 ****************************************************************/
void simpleDOFManager::initialize()
{
    // Get the mesh ids
    if ( d_mesh != nullptr ) {
        d_meshID     = d_mesh->meshID();
        d_isBaseMesh = d_mesh->isBaseMesh();
        // Get the list of global mesh ids (use communication on this->comm)
        d_baseMeshIDs = d_mesh->getBaseMeshIDs();
        std::set<AMP::Mesh::MeshID> set( d_baseMeshIDs.begin(), d_baseMeshIDs.end() );
        d_comm.setGather( set );
        d_baseMeshIDs = std::vector<AMP::Mesh::MeshID>( set.begin(), set.end() );
    } else {
        d_meshID     = AMP::Mesh::MeshID();
        d_isBaseMesh = false;
        d_baseMeshIDs.clear();
    }
    // Create a sorted list of the local and remote types
    d_local_id.clear();
    d_remote_id.clear();
    d_local_id.reserve( d_localIterator.size() );
    d_remote_id.reserve( d_ghostIterator.size() - d_localIterator.size() );
    for ( const auto &elem : d_ghostIterator ) {
        auto id = elem.globalID();
        if ( id.is_local() )
            d_local_id.push_back( id );
        else
            d_remote_id.push_back( id );
    }
    AMP_ASSERT( d_local_id.size() == d_localIterator.size() );
    AMP_ASSERT( d_remote_id.size() == d_ghostIterator.size() - d_localIterator.size() );
    // Sort the elements (they will be sorted by the meshID, then the rank on the
    // comm of the given mesh, then the element type, and finally the local id)
    AMP::Utilities::quicksort( d_local_id );
    AMP::Utilities::quicksort( d_remote_id );
    // Get the number of local elements per processor and the global number of DOFs
    size_t N_local = d_local_id.size() * d_DOFsPerElement;
    d_comm.sumScan<size_t>( &N_local, &d_end, 1 );
    d_begin  = d_end - N_local;
    d_global = d_comm.bcast( d_end, d_comm.getSize() - 1 );
    // Determine the remote DOFs (assuming 1 DOF per node)
    // Note: this must be done after d_local_id is set, d_begin and d_global are set, and remote_ids
    // must be sorted.
    d_remote_dof = getRemoteDOF( d_remote_id );
    AMP_ASSERT( d_remote_dof.size() == d_remote_id.size() );
}


/****************************************************************
 * Subset the DOF manager                                        *
 ****************************************************************/
std::shared_ptr<DOFManager> simpleDOFManager::subset( const std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                                      bool useMeshComm )
{

    // Check if we are dealing with a single mesh for both the internal and desired mesh
    if ( mesh->meshID() == d_meshID ) {
        // The mesh IDs match
        return shared_from_this();
    } else if ( mesh->isBaseMesh() && d_isBaseMesh ) {
        // Both meshes are base meshes and the ids do not match
        return std::shared_ptr<DOFManager>();
    } else if ( d_baseMeshIDs.size() == 1 && mesh->isBaseMesh() ) {
        // The subsetting mesh is a base mesh and we only contain one mesh
        if ( d_baseMeshIDs[0] == mesh->meshID() )
            return shared_from_this();
        else
            return std::shared_ptr<DOFManager>();
    }
    // Check if the desired mesh is a multimesh that contains the current mesh
    if ( std::dynamic_pointer_cast<const AMP::Mesh::MultiMesh>( mesh ) != nullptr ) {
        auto multimesh = std::dynamic_pointer_cast<const AMP::Mesh::MultiMesh>( mesh );
        auto list      = multimesh->getMeshes();
        std::shared_ptr<DOFManager> subset_DOFs;
        bool found_local = false;
        for ( auto &elem : list ) {
            if ( d_meshID == elem->meshID() )
                found_local = true;
        }
        AMP_MPI comm( AMP_COMM_NULL );
        if ( useMeshComm ) {
            comm = AMP_MPI::intersect( d_comm, mesh->getComm() );
        } else {
            comm = d_comm;
        }
        found_local = comm.allReduce( found_local );
        comm.reset();
        if ( found_local )
            return shared_from_this();
    }
    // We were not able to use an efficient subset, use the generic base function
    return DOFManager::subset( mesh, useMeshComm );
}


/****************************************************************
 * Get the DOFs for the element                                  *
 ****************************************************************/
inline void simpleDOFManager::appendDOFs( const AMP::Mesh::MeshElementID &id,
                                          std::vector<size_t> &dofs ) const
{
    // Search for the dof locally
    if ( !d_local_id.empty() ) {
        size_t index = AMP::Utilities::findfirst( d_local_id, id );
        if ( index == d_local_id.size() )
            index--;
        if ( id == d_local_id[index] ) {
            // The id was found
            for ( int j = 0; j < d_DOFsPerElement; j++ )
                dofs.emplace_back( index * d_DOFsPerElement + d_begin + j );
            return;
        }
    }
    // Search for the dof in the remote list
    if ( !d_remote_id.empty() ) {
        size_t index = AMP::Utilities::findfirst( d_remote_id, id );
        if ( index == d_remote_id.size() )
            index--;
        if ( id == d_remote_id[index] ) {
            // The id was found
            for ( int j = 0; j < d_DOFsPerElement; j++ )
                dofs.emplace_back( d_remote_dof[index] * d_DOFsPerElement + j );
        }
    }
}
void simpleDOFManager::getDOFs( const std::vector<AMP::Mesh::MeshElementID> &ids,
                                std::vector<size_t> &dofs ) const
{
    dofs.resize( 0 );
    dofs.reserve( ids.size() * d_DOFsPerElement );
    for ( auto id : ids )
        appendDOFs( id, dofs );
}
void simpleDOFManager::getDOFs( const AMP::Mesh::MeshElementID &id,
                                std::vector<size_t> &dofs ) const
{
    dofs.resize( 0 );
    dofs.reserve( d_DOFsPerElement );
    appendDOFs( id, dofs );
}


/****************************************************************
 * Get the element ID give a dof                                 *
 ****************************************************************/
AMP::Mesh::MeshElement simpleDOFManager::getElement( size_t dof ) const
{
    AMP::Mesh::MeshElementID id;
    if ( dof >= d_begin && dof < d_end ) {
        // We are searching for a local dof
        id = d_local_id[( dof - d_begin ) / d_DOFsPerElement];
    }
    const size_t dof2 = dof / d_DOFsPerElement;
    for ( size_t i = 0; i < d_remote_id.size(); i++ ) {
        if ( d_remote_dof[i] == dof2 )
            id = d_remote_id[i];
    }
    if ( id.isNull() )
        return AMP::Mesh::MeshElement();
    return d_mesh->getElement( id );
}


/****************************************************************
 * Get an entry over the mesh elements associated with the DOFs  *
 ****************************************************************/
AMP::Mesh::MeshIterator simpleDOFManager::getIterator() const { return d_localIterator.begin(); }


/****************************************************************
 * Return the remote DOFs for a vector                           *
 ****************************************************************/
std::vector<size_t> simpleDOFManager::getRemoteDOFs() const
{
    // Create the list of remote DOFs
    size_t N = d_remote_id.size();
    std::vector<size_t> remote_DOFs( N * d_DOFsPerElement, (size_t) -1 );
    for ( size_t i = 0, k = 0; i < N; i++ ) {
        for ( int j = 0; j < d_DOFsPerElement; j++, k++ )
            remote_DOFs[k] = d_remote_dof[i] * d_DOFsPerElement + j;
    }
    AMP::Utilities::quicksort( remote_DOFs );
    return remote_DOFs;
}


/****************************************************************
 * Return the row DOFs                                           *
 ****************************************************************/
std::vector<size_t> simpleDOFManager::getRowDOFs( const AMP::Mesh::MeshElement &obj ) const
{
    // Get a list of all element ids that are part of the row
    auto meshType = d_mesh->getGeomType();
    auto objType  = obj.elementType();
    std::vector<AMP::Mesh::MeshElementID> ids;
    if ( objType == d_type && ( objType == AMP::Mesh::GeomType::Vertex || objType == meshType ) ) {
        // Use the getNeighbors function to get the neighbors of the current element
        auto neighbor_elements = obj.getNeighbors();
        ids.reserve( neighbor_elements.size() + 1 );
        ids.push_back( obj.globalID() );
        for ( auto &neighbor_element : neighbor_elements ) {
            if ( neighbor_element )
                ids.push_back( neighbor_element->globalID() );
        }
    } else if ( objType == d_type ) {
        // We need to use the mesh to get the connectivity of the elements of the same type
        auto parents = d_mesh->getElementParents( obj, meshType );
        for ( auto &parent : parents ) {
            auto children = parent.getElements( objType );
            ids.reserve( ids.size() + children.size() );
            for ( auto &elem : children )
                ids.push_back( elem.globalID() );
        }
        AMP::Utilities::unique( ids );
    } else if ( objType > d_type ) {
        // The desired element type is < the current element type, use getElements
        auto children = obj.getElements( d_type );
        for ( auto &elem : children )
            ids.push_back( elem.globalID() );
    } else if ( objType < d_type ) {
        // The desired element type is < the current element type, use getElementParents
        auto parents = d_mesh->getElementParents( obj, meshType );
        for ( auto &parent : parents )
            ids.push_back( parent.globalID() );
    } else {
        AMP_ERROR( "Internal error" );
    }
    // Get all dofs for each element id
    std::vector<size_t> dofs;
    dofs.reserve( ids.size() * d_DOFsPerElement );
    std::vector<size_t> dofs2( d_DOFsPerElement );
    for ( auto &id : ids ) {
        getDOFs( id, dofs2 );
        for ( auto &elem : dofs2 )
            dofs.push_back( elem );
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
std::vector<size_t>
simpleDOFManager::getRemoteDOF( std::vector<AMP::Mesh::MeshElementID> remote_ids ) const
{
    if ( d_comm.getSize() == 1 )
        return std::vector<size_t>(); // There are no remote DOFs
    // Get the set of mesh ids (must match on all processors)
    AMP_MPI comm = d_mesh->getComm();
    std::set<AMP::Mesh::MeshID> meshIDs;
    for ( auto &remote_id : remote_ids )
        meshIDs.insert( remote_id.meshID() );
    std::vector<AMP::Mesh::MeshID> tmpLocalIDs( meshIDs.begin(), meshIDs.end() );
    auto N = (int) comm.sumReduce<size_t>( tmpLocalIDs.size() );
    if ( N == 0 ) {
        // Nobody has any remote ids to identify
        return std::vector<size_t>();
    }
    AMP::Mesh::MeshID *send_ptr = nullptr;
    if ( !tmpLocalIDs.empty() )
        send_ptr = &tmpLocalIDs[0];
    std::vector<AMP::Mesh::MeshID> tmpGlobalIDs( N );
    int N_recv =
        comm.allGather<AMP::Mesh::MeshID>( send_ptr, tmpLocalIDs.size(), &tmpGlobalIDs[0] );
    AMP_ASSERT( N_recv == N );
    for ( auto &tmpGlobalID : tmpGlobalIDs )
        meshIDs.insert( tmpGlobalID );
    // Get the rank that will own each MeshElement on the current communicator
    std::vector<int> owner_rank( remote_ids.size(), -1 );
    for ( auto meshID : meshIDs ) {
        // Get the mesh with the given meshID
        auto submesh = d_mesh->Subset( meshID );
        // Create a map from the rank of the submesh to the current mesh
        int rank_submesh = -1;
        int root_submesh = comm.getSize();
        int subcommSize  = -1;
        if ( submesh ) {
            AMP_MPI subcomm = submesh->getComm();
            rank_submesh    = subcomm.getRank();
            root_submesh    = comm.getRank();
            subcommSize     = subcomm.getSize();
        }
        root_submesh = comm.minReduce( root_submesh );
        if ( root_submesh == comm.getSize() )
            AMP_ERROR( "Not processors on the current comm exist on the submesh comm" );
        subcommSize = comm.bcast( subcommSize, root_submesh );
        std::vector<int> subrank( comm.getSize() );
        comm.allGather( rank_submesh, &subrank[0] );
        std::vector<int> rank_map( subcommSize, -1 );
        for ( size_t i = 0; i < subrank.size(); i++ ) {
            if ( subrank[i] != -1 )
                rank_map[subrank[i]] = i;
        }
        // Get the rank of the proccessor that will own each meshElement
        for ( size_t i = 0; i < remote_ids.size(); i++ ) {
            if ( remote_ids[i].meshID() == meshID ) {
                int subowner_rank = remote_ids[i].owner_rank();
                AMP_ASSERT( rank_map[subowner_rank] != -1 );
                owner_rank[i] = rank_map[subowner_rank];
            }
        }
    }
    // Check that each element has a vaild owner rank
    int commSize = comm.getSize();
    for ( size_t i = 0; i < remote_ids.size(); i++ )
        AMP_ASSERT( owner_rank[i] >= 0 && owner_rank[i] < commSize );
    // Resort the remote ids according the the owner rank
    auto remote_ids2 = remote_ids;
    AMP::Utilities::quicksort( owner_rank, remote_ids2 );
    // Determine the send count and displacements for each processor
    std::vector<int> send_cnt( comm.getSize(), 0 );
    std::vector<int> send_disp( comm.getSize(), 0 );
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
    // Preform an allToAll to send the remote ids for DOF identification
    std::vector<int> recv_cnt( comm.getSize() );
    comm.allToAll<int>( 1, &send_cnt[0], &recv_cnt[0] );
    std::vector<int> recv_disp( comm.getSize() );
    size_t tot_size = recv_cnt[0];
    recv_disp[0]    = 0;
    for ( int i = 1; i < comm.getSize(); i++ ) {
        tot_size += recv_cnt[i];
        recv_disp[i] = recv_disp[i - 1] + recv_cnt[i - 1];
    }
    std::vector<AMP::Mesh::MeshElementID> recv_id( tot_size + 1 );
    AMP::Mesh::MeshElementID *send_buffer = nullptr;
    if ( !remote_ids2.empty() )
        send_buffer = &remote_ids2[0];
    N = comm.allToAll<AMP::Mesh::MeshElementID>(
        send_buffer, &send_cnt[0], &send_disp[0], &recv_id[0], &recv_cnt[0], &recv_disp[0], true );
    AMP_INSIST( N == (int) tot_size, "Unexpected receive size" );
    recv_id.resize( tot_size );
    // Determine the DOF for each received id
    std::vector<size_t> received_DOF( tot_size );
    for ( size_t i = 0; i < tot_size; i++ ) {
        int j = AMP::Utilities::findfirst( d_local_id, recv_id[i] );
        AMP_ASSERT( d_local_id[j] == recv_id[i] );
        received_DOF[i] = d_begin / d_DOFsPerElement + j;
    }
    // Send the DOFs back to the original processor
    std::vector<size_t> remote_dof;
    remote_dof.resize( remote_ids2.size() + 1, static_cast<size_t>( -1 ) );
    size_t *send_buffer_DOFs = nullptr;
    if ( tot_size > 0 )
        send_buffer_DOFs = &received_DOF[0];
    N = comm.allToAll<size_t>( send_buffer_DOFs,
                               &recv_cnt[0],
                               &recv_disp[0],
                               &remote_dof[0],
                               &send_cnt[0],
                               &send_disp[0],
                               true );
    AMP_INSIST( N == (int) remote_ids2.size(), "Unexpected receive size" );
    remote_dof.resize( remote_ids2.size() );
    // Sort the dofs back to the original order for the remote_ids
    AMP::Utilities::quicksort( remote_ids2, remote_dof );
    for ( size_t i = 0; i < remote_ids.size(); i++ )
        AMP_ASSERT( remote_ids[i] == remote_ids2[i] );
    return remote_dof;
}
} // namespace AMP::Discretization
