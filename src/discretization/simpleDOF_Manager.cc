#include "discretization/simpleDOF_Manager.h"
#include "utils/Utilities.h"

#include "utils/Utilities.h"

namespace AMP {
namespace Discretization {


/****************************************************************
* Constructors                                                  *
****************************************************************/
simpleDOFManager::simpleDOFManager( boost::shared_ptr<AMP::Mesh::Mesh> mesh, AMP::Mesh::GeomType type, int gcw, int DOFsPerObject )
{
    d_mesh = mesh;
    d_type = type;
    d_gcw = gcw;
    DOFsPerElement = DOFsPerObject;
    // Create a sorted list of the local and remote types
    d_local_id.resize(mesh->numLocalElements(d_type));
    d_remote_id.resize(mesh->numGhostElements(d_type,d_gcw));
    AMP::Mesh::MeshIterator pos = d_mesh->getIterator(d_type,d_gcw);
    AMP::Mesh::MeshIterator end = pos.end();
    int i=0;
    int j=0;
    while ( pos != end ) {
        AMP::Mesh::MeshElementID id = pos->globalID();
        if ( id.is_local() ) {
            d_local_id[i] = id;
            i++;
        } else {
            d_remote_id[j] = id;
            j++;
        }
        ++pos;
    }
    // Sort the elements (they will be sorted by the meshID, then the rank on the 
    // comm of the given mesh, then the element type, and finally the local id)
    AMP::Utilities::quicksort(d_local_id);
    AMP::Utilities::quicksort(d_remote_id);
    // Get the number of local elements per processor and the global number of DOFs
    AMP_MPI comm = d_mesh->getComm();
    size_t N_local = d_local_id.size();
    comm.sumScan<size_t>(&N_local,&d_end);
    d_begin = d_end-N_local;
    d_global = comm.bcast(d_end,comm.getSize()-1);
    // Determine the remote DOFs (assuming 1 DOF per node)
    // Note: this must be done after d_local_id is set, d_begin and d_global are set, and remote_ids must be sorted.
    d_remote_dof = getRemoteDOF(d_remote_id);
}


/****************************************************************
* Get the entry indices of nodal values given a mesh element    *
****************************************************************/
void simpleDOFManager::getDOFs( const AMP::Mesh::MeshElement &obj, std::vector <size_t> &dofs, std::vector<size_t> which ) const
{
    std::vector<AMP::Mesh::MeshElement> elements;
    if ( obj.elementType() == d_type )
        elements = std::vector<AMP::Mesh::MeshElement>(1,obj);
    else
        elements = obj.getElements(d_type);
    if ( which.size()==0 ) {
        // Return all dofs
        dofs.resize(elements.size()*DOFsPerElement);
        for (size_t i=0; i<elements.size(); i++) {
            AMP::Mesh::MeshElementID local_id = elements[i].globalID();
            size_t index = AMP::Utilities::findfirst(d_local_id,local_id);
            AMP_INSIST(local_id==d_local_id[index],"Internal Error: id not found");
            for (int j=0; j<DOFsPerElement; j++)
                dofs[i*DOFsPerElement+j] = (index+d_begin)*DOFsPerElement + j;
        }
    } else {
        // Return only the desired dof
        dofs.resize(which.size()*DOFsPerElement);
        for (size_t i=0; i<which.size(); i++) {
            AMP::Mesh::MeshElementID local_id = elements[which[i]].globalID();
            size_t index = AMP::Utilities::findfirst(d_local_id,local_id);
            AMP_INSIST(local_id==d_local_id[index],"Internal Error: id not found");
            for (int j=0; j<DOFsPerElement; j++)
                dofs[i*DOFsPerElement+j] = (index+d_begin)*DOFsPerElement + j;
        }
    }
}
void simpleDOFManager::getDOFs( const AMP::Mesh::MeshElementID &id, std::vector <size_t> &dofs ) const
{
    dofs.resize(DOFsPerElement);
    size_t index = AMP::Utilities::findfirst(d_local_id,id);
    AMP_INSIST(id==d_local_id[index],"Internal Error: id not found");
    for (int j=0; j<DOFsPerElement; j++)
        dofs[j] = (index+d_begin)*DOFsPerElement + j;
}


/****************************************************************
* Get an entry over the mesh elements associated with the DOFs  *
****************************************************************/
AMP::Mesh::MeshIterator simpleDOFManager::getIterator( ) const
{
    return d_mesh->getIterator(d_type,0);
}


/****************************************************************
* Return the first D.O.F. on this core                          *
****************************************************************/
size_t simpleDOFManager::beginDOF( ) const
{
    return d_begin*DOFsPerElement;
}


/****************************************************************
* Return the last D.O.F. on this core                           *
****************************************************************/
size_t simpleDOFManager::endDOF( ) const
{
    return d_end*DOFsPerElement;
}


/****************************************************************
* Return the local number of D.O.F.s                           *
****************************************************************/
size_t simpleDOFManager::numLocalDOF( ) const
{
    return (d_end-d_begin)*DOFsPerElement;
}


/****************************************************************
* Return the global number of D.O.F.s                           *
****************************************************************/
size_t simpleDOFManager::numGlobalDOF( ) const
{
    return d_global*DOFsPerElement;
}


/****************************************************************
* Return the communicator                                       *
****************************************************************/
AMP_MPI simpleDOFManager::getComm( ) const
{
    return d_mesh->getComm();
}


/****************************************************************
* Return the remote DOFs for a vector                           *
****************************************************************/
std::vector<size_t> simpleDOFManager::getRemoteDOFs( ) const
{
    // Create the list of remote DOFs
    std::vector<size_t> remote_DOFs(d_remote_id.size()*DOFsPerElement,(size_t)-1);
    for (size_t i=0; i<d_remote_id.size(); i++) {
        for (int j=0; j<DOFsPerElement; j++)
            remote_DOFs[j+i*DOFsPerElement] = d_remote_dof[i]*DOFsPerElement + j;
    }
    AMP::Utilities::quicksort(remote_DOFs);
    return remote_DOFs;
}


/****************************************************************
* Return the global number of D.O.F.s                           *
****************************************************************/
std::vector<size_t> simpleDOFManager::getRowDOFs( const AMP::Mesh::MeshElement &obj ) const
{
    AMP_INSIST(obj.elementType()==d_type,"Mixing types is not tested/supported yet");
    std::vector< Mesh::MeshElement::shared_ptr > neighbor_elements = obj.getNeighbors();
    std::vector< const Mesh::MeshElement* > elements(neighbor_elements.size()+1,&obj);
    for (size_t i=0; i<neighbor_elements.size(); i++)
        elements[i+1] = neighbor_elements[i].get();
    std::vector<size_t> ids(elements.size()*DOFsPerElement);
    size_t *ids2 = &ids[0];     // Use the pointer directly for speed
    for (size_t i=0; i<elements.size(); i++) {
        AMP::Mesh::MeshElementID  id = elements[i]->globalID();
        if ( id.is_local() ) {
            // We are dealing with a local element
            size_t index = AMP::Utilities::findfirst(d_local_id,id);
            AMP_INSIST(id==d_local_id[index],"Internal Error: id not found");
            for (int j=0; j<DOFsPerElement; j++)
                ids2[i*DOFsPerElement+j] = (index+d_begin)*DOFsPerElement + j;
        } else {
            // We are dealing with a remote element, hopefully we know where it is
            size_t index = AMP::Utilities::findfirst(d_remote_id,id);
            AMP_INSIST(id==d_remote_id[index],"Internal Error: remote id not found");
            for (int j=0; j<DOFsPerElement; j++)
                ids2[i*DOFsPerElement+j] = d_remote_dof[index]*DOFsPerElement + j;
        }
    }
    AMP::Utilities::quicksort(ids);
    return ids;
}


/****************************************************************
* Find the remote DOF given a set of mesh element IDs           *
* Note: for this function to work correctly, the remote ids     *
* must be sorted, and d_local_id must be set                    *
****************************************************************/
std::vector<size_t> simpleDOFManager::getRemoteDOF(std::vector<AMP::Mesh::MeshElementID> remote_ids ) const
{
    AMP_MPI comm = d_mesh->getComm();
    // Get the set of mesh ids (must match on all processors)
    std::set<AMP::Mesh::MeshID> meshIDs;
    for (size_t i=0; i<remote_ids.size(); i++)
        meshIDs.insert(remote_ids[i].meshID());
    std::vector<AMP::Mesh::MeshID> tmpLocalIDs(meshIDs.begin(),meshIDs.end());
    int N = (int) comm.sumReduce<size_t>(tmpLocalIDs.size());
    if ( N==0 ) {
        // Nobody has any remote ids to identify
        return std::vector<size_t>();
    }
    AMP::Mesh::MeshID *send_ptr=NULL;
    if ( tmpLocalIDs.size()>0 )
        send_ptr = &tmpLocalIDs[0];
    std::vector<AMP::Mesh::MeshID> tmpGlobalIDs(N);
    int N_recv = comm.allGather<AMP::Mesh::MeshID>(send_ptr,tmpLocalIDs.size(),&tmpGlobalIDs[0]);
    AMP_ASSERT(N_recv==N);
    for (size_t i=0; i<tmpGlobalIDs.size(); i++)
        meshIDs.insert(tmpGlobalIDs[i]);
    // Get the rank that will own each MeshElement on the current communicator
    std::vector<int> owner_rank(remote_ids.size(),-1);
    for (std::set<AMP::Mesh::MeshID>::iterator it=meshIDs.begin() ; it!=meshIDs.end(); it++) {
        // Get the mesh with the given meshID
        AMP::Mesh::MeshID meshID = *it;
        AMP::Mesh::Mesh::shared_ptr submesh = d_mesh->Subset(meshID);
        // Create a map from the rank of the submesh to the current mesh
        int rank_submesh = -1;
        int root_submesh = comm.getSize();
        int subcommSize = -1;
        if ( submesh.get() != NULL ) {
            AMP_MPI subcomm = submesh->getComm();
            rank_submesh = subcomm.getRank();
            root_submesh = comm.getRank();
            subcommSize = subcomm.getSize();
        }
        root_submesh = comm.minReduce(root_submesh);
        if ( root_submesh==comm.getSize() )
            AMP_ERROR("Not processors on the current comm exist on the submesh comm");
        subcommSize = comm.bcast(subcommSize,root_submesh);
        std::vector<int> subrank(comm.getSize());
        comm.allGather(rank_submesh,&subrank[0]);
        std::vector<int> rank_map(subcommSize,-1);
        for (size_t i=0; i<subrank.size(); i++) {
            if ( subrank[i] != -1 )
                rank_map[subrank[i]] = i;
        }
        // Get the rank of the proccessor that will own each meshElement
        for (size_t i=0; i<remote_ids.size(); i++) {
            if ( remote_ids[i].meshID() == meshID ) {
                int subowner_rank = remote_ids[i].owner_rank();
                AMP_ASSERT(rank_map[subowner_rank]!=-1);
                owner_rank[i] = rank_map[subowner_rank];
            }
        }
    }
    // Check that each element has a vaild owner rank
    int commSize = comm.getSize();
    for (size_t i=0; i<remote_ids.size(); i++)
        AMP_ASSERT(owner_rank[i]>=0&&owner_rank[i]<commSize);
    // Resort the remote ids according the the owner rank
    std::vector<AMP::Mesh::MeshElementID> remote_ids2 = remote_ids;
    AMP::Utilities::quicksort(owner_rank,remote_ids2);
    // Determine the send count and displacements for each processor
    std::vector<int> send_cnt(comm.getSize(),0);
    std::vector<int> send_disp(comm.getSize(),0);
    int rank = 0;
    size_t start = 0;
    size_t index = 0;
    while ( index < owner_rank.size() ) {
        if ( owner_rank[index] < rank ) {
            AMP_ERROR("This should not occur");
        } else if ( owner_rank[index] == rank ) {
            // Move to the next element
            index++;
        } else {
            // Store the number of elements with the given rank, and move to the next rank
            send_disp[rank] = start;
            send_cnt[rank] = index-start;
            start = index;
            rank++;
        }
    }
    send_disp[rank] = start;
    send_cnt[rank] = index-start;
    // Preform an allToAll to send the remote ids for DOF identification
    std::vector<int> recv_cnt(comm.getSize());
    comm.allToAll<int>( 1, &send_cnt[0], &recv_cnt[0] );
    std::vector<int> recv_disp(comm.getSize());
    size_t tot_size = recv_cnt[0];
    recv_disp[0] = 0;
    for (int i=1; i<comm.getSize(); i++) {
        tot_size += recv_cnt[i];
        recv_disp[i] = recv_disp[i-1] + recv_cnt[i-1];
    }
    std::vector<AMP::Mesh::MeshElementID> recv_id(tot_size+1); 
    AMP::Mesh::MeshElementID* send_buffer = NULL;
    if ( remote_ids2.size() > 0 )
        send_buffer = &remote_ids2[0];
    N = comm.allToAll<AMP::Mesh::MeshElementID>( 
        send_buffer, &send_cnt[0], &send_disp[0], 
        &recv_id[0], &recv_cnt[0], &recv_disp[0], true);
    AMP_INSIST(N==(int)tot_size,"Unexpected recieve size");
    recv_id.resize(tot_size);
    // Determine the DOF for each recieved id
    std::vector<size_t> recieved_DOF(tot_size);
    for (size_t i=0; i<tot_size; i++) {
        int j = AMP::Utilities::findfirst(d_local_id,recv_id[i]);
        AMP_ASSERT(d_local_id[j]==recv_id[i]);
        recieved_DOF[i] = d_begin + j;
    }
    // Send the DOFs back to the original processor
    std::vector<size_t> remote_dof;
    remote_dof.resize(remote_ids2.size()+1,static_cast<size_t>(-1));
    size_t* send_buffer_DOFs = NULL;
    if ( tot_size > 0 )
        send_buffer_DOFs = &recieved_DOF[0];
    N = comm.allToAll<size_t>( 
        send_buffer_DOFs, &recv_cnt[0], &recv_disp[0], 
        &remote_dof[0], &send_cnt[0], &send_disp[0], true);
    AMP_INSIST(N==(int)remote_ids2.size(),"Unexpected recieve size");
    remote_dof.resize(remote_ids2.size());
    // Sort the dofs back to the original order for the remote_ids
    AMP::Utilities::quicksort(remote_ids2,remote_dof);
    for (size_t i=0; i<remote_ids.size(); i++)
        AMP_ASSERT(remote_ids[i]==remote_ids2[i]);
    return remote_dof;
}


}
}

