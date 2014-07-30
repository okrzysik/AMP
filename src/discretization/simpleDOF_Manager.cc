#include "discretization/simpleDOF_Manager.h"
#include "discretization/MultiDOF_Manager.h"
#include "ampmesh/MultiMesh.h"
#include "utils/Utilities.h"
#include <vector>
#include <set>


namespace AMP {
namespace Discretization {


/****************************************************************
* Constructors                                                  *
****************************************************************/
DOFManager::shared_ptr  simpleDOFManager::create( boost::shared_ptr<AMP::Mesh::Mesh> mesh, 
    AMP::Mesh::GeomType type, int gcw, int DOFsPerObject, bool split )
{
    if ( mesh.get()==NULL )
        return DOFManager::shared_ptr();
    if ( split && boost::dynamic_pointer_cast<AMP::Mesh::MultiMesh>(mesh).get()!=NULL ) {
        // We want to split the DOFs by the mesh
        std::vector<AMP::Mesh::MeshID> meshIDs = mesh->getLocalBaseMeshIDs();
        std::vector<DOFManager::shared_ptr> managers;
        for (size_t i=0; i<meshIDs.size(); i++) {
            AMP::Mesh::Mesh::shared_ptr subMesh = mesh->Subset(meshIDs[i]);
            if ( subMesh.get() != NULL )
                managers.push_back( create( subMesh, type, gcw, DOFsPerObject, false ) );
        }
        boost::shared_ptr<multiDOFManager> rtn( new multiDOFManager( mesh->getComm(), managers ) );
        return rtn;
    } 
    // We are ready to create the simpleDOFManager
    boost::shared_ptr<simpleDOFManager> rtn( new simpleDOFManager() );
    rtn->d_mesh = mesh;
    rtn->d_type = type;
    rtn->d_comm = mesh->getComm();
    rtn->DOFsPerElement = DOFsPerObject;
    rtn->d_localIterator = mesh->getIterator(type,0);
    rtn->d_ghostIterator = mesh->getIterator(type,gcw);
    rtn->initialize();
    return rtn;
}
DOFManager::shared_ptr  simpleDOFManager::create( boost::shared_ptr<AMP::Mesh::Mesh> mesh, const AMP::Mesh::MeshIterator it1, const AMP::Mesh::MeshIterator it2, int DOFsPerElement )
{
    // Check the iterators
    AMP::Mesh::MeshIterator intersection = AMP::Mesh::Mesh::getIterator( AMP::Mesh::Intersection, it1, it2 );
    AMP_INSIST(intersection.size()==it2.size(),"it1 must include it2");
    AMP::Mesh::MeshIterator tmp = it2.begin();
    for (size_t i=0; i<tmp.size(); i++) {
        AMP::Mesh::MeshElementID id = tmp->globalID();
        AMP_INSIST(id.is_local(),"it2 may not contain any ghost elements");
        ++tmp;
    }
    tmp = it1.begin();
    AMP::Mesh::GeomType type = tmp->globalID().type();
    for (size_t i=0; i<tmp.size(); i++) {
        AMP::Mesh::MeshElementID id = tmp->globalID();
        AMP_INSIST(id.type()==type,"All elements in the iterator must be the same type");
        ++tmp;
    }
    // Create the simpleDOFManager
    boost::shared_ptr<simpleDOFManager> rtn( new simpleDOFManager() );
    rtn->d_mesh = mesh;
    rtn->d_type = type;
    rtn->d_comm = mesh->getComm();
    rtn->DOFsPerElement = DOFsPerElement;
    rtn->d_ghostIterator = it1;
    rtn->d_localIterator = it2;
    rtn->initialize();
    return rtn;
}
DOFManager::shared_ptr  simpleDOFManager::create( const AMP::Mesh::MeshIterator it, int DOFsPerElement )
{
    // Check the iterator
    AMP::Mesh::MeshIterator tmp = it.begin();
    AMP::Mesh::GeomType type = tmp->globalID().type();
    for (size_t i=0; i<tmp.size(); i++) {
        AMP::Mesh::MeshElementID id = tmp->globalID();
        AMP_INSIST(id.type()==type,"All elements in the iterator must be the same type");
        ++tmp;
    }
    // Create the simpleDOFManager
    boost::shared_ptr<simpleDOFManager> rtn( new simpleDOFManager() );
    rtn->d_mesh = AMP::Mesh::Mesh::shared_ptr();
    rtn->d_type = type;
    rtn->d_comm = AMP_MPI(AMP_COMM_SELF);
    rtn->DOFsPerElement = DOFsPerElement;
    rtn->d_ghostIterator = it;
    rtn->d_localIterator = it;
    rtn->initialize();
    return rtn;
}


/****************************************************************
* Deconstructor                                                 *
****************************************************************/
simpleDOFManager::~simpleDOFManager( )
{
}


/****************************************************************
* Initialize the data                                           *
****************************************************************/
void simpleDOFManager::initialize()
{
    // Create a sorted list of the local and remote types
    d_local_id.resize(d_localIterator.size());
    d_remote_id.resize(d_ghostIterator.size()-d_localIterator.size());
    AMP::Mesh::MeshIterator pos = d_ghostIterator.begin();
    AMP::Mesh::MeshIterator end = d_ghostIterator.end();
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
    size_t N_local = d_local_id.size()*DOFsPerElement;
    d_comm.sumScan<size_t>(&N_local,&d_end);
    d_begin = d_end-N_local;
    d_global = d_comm.bcast(d_end,d_comm.getSize()-1);
    // Determine the remote DOFs (assuming 1 DOF per node)
    // Note: this must be done after d_local_id is set, d_begin and d_global are set, and remote_ids must be sorted.
    d_remote_dof = getRemoteDOF(d_remote_id);
}


/****************************************************************
* Subset the DOF manager                                        *
****************************************************************/
boost::shared_ptr<DOFManager>  simpleDOFManager::subset( const AMP::Mesh::Mesh::shared_ptr mesh, bool useMeshComm )
{

    // Check if we are dealing with a single mesh for both the internal and desired mesh
    if ( mesh->isBaseMesh() && d_mesh->isBaseMesh() ) {
        if ( mesh->meshID() == d_mesh->meshID() )
            return shared_from_this();
        else
            return boost::shared_ptr<DOFManager>();
    } 
    // Check if the desired mesh is a multimesh that contains the current mesh
    std::vector<AMP::Mesh::MeshID> ids = mesh->getLocalMeshIDs();
    bool found_local = false;
    for (size_t i=0; i<ids.size(); i++) {
        if ( ids[i] == d_mesh->meshID() )
            found_local = true;
    }
    AMP_MPI comm(AMP_COMM_NULL);
    if ( useMeshComm ) {
        AMP_MPI comm = AMP_MPI::intersect( d_comm, mesh->getComm() );
    } else {
        comm = d_comm;
    }
    found_local = comm.allReduce( found_local );
    if ( found_local ) 
        return shared_from_this();
    // We were not able to use an efficient subset, use the generic base function
    return DOFManager::subset( mesh, useMeshComm );
}


/****************************************************************
* Get the DOFs for the element                                  *
****************************************************************/
void simpleDOFManager::getDOFs( const AMP::Mesh::MeshElementID &id, std::vector <size_t> &dofs ) const
{
    dofs.resize(0);
    // Search for the dof locally
    size_t index = AMP::Utilities::findfirst(d_local_id,id);
    if ( index == d_local_id.size() ) { index--; }
    if ( id==d_local_id[index] ) {
        // The id was found
        dofs.resize(DOFsPerElement);
        for (int j=0; j<DOFsPerElement; j++)
            dofs[j] = index*DOFsPerElement + d_begin + j;
        return;
    } 
    // Search for the dof in the remote list
    if ( !d_remote_id.empty() && dofs.empty() ) {
        index = AMP::Utilities::findfirst(d_remote_id,id);
        if ( index == d_remote_id.size() ) { index--; }
        if ( id==d_remote_id[index] ) {
            // The id was found
            dofs.resize(DOFsPerElement);
            for (int j=0; j<DOFsPerElement; j++)
                dofs[j] = d_remote_dof[index]*DOFsPerElement + j;
            return;
        } 
    }
}


/****************************************************************
* Get an entry over the mesh elements associated with the DOFs  *
****************************************************************/
AMP::Mesh::MeshIterator simpleDOFManager::getIterator( ) const
{
    return d_localIterator.begin();
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
* Return the row DOFs                                           *
****************************************************************/
std::vector<size_t> simpleDOFManager::getRowDOFs( const AMP::Mesh::MeshElement &obj ) const
{
    // Get a list of all element ids that are part of the row
    AMP::Mesh::GeomType meshType = d_mesh->getGeomType();
    AMP::Mesh::GeomType objType = obj.elementType();
    std::vector<AMP::Mesh::MeshElementID> ids;
    if ( objType==d_type && ( objType==AMP::Mesh::Vertex || objType==meshType ) ) {
        // Use the getNeighbors function to get the neighbors of the current element
        std::vector<Mesh::MeshElement::shared_ptr> neighbor_elements = obj.getNeighbors();
        ids.reserve(neighbor_elements.size()+1);
        ids.push_back(obj.globalID());
        for (size_t i=0; i<neighbor_elements.size(); i++) {
            if ( neighbor_elements[i].get() != NULL )
                ids.push_back(neighbor_elements[i]->globalID());
        }
    } else if ( objType==d_type ) {
        // We need to use the mesh to get the connectivity of the elements of the same type
        std::vector<AMP::Mesh::MeshElement> parents = d_mesh->getElementParents(obj,meshType);
        for (size_t i=0; i<parents.size(); i++) {
            std::vector<AMP::Mesh::MeshElement> children = parents[i].getElements(objType);
            ids.reserve(ids.size()+children.size());
            for (size_t j=0; j<children.size(); j++)
                ids.push_back(children[j].globalID());
        }
        AMP::Utilities::unique(ids);
    } else if ( objType>d_type ) {
        // The desired element type is < the current element type, use getElements
        std::vector<AMP::Mesh::MeshElement> children = obj.getElements(d_type);
        for (size_t i=0; i<children.size(); i++)
            ids.push_back(children[i].globalID());
    } else if ( objType<d_type ) {
        // The desired element type is < the current element type, use getElementParents
        std::vector<AMP::Mesh::MeshElement> parents = d_mesh->getElementParents(obj,meshType);
        for (size_t i=0; i<parents.size(); i++)
            ids.push_back(parents[i].globalID());
    } else {
        AMP_ERROR("Internal error");
    }
    // Get all dofs for each element id
    std::vector<size_t> dofs;
    dofs.reserve(ids.size()*DOFsPerElement);
    std::vector<size_t> dofs2(DOFsPerElement);
    for (size_t i=0; i<ids.size(); i++) {
        getDOFs( ids[i], dofs2 );
        for (size_t j=0; j<dofs2.size(); j++)
            dofs.push_back(dofs2[j]);
    }
    // Sort the row dofs
    AMP::Utilities::quicksort(dofs);
    return dofs;
}


/****************************************************************
* Find the remote DOF given a set of mesh element IDs           *
* Note: for this function to work correctly, the remote ids     *
* must be sorted, and d_local_id must be set                    *
****************************************************************/
std::vector<size_t> simpleDOFManager::getRemoteDOF(std::vector<AMP::Mesh::MeshElementID> remote_ids ) const
{
    if ( d_comm.getSize()==1 )
        return std::vector<size_t>();     // There are no remote DOFs
    // Get the set of mesh ids (must match on all processors)
    AMP_MPI comm = d_mesh->getComm();
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
    if ( !tmpLocalIDs.empty() )
        send_ptr = &tmpLocalIDs[0];
    std::vector<AMP::Mesh::MeshID> tmpGlobalIDs(N);
    int N_recv = comm.allGather<AMP::Mesh::MeshID>(send_ptr,tmpLocalIDs.size(),&tmpGlobalIDs[0]);
    AMP_ASSERT(N_recv==N);
    for (size_t i=0; i<tmpGlobalIDs.size(); i++)
        meshIDs.insert(tmpGlobalIDs[i]);
    // Get the rank that will own each MeshElement on the current communicator
    std::vector<int> owner_rank(remote_ids.size(),-1);
    for (std::set<AMP::Mesh::MeshID>::iterator it=meshIDs.begin() ; it!=meshIDs.end(); ++it) {
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
    if ( !remote_ids2.empty() )
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
        recieved_DOF[i] = d_begin/DOFsPerElement + j;
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

