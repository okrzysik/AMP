#include "discretization/simpleDOF_Manager.h"
#include "utils/Utilities.h"

#include "vectors/petsc/ManagedPetscVector.h"
#include "vectors/trilinos/EpetraVectorEngine.h"
#include "matrices/petsc/ManagedPetscMatrix.h"
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
void simpleDOFManager::getDOFs( const AMP::Mesh::MeshElement &obj, std::vector <unsigned int> &dofs, std::vector<unsigned int> which ) const
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
void simpleDOFManager::getDOFs( const AMP::Mesh::MeshElementID &id, std::vector <unsigned int> &dofs ) const
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
* Create a vector                                               *
****************************************************************/
AMP::LinearAlgebra::Vector::shared_ptr simpleDOFManager::createVector( AMP::LinearAlgebra::Variable::shared_ptr variable )
{
    // Check the inputs
    AMP::Mesh::GeomType type = (AMP::Mesh::GeomType) variable->variableID();
    boost::shared_ptr<AMP::Discretization::DOFManager> DOFs = shared_from_this();
    if ( type != d_type )
        AMP_ERROR("The variableID must match the element type specified at construction");
    if ( (int) variable->DOFsPerObject() != DOFsPerElement )
        AMP_ERROR("The variableID must have the same number of DOFs per object as the DOF Manager");
    // Create the communication list
    AMP::LinearAlgebra::CommunicationList::shared_ptr comm_list;
    if ( d_gcw == 0 ) {
        comm_list = AMP::LinearAlgebra::CommunicationList::createEmpty( DOFs->numLocalDOF(), DOFs->getComm() );
    } else {
        AMP::LinearAlgebra::CommunicationListParameters::shared_ptr params( new AMP::LinearAlgebra::CommunicationListParameters );
        params->d_comm = DOFs->getComm();
        params->d_localsize = DOFs->numLocalDOF();
        params->d_remote_DOFs = DOFs->getRemoteDOFs();
        comm_list = AMP::LinearAlgebra::CommunicationList::shared_ptr( new AMP::LinearAlgebra::CommunicationList(params) );
    }
    // Create the vector parameters
    boost::shared_ptr<AMP::LinearAlgebra::ManagedPetscVectorParameters> mvparams(
        new AMP::LinearAlgebra::ManagedPetscVectorParameters() );
    boost::shared_ptr<AMP::LinearAlgebra::EpetraVectorEngineParameters> eveparams(
        new AMP::LinearAlgebra::EpetraVectorEngineParameters( DOFs->numLocalDOF(), DOFs->numGlobalDOF(), DOFs->getComm() ) );
    int i = 0;
    for (size_t local_start=DOFs->beginDOF(); local_start<DOFs->endDOF(); local_start++, i++ ) {
        eveparams->addMapping( i, local_start );
    }
    AMP::LinearAlgebra::VectorEngine::BufferPtr t_buffer ( new AMP::LinearAlgebra::VectorEngine::Buffer( DOFs->numLocalDOF() ) );
    AMP::LinearAlgebra::VectorEngine::shared_ptr epetra_engine( new AMP::LinearAlgebra::EpetraVectorEngine( eveparams, t_buffer ) );
    mvparams->d_Engine = epetra_engine;
    mvparams->d_CommList = comm_list;
    mvparams->d_DOFManager = DOFs;
    // Create the vector
    AMP::LinearAlgebra::Vector::shared_ptr vector = AMP::LinearAlgebra::Vector::shared_ptr( new AMP::LinearAlgebra::ManagedPetscVector(mvparams) );
    vector->setVariable(variable);
    return vector;
}


/****************************************************************
* Create a matrix                                               *
****************************************************************/
AMP::LinearAlgebra::Matrix::shared_ptr simpleDOFManager::createMatrix( 
    AMP::LinearAlgebra::Vector::shared_ptr operandVec, 
    AMP::LinearAlgebra::Vector::shared_ptr resultVec )
{

    // Get the DOFs
    AMP::Discretization::DOFManager::shared_ptr operandDOF = operandVec->getDOFManager();
    AMP::Discretization::DOFManager::shared_ptr resultDOF = resultVec->getDOFManager();
    if ( operandDOF->getComm().compare(resultVec->getComm()) == 0 )
        AMP_ERROR("operandDOF and resultDOF on different comm groups is NOT tested, and needs to be fixed");

    // Create the matrix parameters
    boost::shared_ptr<AMP::LinearAlgebra::ManagedPetscMatrixParameters> params( 
        new AMP::LinearAlgebra::ManagedPetscMatrixParameters ( resultDOF->numLocalDOF(),
        resultDOF->numGlobalDOF(), 0, operandDOF->numGlobalDOF(), 0, d_mesh->getComm() ) );

    // Add the rows to the matrix parameters
    AMP::Mesh::MeshIterator cur_elem = resultDOF->getIterator();
    AMP::Mesh::MeshIterator end_elem = cur_elem.end();
    int columns[1000];   // Preallocate for the columns for speed
    while ( cur_elem != end_elem) {
        AMP::Mesh::MeshElement obj = *cur_elem;
        // Get the result DOFs associated with the given element
        std::vector<unsigned int> ids;
        resultDOF->getDOFs(obj,ids);
        // Get the operand DOFs associated with the given element
        std::vector<size_t> row = operandDOF->getRowDOFs(obj);
        size_t nnz = row.size();
        for (size_t i=0; i<row.size(); i++)
            columns[i] = (int) row[i];
        // Add the rows
        for (size_t i=0; i<ids.size(); i++) {
            int globalRowID = ids[i];
            int localRowID = globalRowID - resultDOF->beginDOF();
            params->addMapping( localRowID, globalRowID );
            params->setEntriesInRow( localRowID, nnz );
        }
        // Add the columns
        params->addColumns( nnz, columns );
        // Increment the iterator (pre-increment for speed)
        ++cur_elem;
    }

    // Get the communication lists for the vectors
    params->d_CommListLeft = resultVec->getCommunicationList();
    params->d_CommListRight = operandVec->getCommunicationList();

    // Create the matrix
    boost::shared_ptr<AMP::LinearAlgebra::ManagedPetscMatrix>  newMatrix( new AMP::LinearAlgebra::ManagedPetscMatrix(params) );

    // Initialize the matrix
    cur_elem = resultDOF->getIterator();
    end_elem = cur_elem.end();
    double  values[1000];
    for (size_t i=0 ; i<1000; i++) { values[i] = 0.0; }
    while ( cur_elem != end_elem) {
        AMP::Mesh::MeshElement obj = *cur_elem;
        // Get the result DOFs associated with the given element
        std::vector<unsigned int> ids;
        resultDOF->getDOFs(obj,ids);
        // Get the operand DOFs associated with the given element
        std::vector<size_t> row = operandDOF->getRowDOFs(obj);
        size_t nnz = row.size();
        for (size_t i=0; i<row.size(); i++)
            columns[i] = (int) row[i];
        // Add the rows
        for (size_t i=0; i<ids.size(); i++) {
            int globalRowID = ids[i];
            newMatrix->createValuesByGlobalID( 1, nnz, &globalRowID, columns, values );
        }
        ++cur_elem;
    }
    newMatrix->castTo<AMP::LinearAlgebra::EpetraMatrix>().setEpetraMaps( resultVec, operandVec );
    newMatrix->makeConsistent();

    return newMatrix;
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

