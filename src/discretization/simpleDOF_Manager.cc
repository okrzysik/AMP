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
    d_local_id.resize(mesh->numLocalElements(type));
    d_remote_id.resize(mesh->numGhostElements(type,gcw));
    AMP::Mesh::MeshIterator pos = mesh->getIterator(type,gcw);
    AMP::Mesh::MeshIterator end = pos.end();
    int i=0;
    int j=0;
    while ( pos != end ) {
        AMP::Mesh::MeshElementID id = pos->globalID();
        if ( id.is_local ) {
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
* Note:  this function is likely temporary, we are assuming     *
* all data will be stored on nodes.                             *
****************************************************************/
void simpleDOFManager::getDOFs( const AMP::Mesh::MeshElement &obj, std::vector <unsigned int> &ids, std::vector<unsigned int> which ) const
{
    std::vector<AMP::Mesh::MeshElement> elements = obj.getElements(d_type);
    if ( which.size()==0 ) {
        // Return all dofs
        ids.resize(elements.size());
        for (size_t i=0; i<elements.size(); i++) {
            AMP::Mesh::MeshElementID local_id = elements[i].globalID();
            ids[i] = AMP::Utilities::findfirst(d_local_id,local_id);
            AMP_INSIST(local_id==d_local_id[ids[i]],"Internal Error: id not found");
        }
    } else {
        // Return only the desired dof
        ids.resize(which.size());
        for (size_t i=0; i<which.size(); i++) {
            AMP::Mesh::MeshElementID local_id = elements[which[i]].globalID();
            ids[i] = AMP::Utilities::findfirst(d_local_id,local_id);
            AMP_INSIST(local_id==d_local_id[ids[i]],"Internal Error: id not found");
        }
    }
}


/****************************************************************
* Return the first D.O.F. on this core                          *
****************************************************************/
size_t simpleDOFManager::beginDOF( )
{
    return d_begin*DOFsPerElement;
}


/****************************************************************
* Return the last D.O.F. on this core                           *
****************************************************************/
size_t simpleDOFManager::endDOF( )
{
    return d_end*DOFsPerElement;
}


/****************************************************************
* Return the local number of D.O.F.s                           *
****************************************************************/
size_t simpleDOFManager::numLocalDOF( )
{
    return (d_end-d_begin)*DOFsPerElement;
}


/****************************************************************
* Return the global number of D.O.F.s                           *
****************************************************************/
size_t simpleDOFManager::numGlobalDOF( )
{
    return d_global*DOFsPerElement;
}


/****************************************************************
* Create a vector                                               *
****************************************************************/
AMP::LinearAlgebra::Vector::shared_ptr simpleDOFManager::createVector( AMP::LinearAlgebra::Variable::shared_ptr variable )
{
    // Check the inputs
    AMP::Mesh::GeomType type = (AMP::Mesh::GeomType) variable->variableID();
    if ( type != d_type )
        AMP_ERROR("The variableID must match the element type specified at construction");
    if ( (int) variable->DOFsPerObject() != DOFsPerElement )
        AMP_ERROR("The variableID must have the same number of DOFs per object as the DOF Manager");
    // Get the number or local DOFs and the start and end DOF indicies
    size_t N_global = d_mesh->numGlobalElements(type)*DOFsPerElement;
    size_t N_local = d_mesh->numLocalElements(type)*DOFsPerElement;
    size_t N_begin = d_begin*DOFsPerElement;
    size_t N_end = d_end*variable->DOFsPerObject();
    // Create the list of remote DOFs
    std::vector<unsigned int> remote_DOFs(d_remote_id.size()*DOFsPerElement,(unsigned int)-1);
    for (size_t i=0; i<d_remote_id.size(); i++) {
        for (int j=0; j<DOFsPerElement; j++)
            remote_DOFs[j+i*DOFsPerElement] = d_remote_dof[i]*DOFsPerElement + j;
    }
    // Create the communication list
    AMP::LinearAlgebra::CommunicationList::shared_ptr comm_list;
    if ( d_gcw == 0 ) {
        comm_list = AMP::LinearAlgebra::CommunicationList::createEmpty( N_local, d_mesh->getComm() );
    } else {
        AMP::LinearAlgebra::CommunicationListParameters::shared_ptr params( new AMP::LinearAlgebra::CommunicationListParameters );
        params->d_comm = d_mesh->getComm();
        params->d_localsize = N_local;
        params->d_remote_DOFs = remote_DOFs;
        comm_list = AMP::LinearAlgebra::CommunicationList::shared_ptr( new AMP::LinearAlgebra::CommunicationList(params) );
    }
    // Create the vector parameters
    boost::shared_ptr<AMP::LinearAlgebra::ManagedPetscVectorParameters> mvparams(
        new AMP::LinearAlgebra::ManagedPetscVectorParameters() );
    boost::shared_ptr<AMP::LinearAlgebra::EpetraVectorEngineParameters> eveparams(
        new AMP::LinearAlgebra::EpetraVectorEngineParameters( N_local, N_global, d_mesh->getComm() ) );
    int i = 0;
    for (size_t local_start=N_begin; local_start<N_end; local_start++, i++ ) {
        eveparams->addMapping ( i , local_start );
    }
    AMP::LinearAlgebra::VectorEngine::BufferPtr t_buffer ( new AMP::LinearAlgebra::VectorEngine::Buffer( N_local ) );
    AMP::LinearAlgebra::VectorEngine::shared_ptr epetra_engine( new AMP::LinearAlgebra::EpetraVectorEngine( eveparams, t_buffer ) );
    mvparams->d_Engine = epetra_engine;
    mvparams->d_CommList = comm_list;
    mvparams->d_DOFManager = shared_from_this();
    // Create the vector
    AMP::LinearAlgebra::Vector::shared_ptr vector = AMP::LinearAlgebra::Vector::shared_ptr( new AMP::LinearAlgebra::ManagedPetscVector(mvparams) );
    return vector;
}


/****************************************************************
* Create a matrix                                               *
****************************************************************/
AMP::LinearAlgebra::Matrix::shared_ptr simpleDOFManager::createMatrix( 
    AMP::LinearAlgebra::Variable::shared_ptr operandVar, 
    AMP::LinearAlgebra::Variable::shared_ptr resultVar )
{
    /*// Create the vectors (note: only square matricies are currently implimented)
    DOFManager *operandDOF = this;
    DOFManager *resultDOF = this;
    AMP::LinearAlgebra::Vector::shared_ptr  operandVec = operandDOF->createVector(operandVar);
    AMP::LinearAlgebra::Vector::shared_ptr  resultVec = resultDOF->createVector(resultVar);

    // Create the matrix parameters
    boost::shared_ptr<AMP::LinearAlgebra::ManagedPetscMatrixParameters> params( 
        new AMP::LinearAlgebra::ManagedPetscMatrixParameters( resultDOF->numLocalDOF()/resultVar->DOFsPerObject(),
                                                   resultDOF->numGlobalDOF()/resultVar->DOFsPerObject(),
                                                   0,
                                                   operandDOF->numGlobalDOF()/operandVar->DOFsPerObject(),
                                                   0,
                                                   d_mesh->getComm() ) );

    int multiplier = operandVar->DOFsPerObject();
    int divisor = resultVar->DOFsPerObject();
    NodalRowMap rowMap = pResultDofMap->getCommunicationList()->castTo<NodalRowMap>();
    size_t numLocalElements = pResultDofMap->numLocalElements();
    for ( unsigned int i = 0 ; i != numLocalElements; i++ )
    {
      params->addMapping ( i , pResultDofMap->beginDOF() + i );
      size_t nnz = rowMap.getNNZ ( i );
      params->setEntriesInRow ( i , nnz * multiplier/divisor );
      params->addColumns ( nnz * multiplier/divisor , (int *)rowMap.getColumns ( i*multiplier/divisor ) );
    }

    params->d_CommListLeft = d_vDOFMapCache[var_result->variableID()]->getCommunicationList();
    params->d_CommListRight = d_vDOFMapCache[var_operand->variableID()]->getCommunicationList();
    AMP::LinearAlgebra::Matrix::shared_ptr  newMatrix = AMP::LinearAlgebra::Matrix::shared_ptr ( new AMP::LinearAlgebra::ManagedPetscMatrix( params ) );
    size_t  mat_id = (var_result->variableID() << 10) + var_operand->variableID();
    d_vMatrixCache[mat_id] = newMatrix;

    double  values[1000];  // A little bit of a hack...
    for ( size_t i = 0 ; i != 1000 ; i++ )
        values[i] = 0.0;
    NodalRowMap rowMap2 = pOperandDofMap->getCommunicationList()->castTo<NodalRowMap>();
    for ( size_t i=0 ; i!=numLocalElements; i++ )
    {
      int cur_row_id = (int)i + pResultDofMap->beginDOF();
      int related_col = i * multiplier/divisor;
      newMatrix->castTo<AMP::LinearAlgebra::ManagedMatrix>().createValuesByGlobalID ( 1 ,
                                          rowMap2.getNNZ (related_col) ,
                                         &cur_row_id ,
                                   (int *)rowMap2.getColumns(related_col) ,
                                          values );
    }
    newMatrix->castTo<AMP::LinearAlgebra::EpetraMatrix>().setEpetraMaps ( pResultVec , pOperandVec );
    newMatrix->makeConsistent ();
*/
    AMP_ERROR("Not implimented yet");
    return AMP::LinearAlgebra::Matrix::shared_ptr();

}


/****************************************************************
* Find the remote DOF given a set of mesh element IDs           *
* Note: for this function to work correctly, the remote ids     *
* must be sorted, and d_local_id must be set                    *
****************************************************************/
std::vector<size_t> simpleDOFManager::getRemoteDOF(std::vector<AMP::Mesh::MeshElementID> remote_ids )
{
    // Get the rank that will own each MeshElement on the current communicator
    std::set<size_t> meshIDs;
    for (size_t i=0; i<remote_ids.size(); i++)
        meshIDs.insert(remote_ids[i].meshID);
    AMP_MPI comm = d_mesh->getComm();
    std::vector<int> owner_rank(remote_ids.size(),-1);
    for (std::set<size_t>::iterator it=meshIDs.begin() ; it!=meshIDs.end(); it++) {
        // Get the mesh with the given meshID
        size_t meshID = *it;
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
            if ( remote_ids[i].meshID == meshID ) {
                int subowner_rank = remote_ids[i].owner_rank;
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
    int N = comm.allToAll<AMP::Mesh::MeshElementID>( 
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
    std::vector<size_t> remote_dof(d_remote_id.size()+1,static_cast<size_t>(-1));
    size_t* send_buffer_DOFs = NULL;
    if ( tot_size > 0 )
        send_buffer_DOFs = &recieved_DOF[0];
    N = comm.allToAll<size_t>( 
        send_buffer_DOFs, &recv_cnt[0], &recv_disp[0], 
        &remote_dof[0], &send_cnt[0], &send_disp[0], true);
    AMP_INSIST(N==(int)d_remote_id.size(),"Unexpected recieve size");
    remote_dof.resize(d_remote_id.size());
    // Sort the dofs back to the original order for the remote_ids
    AMP::Utilities::quicksort(remote_ids2,remote_dof);
    for (size_t i=0; i<remote_ids.size(); i++)
        AMP_ASSERT(remote_ids[i]==remote_ids2[i]);
    return remote_dof;
}


}
}

