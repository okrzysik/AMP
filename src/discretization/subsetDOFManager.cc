#include "discretization/subsetDOFManager.h"

#include "ProfilerApp.h"
#include "ampmesh/MultiIterator.h"
#include "utils/Utilities.h"


namespace AMP {
namespace Discretization {


/****************************************************************
* Constructors                                                  *
****************************************************************/
DOFManager::shared_ptr subsetDOFManager::create( AMP::shared_ptr<const DOFManager> parentDOFManager,
                                                 const std::vector<size_t> &dofs,
                                                 const AMP::Mesh::MeshIterator &iterator,
                                                 const AMP_MPI &comm_in )
{
    // Limit the new comm to be <= the parent comm
    if ( parentDOFManager.get() == nullptr || comm_in.isNull() )
        return DOFManager::shared_ptr();
    PROFILE_START( "subsetDOFManager", 2 );
    AMP_MPI comm = AMP_MPI::intersect( parentDOFManager->getComm(), comm_in );
    // Set the basic info
    AMP::shared_ptr<subsetDOFManager> subsetDOF( new subsetDOFManager() );
    subsetDOF->d_comm             = comm;
    subsetDOF->d_iterator         = iterator;
    subsetDOF->d_parentDOFManager = parentDOFManager;
    // Get the parent DOFs
    subsetDOF->d_parentBegin  = parentDOFManager->beginDOF();
    subsetDOF->d_parentEnd    = parentDOFManager->endDOF();
    subsetDOF->d_parentGlobal = parentDOFManager->numGlobalDOF();
    // Copy the local list of DOFs
    subsetDOF->d_localDOFs.reserve( dofs.size() );
    size_t begin_dof = parentDOFManager->beginDOF();
    size_t end_dof   = parentDOFManager->endDOF();
    for ( auto &dof : dofs ) {
        if ( dof >= begin_dof && dof < end_dof )
            subsetDOF->d_localDOFs.push_back( dof );
    }
    AMP::Utilities::unique( subsetDOF->d_localDOFs );
    // Get the begin and global DOFs for the subset
    size_t N_local = dofs.size();
    subsetDOF->d_comm.sumScan( &N_local, &( subsetDOF->d_end ), 1 );
    subsetDOF->d_begin = subsetDOF->d_end - N_local;
    subsetDOF->d_global =
        subsetDOF->d_comm.bcast( subsetDOF->d_end, subsetDOF->d_comm.getSize() - 1 );
    // Return if the subset DOF is empty
    if ( subsetDOF->d_global == 0 ) {
        PROFILE_STOP2( "subsetDOFManager", 2 );
        return DOFManager::shared_ptr();
    }
    // Return if the subset DOF == parent DOF
    if ( subsetDOF->d_global == parentDOFManager->numGlobalDOF() ) {
        PROFILE_STOP2( "subsetDOFManager", 2 );
        return AMP::const_pointer_cast<DOFManager>( parentDOFManager );
    }
    // Determine which remote DOFs we will need to keep
    size_t *send_data = nullptr;
    if ( N_local > 0 )
        send_data = &( subsetDOF->d_localDOFs[0] );
    std::vector<int> N_remote( subsetDOF->d_comm.getSize(), 0 );
    std::vector<int> N_disp( subsetDOF->d_comm.getSize(), 0 );
    std::vector<size_t> recv_data( subsetDOF->d_global );
    subsetDOF->d_comm.allGather( (int) N_local, &N_remote[0] );
    N_disp[0] = 0;
    for ( int i   = 1; i < subsetDOF->d_comm.getSize(); i++ )
        N_disp[i] = N_disp[i - 1] + N_remote[i - 1];
    subsetDOF->d_comm.allGather(
        send_data, (int) N_local, &recv_data[0], &N_remote[0], &N_disp[0], true );
    AMP::Utilities::quicksort( recv_data );
    std::vector<size_t> remoteDOFs = subsetDOF->d_parentDOFManager->getRemoteDOFs();
    subsetDOF->d_remoteParentDOFs  = std::vector<size_t>();
    subsetDOF->d_remoteSubsetDOFs  = std::vector<size_t>();
    subsetDOF->d_remoteParentDOFs.reserve( remoteDOFs.size() );
    subsetDOF->d_remoteSubsetDOFs.reserve( remoteDOFs.size() );
    size_t k = 0;
    for ( auto &remoteDOF : remoteDOFs ) {
        size_t index = AMP::Utilities::findfirst( recv_data, remoteDOF );
        if ( index == recv_data.size() ) {
            index--;
        }
        if ( recv_data[index] == remoteDOF ) {
            subsetDOF->d_remoteParentDOFs.push_back( remoteDOF );
            subsetDOF->d_remoteSubsetDOFs.push_back( index );
            k++;
        }
    }
    PROFILE_STOP( "subsetDOFManager", 2 );
    if ( subsetDOF->numGlobalDOF() == 0 )
        return DOFManager::shared_ptr();
    return subsetDOF;
}


/****************************************************************
* Deconstructor                                                 *
****************************************************************/
subsetDOFManager::~subsetDOFManager() {}


/****************************************************************
* Get the dofs for the given element                            *
****************************************************************/
void subsetDOFManager::getDOFs( const std::vector<AMP::Mesh::MeshElementID> &ids,
                                std::vector<size_t> &dofs ) const
{
    // Get the parent DOFs
    std::vector<size_t> parentDOFs;
    d_parentDOFManager->getDOFs( ids, parentDOFs );
    if ( parentDOFs.empty() ) {
        dofs.resize( 0 );
        return;
    }
    // Get the subset DOFs
    std::vector<size_t> subsetDOFs = getSubsetDOF( parentDOFs );
    // Remove any DOFs == -1
    std::vector<size_t>::iterator cur = subsetDOFs.begin();
    while ( cur != subsetDOFs.end() ) {
        if ( *cur >= d_global ) {
            cur = subsetDOFs.erase( cur );
        } else {
            ++cur;
        }
    }
    dofs.resize( subsetDOFs.size() );
    for ( size_t i = 0; i < subsetDOFs.size(); i++ )
        dofs[i]    = subsetDOFs[i];
}
void subsetDOFManager::getDOFs( const AMP::Mesh::MeshElementID &id,
                                std::vector<size_t> &dofs ) const
{
    getDOFs( std::vector<AMP::Mesh::MeshElementID>( 1, id ), dofs );
}

/****************************************************************
* Get the element ID give a dof                                 *
****************************************************************/
AMP::Mesh::MeshElement subsetDOFManager::getElement( size_t dof ) const
{
    std::vector<size_t> dof2 = getParentDOF( { dof } );
    return d_parentDOFManager->getElement( dof2[0] );
}


/****************************************************************
* Get an entry over the mesh elements associated with the DOFs  *
* Note: if any sub-DOFManagers are the same, then this will     *
* iterate over repeated elements.                               *
****************************************************************/
AMP::Mesh::MeshIterator subsetDOFManager::getIterator() const { return d_iterator; }


/****************************************************************
* Return the remote DOFs for a vector                           *
****************************************************************/
std::vector<size_t> subsetDOFManager::getRemoteDOFs() const { return d_remoteSubsetDOFs; }


/****************************************************************
* Return the global number of D.O.F.s                           *
****************************************************************/
std::vector<size_t> subsetDOFManager::getRowDOFs( const AMP::Mesh::MeshElement &obj ) const
{
    std::vector<size_t> parentDOFs = d_parentDOFManager->getRowDOFs( obj );
    std::vector<size_t> subsetDOFs = getSubsetDOF( parentDOFs );
    size_t index                   = 0;
    for ( size_t i = 0; i < subsetDOFs.size(); i++ ) {
        if ( subsetDOFs[i] < d_global ) {
            subsetDOFs[index] = subsetDOFs[i];
            index++;
        }
    }
    subsetDOFs.resize( index );
    return subsetDOFs;
}


/****************************************************************
* Function to convert DOFs                                      *
****************************************************************/
std::vector<size_t> subsetDOFManager::getParentDOF( const std::vector<size_t> &subsetDOFs ) const
{
    std::vector<size_t> parentDOFs( subsetDOFs.size() );
    for ( size_t i = 0; i < subsetDOFs.size(); i++ ) {
        size_t DOF = subsetDOFs[i];
        AMP_ASSERT( DOF < d_global );
        if ( DOF >= d_begin && DOF < d_end ) {
            // The DOF is local
            parentDOFs[i] = d_localDOFs[DOF - d_begin];
        } else {
            // The DOF is a remote DOF
            size_t index = AMP::Utilities::findfirst( d_remoteSubsetDOFs, DOF );
            AMP_ASSERT( d_remoteSubsetDOFs[index] == DOF );
            parentDOFs[i] = d_remoteParentDOFs[index];
        }
    }
    return parentDOFs;
}
std::vector<size_t> subsetDOFManager::getSubsetDOF( const std::vector<size_t> &parentDOFs ) const
{
    std::vector<size_t> subsetDOFs( parentDOFs.size(), (size_t) -1 );
    for ( size_t i = 0; i < parentDOFs.size(); i++ ) {
        size_t DOF = parentDOFs[i];
        AMP_ASSERT( DOF < d_parentGlobal );
        if ( DOF >= d_parentBegin && DOF < d_parentEnd ) {
            // The DOF is local
            size_t index = AMP::Utilities::findfirst( d_localDOFs, DOF );
            if ( index == d_localDOFs.size() ) {
                index--;
            }
            if ( d_localDOFs[index] == DOF )
                subsetDOFs[i] = index + d_begin;
        } else if ( !d_remoteParentDOFs.empty() ) {
            // The DOF is a remote DOF
            size_t index = AMP::Utilities::findfirst( d_remoteParentDOFs, DOF );
            if ( index == d_remoteParentDOFs.size() ) {
                index--;
            }
            if ( d_remoteParentDOFs[index] == DOF )
                subsetDOFs[i] = d_remoteSubsetDOFs[index];
        }
    }
    return subsetDOFs;
}
std::vector<size_t> subsetDOFManager::getLocalParentDOFs() const { return d_localDOFs; }


/****************************************************************
* Function to return the DOFManagers                            *
****************************************************************/
AMP::shared_ptr<const DOFManager> subsetDOFManager::getDOFManager() const
{
    return d_parentDOFManager;
}
}
}
