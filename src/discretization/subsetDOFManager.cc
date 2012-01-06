#include "discretization/subsetDOFManager.h"

#include "ampmesh/MultiIterator.h"
#include "utils/Utilities.h"


namespace AMP {
namespace Discretization {


/****************************************************************
* Constructors                                                  *
****************************************************************/
subsetDOFManager::subsetDOFManager( DOFManager::shared_ptr parentDOFManager, const std::vector <size_t> &dofs )
{
    d_parentDOFManager = parentDOFManager;
    d_comm = d_parentDOFManager->getComm();
    // Copy the local list of DOFs
    d_localDOFs = dofs;
    AMP::Utilities::quicksort(d_localDOFs);
    size_t N_local = dofs.size();
    d_comm.sumScan(&N_local,&d_end,1);
    d_begin = d_end - N_local;
    d_global = d_comm.bcast(d_end,d_comm.getSize()-1);
    // Get the parent begin, end and global DOFs
    d_parentBegin = d_parentDOFManager->beginDOF();
    d_parentEnd = d_parentDOFManager->endDOF();
    d_parentGlobal = d_parentDOFManager->numGlobalDOF();
    // Return if the subset DOF is empty
    if ( d_global==0 )
        return;
    // Determine which remote DOFs we will need to keep
    size_t *send_data = NULL;
    if ( N_local > 0 )
        send_data = &d_localDOFs[0];
    int *N_remote = new int[d_comm.getSize()];
    int *N_disp = new int[d_comm.getSize()];
    std::vector<size_t> recv_data(d_global);
    d_comm.allGather( (int) N_local, N_remote );
    N_disp[0] = 0;
    for (int i=1; i<d_comm.getSize(); i++)
        N_disp[i] = N_disp[i-1] + N_remote[i-1];
    d_comm.allGather( send_data, (int) N_local, &recv_data[0], N_remote, N_disp, true );
    AMP::Utilities::quicksort( recv_data );
    std::vector<size_t> remoteDOFs = d_parentDOFManager->getRemoteDOFs();
    d_remoteParentDOFs = std::vector<size_t>();
    d_remoteSubsetDOFs = std::vector<size_t>();
    d_remoteParentDOFs.reserve(remoteDOFs.size());
    d_remoteSubsetDOFs.reserve(remoteDOFs.size());
    size_t k = 0;
    for (size_t i=0; i<remoteDOFs.size(); i++ ) {
        size_t index = AMP::Utilities::findfirst(recv_data,remoteDOFs[i]);
        if ( recv_data[index]==remoteDOFs[i] ) {
            d_remoteParentDOFs[k] = remoteDOFs[i];
            d_remoteSubsetDOFs[k] = index;
            k++;
        }
    }
    delete [] N_remote;
    delete [] N_disp;
}


/****************************************************************
* Get the entry indices of nodal values given a mesh element    *
****************************************************************/
void subsetDOFManager::getDOFs( const AMP::Mesh::MeshElement &obj, std::vector <size_t> &dofs, std::vector<size_t> which ) const
{
    std::vector<size_t> parentDOFs;
    d_parentDOFManager->getDOFs( obj, parentDOFs, which );
    std::vector<size_t> subsetDOFs = getSubsetDOF( parentDOFs );
    std::vector<size_t>::iterator cur = subsetDOFs.begin();
    std::vector<size_t>::iterator end = subsetDOFs.end();
    while ( cur != end ) {
        if ( *cur >= d_global )
            subsetDOFs.erase(cur);
        ++cur;
    }
    dofs.resize(subsetDOFs.size());
    for (size_t i=0; i<subsetDOFs.size(); i++)
        dofs[i] = subsetDOFs[i];
}
void subsetDOFManager::getDOFs( const AMP::Mesh::MeshElementID &id, std::vector <size_t> &dofs ) const
{
    std::vector<size_t> parentDOFs;
    d_parentDOFManager->getDOFs( id, parentDOFs );
    std::vector<size_t> subsetDOFs = getSubsetDOF( parentDOFs );
    std::vector<size_t>::iterator cur = subsetDOFs.begin();
    std::vector<size_t>::iterator end = subsetDOFs.end();
    while ( cur != end ) {
        if ( *cur >= d_global )
            subsetDOFs.erase(cur);
        ++cur;
    }
    dofs.resize(subsetDOFs.size());
    for (size_t i=0; i<subsetDOFs.size(); i++)
        dofs[i] = subsetDOFs[i];
}


/****************************************************************
* Get an entry over the mesh elements associated with the DOFs  *
* Note: if any sub-DOFManagers are the same, then this will     *
* iterate over repeated elements.                               *
****************************************************************/
AMP::Mesh::MeshIterator subsetDOFManager::getIterator( ) const
{
    AMP_ERROR("Not programmed yet");
    return AMP::Mesh::MeshIterator();
}


/****************************************************************
* Return the remote DOFs for a vector                           *
****************************************************************/
std::vector<size_t> subsetDOFManager::getRemoteDOFs( ) const
{
    return d_remoteSubsetDOFs;
}


/****************************************************************
* Return the global number of D.O.F.s                           *
****************************************************************/
std::vector<size_t> subsetDOFManager::getRowDOFs( const AMP::Mesh::MeshElement &obj ) const
{
    std::vector<size_t> parentDOFs = d_parentDOFManager->getRowDOFs( obj );
    std::vector<size_t> subsetDOFs = getSubsetDOF( parentDOFs );
    std::vector<size_t>::iterator cur = subsetDOFs.begin();
    std::vector<size_t>::iterator end = subsetDOFs.end();
    while ( cur != end ) {
        if ( *cur >= d_global )
            subsetDOFs.erase(cur);
        ++cur;
    }
    return subsetDOFs;
}


/****************************************************************
* Function to convert DOFs                                      *
****************************************************************/
std::vector<size_t> subsetDOFManager::getParentDOF( const std::vector<size_t> &subsetDOFs ) const
{
    std::vector<size_t> parentDOFs(subsetDOFs.size());
    for (size_t i=0; i<subsetDOFs.size(); i++) {
        size_t DOF = subsetDOFs[i];
        AMP_ASSERT(DOF<d_global);
        if ( DOF>=d_begin && DOF<d_end ) {
            // The DOF is local
            parentDOFs[i] = d_localDOFs[DOF-d_begin];
        } else {
            // The DOF is a remote DOF
            size_t index = AMP::Utilities::findfirst(d_remoteSubsetDOFs,DOF);
            AMP_ASSERT(d_remoteSubsetDOFs[index]==DOF);
            parentDOFs[i] = d_remoteParentDOFs[index];
        }
    }
    return parentDOFs;
}
std::vector<size_t> subsetDOFManager::getSubsetDOF( const std::vector<size_t> &parentDOFs ) const
{
    std::vector<size_t> subsetDOFs(subsetDOFs.size(),(size_t)-1);
    for (size_t i=0; i<parentDOFs.size(); i++) {
        size_t DOF = parentDOFs[i];
        AMP_ASSERT(DOF<d_parentGlobal);
        if ( DOF>=d_parentBegin && DOF<d_parentEnd ) {
            // The DOF is local
            size_t index = AMP::Utilities::findfirst(d_localDOFs,DOF);
            if ( d_localDOFs[index] == DOF )
                subsetDOFs[i] = index + d_begin;
        } else {
            // The DOF is a remote DOF
            size_t index = AMP::Utilities::findfirst(d_remoteParentDOFs,DOF);
            if ( d_remoteParentDOFs[index] == DOF )
                subsetDOFs[i] = d_remoteSubsetDOFs[index];
        }
    }
    return parentDOFs;
}
std::vector<size_t> subsetDOFManager::getLocalParentDOFs( ) const
{
    return d_localDOFs;
}


/****************************************************************
* Function to return the DOFManagers                            *
****************************************************************/
DOFManager::shared_ptr  subsetDOFManager::getDOFManager() const
{
    return d_parentDOFManager;
}


}
}

