#include "discretization/MultiDOF_Manager.h"

#include "ampmesh/MultiIterator.h"
#include "utils/Utilities.h"


namespace AMP {
namespace Discretization {


/****************************************************************
* Constructors                                                  *
****************************************************************/
multiDOFManager::multiDOFManager ( AMP_MPI globalComm, std::vector<DOFManager::shared_ptr> managers )
{
    d_managers = managers;
    d_comm = globalComm;
    // Compute the total begin, end, and global size
    d_globalSize = std::vector<size_t>(managers.size(),0);
    d_localSize = std::vector<size_t>(managers.size(),0);
    size_t local_size = 0;
    for (size_t i=0; i<managers.size(); i++) {
        d_globalSize[i] = managers[i]->numGlobalDOF();
        d_localSize[i] = managers[i]->numLocalDOF();
        local_size += d_localSize[i];
    }
    d_comm.sumScan(&local_size,&d_end,1);
    d_begin = d_end - local_size;
    d_global = d_comm.bcast(d_end,d_comm.getSize()-1);
    // Compute the relationships between the DOFs
    d_subToGlobalDOF = std::vector< std::vector< std::pair<size_t,size_t> > >(managers.size());
    d_globalToSubDOF = std::vector< std::vector< std::pair<size_t,size_t> > >(managers.size());
    size_t manager_begin = d_begin;
    std::pair<size_t,size_t>  mypair;
    for (size_t i=0; i<managers.size(); i++) {
        AMP_MPI comm = managers[i]->getComm();
        comm.barrier();     // Make sure all necessary processors are at the same point
        mypair.first = managers[i]->beginDOF();
        mypair.second = manager_begin;
        d_subToGlobalDOF[i] = std::vector< std::pair<size_t,size_t> >(comm.getSize());
        comm.allGather(mypair,&d_subToGlobalDOF[i][0]);
        AMP::Utilities::quicksort(d_subToGlobalDOF[i]);
        d_globalToSubDOF[i] = std::vector< std::pair<size_t,size_t> >(comm.getSize());
        for (size_t j=0; j<d_subToGlobalDOF[i].size(); j++) {
            mypair.first = d_subToGlobalDOF[i][j].second;
            mypair.second = d_subToGlobalDOF[i][j].first;
            d_globalToSubDOF[i][j] = mypair;
        }
        AMP::Utilities::quicksort(d_globalToSubDOF[i]);
        manager_begin += d_localSize[i];
    }
}


/****************************************************************
* Get the entry indices of nodal values given a mesh element    *
****************************************************************/
void multiDOFManager::getDOFs( const AMP::Mesh::MeshElement &obj, std::vector <size_t> &dofs, std::vector<size_t> which ) const
{
    dofs.resize(0);
    std::vector<size_t> local_dofs;
    for (size_t i=0; i<d_managers.size(); i++) {
        d_managers[i]->getDOFs( obj, local_dofs, which );
        if ( local_dofs.size() > 0 ) {
            std::vector<size_t> tmp_dofs = getGlobalDOF( d_managers[i], local_dofs );
            dofs.insert(dofs.end(),tmp_dofs.begin(),tmp_dofs.end());
        }
    }
}
void multiDOFManager::getDOFs( const AMP::Mesh::MeshElementID &id, std::vector <size_t> &dofs ) const
{
    dofs.resize(0);
    std::vector<size_t> local_dofs;
    for (size_t i=0; i<d_managers.size(); i++) {
        d_managers[i]->getDOFs( id, local_dofs );
        if ( local_dofs.size() > 0 ) {
            std::vector<size_t> tmp_dofs = getGlobalDOF( d_managers[i], local_dofs );
            dofs.insert(dofs.end(),tmp_dofs.begin(),tmp_dofs.end());
        }
    }
}


/****************************************************************
* Get an entry over the mesh elements associated with the DOFs  *
* Note: if any sub-DOFManagers are the same, then this will     *
* iterate over repeated elements.                               *
****************************************************************/
AMP::Mesh::MeshIterator multiDOFManager::getIterator( ) const
{
    std::vector<boost::shared_ptr<AMP::Mesh::MeshIterator> >  iterators(d_managers.size());
    for (size_t i=0; i<d_managers.size(); i++)
        iterators[i] = boost::shared_ptr<AMP::Mesh::MeshIterator>( new AMP::Mesh::MeshIterator(d_managers[i]->getIterator()) );
    return AMP::Mesh::MultiIterator( iterators );    
}


/****************************************************************
* Return the remote DOFs for a vector                           *
****************************************************************/
std::vector<size_t> multiDOFManager::getRemoteDOFs( ) const
{
    std::vector<size_t> global_dofs;
    std::vector<size_t> local_dofs;
    for (size_t i=0; i<d_managers.size(); i++) {
        d_managers[i]->getRemoteDOFs( );
        if ( local_dofs.size() > 0 ) {
            std::vector<size_t> tmp_dofs = getGlobalDOF( d_managers[i], local_dofs );
            global_dofs.insert(global_dofs.end(),tmp_dofs.begin(),tmp_dofs.end());
        }
    }
    return global_dofs;
}


/****************************************************************
* Return the global number of D.O.F.s                           *
****************************************************************/
std::vector<size_t> multiDOFManager::getRowDOFs( const AMP::Mesh::MeshElement &obj ) const
{
    std::vector<size_t> global_dofs;
    std::vector<size_t> local_dofs;
    for (size_t i=0; i<d_managers.size(); i++) {
        d_managers[i]->getRowDOFs( obj );
        if ( local_dofs.size() > 0 ) {
            std::vector<size_t> tmp_dofs = getGlobalDOF( d_managers[i], local_dofs );
            global_dofs.insert(global_dofs.end(),tmp_dofs.begin(),tmp_dofs.end());
        }
    }
    return global_dofs;
}


/****************************************************************
* Function to convert DOFs                                      *
****************************************************************/
std::vector<size_t> multiDOFManager::getGlobalDOF( DOFManager::shared_ptr manager, std::vector<size_t> &subDOFs ) const
{
    std::vector<size_t> globalDOFs;
    for (size_t i=0; i<d_managers.size(); i++) {
        if ( d_managers[i]==manager ) {
            globalDOFs = std::vector<size_t>(subDOFs.size(),-1);
            std::pair<size_t,size_t> search(0,~((size_t)0));
            for (size_t j=0; j<subDOFs.size(); j++) {
                search.first = subDOFs[j];
                size_t index = AMP::Utilities::findfirst(d_subToGlobalDOF[i],search);
                index--;
                globalDOFs[j] = subDOFs[j] - d_subToGlobalDOF[i][j].first + d_subToGlobalDOF[i][j].second; 
            }
        }
    }
    return globalDOFs;
}
std::vector<size_t> multiDOFManager::getSubDOF( DOFManager::shared_ptr manager, std::vector<size_t> &globalDOFs ) const
{
    std::vector<size_t> subDOFs;
    for (size_t i=0; i<d_managers.size(); i++) {
        if ( d_managers[i]==manager ) {
            subDOFs = std::vector<size_t>(globalDOFs.size(),-1);
            std::pair<size_t,size_t> search(0,~((size_t)0));
            for (size_t j=0; j<globalDOFs.size(); j++) {
                search.first = globalDOFs[i];
                size_t index = AMP::Utilities::findfirst(d_globalToSubDOF[i],search);
                index--;
                subDOFs[j] = globalDOFs[j] - d_subToGlobalDOF[i][j].first + d_globalToSubDOF[i][j].second;
                AMP_ERROR("This is not tested and likely does not work properly");                
                if ( subDOFs[j]<0 || subDOFs[i]>=d_globalSize[i] )
                    subDOFs[j] = ~((size_t)0);
            }
        }
    }
    return subDOFs;
}


/****************************************************************
* Function to return the DOFManagers                            *
****************************************************************/
std::vector<DOFManager::shared_ptr>  multiDOFManager::getDOFManagers() const
{
    return d_managers;
}


}
}

