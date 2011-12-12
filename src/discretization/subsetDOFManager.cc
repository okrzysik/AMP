#include "discretization/subsetDOFManager.h"

#include "ampmesh/MultiIterator.h"
#include "utils/Utilities.h"


namespace AMP {
namespace Discretization {


/****************************************************************
* Constructors                                                  *
****************************************************************/



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
    return getSubsetDOF( d_remoteDOFs );
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
std::vector<size_t> subsetDOFManager::getParentDOF( const std::vector<size_t> &subDOFs ) const
{
    AMP_ERROR("Not programmed");
    return std::vector<size_t>();
}
std::vector<size_t> subsetDOFManager::getSubsetDOF( const std::vector<size_t> &globalDOFs ) const
{
    AMP_ERROR("Not programmed");
    return std::vector<size_t>();
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

