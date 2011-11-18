#include "discretization/DOF_Manager.h"
#include "utils/Utilities.h"


namespace AMP {
namespace Discretization {


/****************************************************************
* Get the entry indices of DOFs given a mesh element            *
****************************************************************/
void DOFManager::getDOFs( const AMP::Mesh::MeshElement &obj, std::vector <unsigned int> &ids, std::vector<unsigned int> ) const
{
    AMP_ERROR("getDOFs is not implimented for the base class");
}


/****************************************************************
* Get an entry over the mesh elements associated with the DOFs  *
****************************************************************/
AMP::Mesh::MeshIterator DOFManager::getIterator( ) const
{
    AMP_ERROR("getIterator is not implimented for the base class");
    return AMP::Mesh::MeshIterator();
}


/****************************************************************
* Return the first D.O.F. on this core                          *
****************************************************************/
size_t DOFManager::beginDOF( ) const
{
    AMP_ERROR("beginDOF is not implimented for the base class");
    return 0;
}


/****************************************************************
* Return the last D.O.F. on this core                           *
****************************************************************/
size_t DOFManager::endDOF( ) const
{
    AMP_ERROR("endDOF is not implimented for the base class");
    return 0;
}


/****************************************************************
* Return the local number of D.O.F.s                           *
****************************************************************/
size_t DOFManager::numLocalDOF( ) const
{
    AMP_ERROR("numGlobalDOF is not implimented for the base class");
    return 0;
}


/****************************************************************
* Return the global number of D.O.F.s                           *
****************************************************************/
size_t DOFManager::numGlobalDOF( ) const
{
    AMP_ERROR("numGlobalDOF is not implimented for the base class");
    return 0;
}


/****************************************************************
* Return the global number of D.O.F.s                           *
****************************************************************/
AMP_MPI DOFManager::getComm( ) const
{
    AMP_ERROR("getComm is not implimented for the base class");
    return AMP_MPI(AMP_COMM_NULL);
}


/****************************************************************
* Return the global number of D.O.F.s                           *
****************************************************************/
std::vector<size_t> DOFManager::getRemoteDOFs( ) const
{
    AMP_ERROR("getRemoteDOFs is not implimented for the base class");
    return std::vector<size_t>();
}


/****************************************************************
* Return the global number of D.O.F.s                           *
****************************************************************/
std::vector<size_t> DOFManager::getRowDOFs( const AMP::Mesh::MeshElement &obj ) const
{
    AMP_ERROR("getRowDOFs is not implimented for the base class");
    return std::vector<size_t>();
}


}
}

