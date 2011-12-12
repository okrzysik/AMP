#include "discretization/DOF_Manager.h"
#include "utils/Utilities.h"


namespace AMP {
namespace Discretization {


/****************************************************************
* Constructors                                                  *
****************************************************************/
DOFManager::DOFManager ( size_t N_local, AMP_MPI comm )
{
    d_comm = comm;
    d_comm.sumScan(&N_local,&d_end,1);
    d_begin = d_end - N_local;
    d_global = d_comm.bcast(d_end,d_comm.getSize()-1);
}


/****************************************************************
* Get the entry indices of DOFs given a mesh element            *
****************************************************************/
void DOFManager::getDOFs( const AMP::Mesh::MeshElement &obj, std::vector <size_t> &dofs, std::vector<size_t> ) const
{
    AMP_ERROR("getDOFs is not implimented for the base class");
}
void DOFManager::getDOFs( const AMP::Mesh::MeshElementID &id, std::vector <size_t> &dofs ) const
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
    return d_begin;
}


/****************************************************************
* Return the last D.O.F. on this core                           *
****************************************************************/
size_t DOFManager::endDOF( ) const
{
    return d_end;
}


/****************************************************************
* Return the local number of D.O.F.s                           *
****************************************************************/
size_t DOFManager::numLocalDOF( ) const
{
    return (d_end-d_begin);
}


/****************************************************************
* Return the global number of D.O.F.s                           *
****************************************************************/
size_t DOFManager::numGlobalDOF( ) const
{
    return d_global;
}


/****************************************************************
* Return the communicator                                       *
****************************************************************/
AMP_MPI DOFManager::getComm( ) const
{
    return d_comm;
}


/****************************************************************
* Return the global number of D.O.F.s                           *
****************************************************************/
std::vector<size_t> DOFManager::getRemoteDOFs( ) const
{
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

