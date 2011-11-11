#include "descritization/DOF_Manager.h"
#include "utils/Utilities.h"


namespace AMP {
namespace Discretization {


/****************************************************************
* Functions that are not implimented for the base class         *
****************************************************************/
void DOFManager::getDOFs( const AMP::Mesh::MeshElement &obj, std::vector <unsigned int> &ids, std::vector<unsigned int> ) const
{
    AMP_ERROR("getDOFs is not implimented for the base class");
}


/****************************************************************
* Return the first D.O.F. on this core                          *
****************************************************************/
size_t DOFManager::beginDOF( )
{
    AMP_ERROR("beginDOF is not implimented for the base class");
    return 0;
}


/****************************************************************
* Return the last D.O.F. on this core                           *
****************************************************************/
size_t DOFManager::endDOF( )
{
    AMP_ERROR("endDOF is not implimented for the base class");
    return 0;
}


/****************************************************************
* Return the local number of D.O.F.s                           *
****************************************************************/
size_t DOFManager::numLocalDOF( )
{
    AMP_ERROR("numGlobalDOF is not implimented for the base class");
    return 0;
}


/****************************************************************
* Return the global number of D.O.F.s                           *
****************************************************************/
size_t DOFManager::numGlobalDOF( )
{
    AMP_ERROR("numGlobalDOF is not implimented for the base class");
    return 0;
}


}
}

