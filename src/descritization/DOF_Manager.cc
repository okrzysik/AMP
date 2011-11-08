#include "descritization/DOF_Manager.h"
#include "utils/Utilities.h"


namespace AMP {
namespace Discretization {


/****************************************************************
* Functions that are not implimented for the base class         *
****************************************************************/
void DOFManager::getDOFs( const AMP::Mesh::MeshElement &obj, std::vector <unsigned int> &ids, unsigned int which ) const
{
    AMP_ERROR("getDOFs is not implimented for the base class");
}


/****************************************************************
* Retrun the first D.O.F. on this core                          *
****************************************************************/
size_t DOFManager::beginDOF( )
{
    AMP_ERROR("beginDOF is not implimented for the base class");
    return 0;
}


/****************************************************************
* Retrun the last D.O.F. on this core                           *
****************************************************************/
size_t DOFManager::endDOF( )
{
    AMP_ERROR("endDOF is not implimented for the base class");
    return 0;
}


/****************************************************************
* Retrun the first D.O.F. on this core                          *
****************************************************************/
AMP::LinearAlgebra::Vector::shared_ptr DOFManager::createVector( AMP::LinearAlgebra::Variable::shared_ptr variable )
{
    AMP_ERROR("createVector is not implimented for the base class");
    return AMP::LinearAlgebra::Vector::shared_ptr();
}


}
}

