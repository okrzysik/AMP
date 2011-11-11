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


/****************************************************************
* Create a vector                                               *
****************************************************************/
AMP::LinearAlgebra::Vector::shared_ptr DOFManager::createVector( AMP::LinearAlgebra::Variable::shared_ptr variable )
{
    AMP_ERROR("createVector is not implimented for the base class");
    return AMP::LinearAlgebra::Vector::shared_ptr();
}


/****************************************************************
* Create a matrix                                               *
****************************************************************/
AMP::LinearAlgebra::Matrix::shared_ptr DOFManager::createMatrix( 
    AMP::LinearAlgebra::Variable::shared_ptr operand, 
    AMP::LinearAlgebra::Variable::shared_ptr result )
{
    AMP_ERROR("createMatrix is not implimented for the base class");
    return AMP::LinearAlgebra::Matrix::shared_ptr();
}


}
}

