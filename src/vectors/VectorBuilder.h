#ifdef USE_AMP_DISCRETIZATION
#ifndef included_AMP_VectorBuider
#define included_AMP_VectorBuider

#include <string>
#include "vectors/Vector.h"
#include "discretization/DOF_Manager.h"


namespace AMP {
namespace LinearAlgebra {


/**
 * \brief  This function will create a vector from an arbitrary DOFManager
 * \details  This function is responsible for creating vectors from a DOFManager and variable.
 * \param DOFs      DOFManager to use for constucting the vector
 * \param variable  Variable for the vector
 */
AMP::LinearAlgebra::Vector::shared_ptr  createVector( AMP::Discretization::DOFManager::shared_ptr DOFs, AMP::LinearAlgebra::Variable::shared_ptr variable );


}
}

#endif
#endif
