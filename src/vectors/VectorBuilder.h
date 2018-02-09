#ifdef USE_AMP_DISCRETIZATION
#ifndef included_AMP_VectorBuider
#define included_AMP_VectorBuider

#include "AMP/discretization/DOF_Manager.h"
#include "AMP/vectors/Vector.h"
#include <string>


namespace AMP {
namespace LinearAlgebra {


/**
 * \brief  This function will create a vector from an arbitrary DOFManager
 * \details  This function is responsible for creating vectors from a DOFManager and variable.
 * \param DOFs      DOFManager to use for constucting the vector
 * \param variable  Variable for the vector
 * \param split     If we are given a multiDOFManager, do we want to split the vector
 *                      based on the individual DOFManagers to create a MultiVector
 */
AMP::LinearAlgebra::Vector::shared_ptr
createVector( AMP::Discretization::DOFManager::shared_ptr DOFs,
              AMP::LinearAlgebra::Variable::shared_ptr variable,
              bool split = true );
} // namespace LinearAlgebra
} // namespace AMP

#endif
#endif
