#ifdef USE_AMP_DISCRETIZATION
#ifndef included_AMP_VectorBuider
#define included_AMP_VectorBuider

#include "AMP/discretization/DOF_Manager.h"
#include "AMP/vectors/Vector.h"
#include <string>

extern "C" {
typedef struct _p_Vec *Vec;
}

namespace AMP {
namespace LinearAlgebra {


/**
 * \brief  This function will create a vector from an arbitrary DOFManager
 * \details  This function is responsible for creating vectors from a DOFManager and variable.
 * \param[in] DOFs          DOFManager to use for constucting the vector
 * \param[in] variable      Variable for the vector
 * \param[in] split         If we are given a multiDOFManager, do we want to split the vector
 *                              based on the individual DOFManagers to create a MultiVector
 */
AMP::LinearAlgebra::Vector::shared_ptr
createVector( AMP::Discretization::DOFManager::shared_ptr DOFs,
              AMP::LinearAlgebra::Variable::shared_ptr variable,
              bool split = true );


#if defined( USE_EXT_PETSC )
/**
 * \brief  Create a vector from an arbitrary PETSc Vec
 * \details  This function creates a vector from an arbitrary PETSc Vec
 * \param[in] v             PETSc Vec
 * \param[in] deleteable    If true, ~Vector() will call VecDestroy()
 * \param[in] comm          The communicator associated with the Vec (optional)
 */
std::shared_ptr<Vector> createVector( Vec v, bool deleteable, AMP_MPI comm = AMP_MPI() );
#endif


} // namespace LinearAlgebra
} // namespace AMP

#endif
#endif
