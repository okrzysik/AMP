#ifndef included_AMP_NativePetscVector
#define included_AMP_NativePetscVector

#include "AMP/utils/AMP_MPI.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/petsc/PetscVector.h"
#include "AMP/vectors/petsc/NativePetscVectorData.h"

namespace AMP {
namespace LinearAlgebra {


/*! \struct Vec
    \brief PETSc vector
*/


/** \class NativePetscVector
 * \brief An AMP Vector that uses PETSc for parallel data management, linear algebra,
 * etc.
 * \details  This is an AMP wrapper to PETSc.  This is different from ManagedPetscVector
 * in that this class does not replace calls to Vec*.  Rather, it wraps these calls.
 * This class is used when PETSc is chosen as the default linear algebra engine.
 *
 * This class is not to be used directly, just through base class interfaces.
 * \see PetscVector
 * \see ManagedPetscVector
 */
class NativePetscVector : public Vector
{
public:
    /** \brief Construct a wrapper for a PETSc Vec from a set of parameters
     * \param[in] params The parameters describing the Vec
     */
    explicit NativePetscVector( VectorParameters::shared_ptr params );
    explicit NativePetscVector( std::shared_ptr<VectorData> data );

    //! Destructor
    virtual ~NativePetscVector();

    std::string type() const override { return "Native PETSc Vector"; }

    using Vector::cloneVector;
    Vector::shared_ptr cloneVector( const Variable::shared_ptr ) const override;
    void swapVectors( Vector &other ) override;
    void aliasVector( Vector & ) override;
    void assemble() override;

    std::shared_ptr<ParameterBase> getParameters() override;

};


} // namespace LinearAlgebra
} // namespace AMP

#endif
