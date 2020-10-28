#ifndef included_AMP_ManagedPetscVector
#define included_AMP_ManagedPetscVector

#include "AMP/vectors/Vector.h"
#include "AMP/vectors/petsc/PetscHelpers.h"
#include "AMP/vectors/petsc/PetscVector.h"


namespace AMP {
namespace LinearAlgebra {


/** \class
 * ManagedPetscVector/projects/AMP/build/debug/AMP/include/AMP/vectors/petsc/ManagedPetscVector.h
 * \brief A class that provides a PETSc vector interfaced to a ManagedVector.
 * \details  This class provides a PETSc Vec specially configured to call
 * through ManagedVector.
 *
 * In general, the end user will not have to use this class.  This class is
 * returned by PetscVector::view() and PetscVector::constView();
 *
 * \see PetscVector
 */
class ManagedPetscVector : public Vector, public PetscVector
{
public:
    /** \brief Construct a view of another vector
     * \param[in] alias The vector to view
     */
    explicit ManagedPetscVector( Vector::shared_ptr alias );

    /** \brief Destructor
     */
    virtual ~ManagedPetscVector();

public: // These are adequately documented in a base class
    void swapVectors( Vector &other ) override;
    std::unique_ptr<Vector> rawClone( const Variable::shared_ptr p ) const override;

    Vector::shared_ptr subsetVectorForVariable( Variable::const_shared_ptr name ) override;
    Vector::const_shared_ptr
    constSubsetVectorForVariable( Variable::const_shared_ptr name ) const override;

    std::string type() const override
    {
        return "Managed PETSc Vector" + d_VectorData->VectorDataName();
    }


    std::shared_ptr<Vector> getManagedVec() override { return shared_from_this(); }
    std::shared_ptr<const Vector> getManagedVec() const override { return shared_from_this(); }
};


} // namespace LinearAlgebra
} // namespace AMP

#endif
