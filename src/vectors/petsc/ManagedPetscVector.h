#ifndef included_AMP_ManagedPetscVector
#define included_AMP_ManagedPetscVector

#include "AMP/vectors/data/DataChangeListener.h"
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
class ManagedPetscVector : public Vector, public PetscVector, public DataChangeListener
{
private:
    std::shared_ptr<PETSC::PetscVectorWrapper> d_wrapper;

public:
    /** \brief Construct a view of another vector
     * \param[in] alias The vector to view
     */
    explicit ManagedPetscVector( Vector::shared_ptr alias );

    /** \brief Method to create a duplicate of this vector for VecDuplicate
     * \return Raw pointer to a new vector.  This does not copy data
     */
    ManagedPetscVector *petscDuplicate();

    /** \brief Identifies whether this vector was created through the
     * VecDuplicate interface
     * \return true if constructed with VecDuplicate.  False otherwise
     */
    bool constructedWithPetscDuplicate();

    /** \brief Destructor
     */
    virtual ~ManagedPetscVector();

    /** \brief Create an exact clone of this vector.
     * \return A raw pointer to a clone of this vector
     */
    ManagedPetscVector *rawClone() const;

    //! Check if the two vectors are alias of each other
    bool isAnAliasOf( const ManagedPetscVector &rhs ) const;

    //! Get the PETSc vector
    Vec &getVec() override { return d_wrapper->getVec(); }

    //! Get the PETSc vector
    const Vec &getVec() const override { return d_wrapper->getVec(); }

public: // These are adequately documented in a base class
    void swapVectors( Vector &other ) override;
    using Vector::cloneVector;
    Vector::shared_ptr cloneVector( const Variable::shared_ptr p ) const override;

    Vector::shared_ptr subsetVectorForVariable( Variable::const_shared_ptr name ) override;
    Vector::const_shared_ptr
    constSubsetVectorForVariable( Variable::const_shared_ptr name ) const override;

    virtual std::string type() const override
    {
        return "Managed PETSc Vector" + d_VectorData->VectorDataName();
    }

    virtual bool petscHoldsView() const override;

    void receiveDataChanged() override;
};


} // namespace LinearAlgebra
} // namespace AMP

#endif
