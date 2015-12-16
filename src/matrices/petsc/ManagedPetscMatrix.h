#ifndef included_AMP_ManagedPetscMatrix
#define included_AMP_ManagedPetscMatrix

// Petsc files
extern "C" {
#include "petscmat.h"
}

// AMP files
#include "matrices/petsc/PetscMatrix.h"
#include "matrices/trilinos/ManagedEpetraMatrix.h"

namespace AMP {
namespace LinearAlgebra {

/** \class ManagedPetscMatrix
  * \brief  A PETSc matrix that allows PETSc to call back into AMP
  * \details  Rather than hold onto a PETSc construct Mat pointer,
  * this class creates a new one and replaces methods in the Mat
  * structure with AMP methods.
  */
class ManagedPetscMatrix : public PetscMatrix, public ManagedEpetraMatrix
{
protected:
    /** \brief Unused constructor
      */
    ManagedPetscMatrix();

    /** \brief Unused copy constructor
      */
    ManagedPetscMatrix( const ManagedPetscMatrix &rhs );

    /** \brief Method that replaces PETSc methods with AMP methods
      */
    void initPetscMat();

public:
    /** \brief Construct a ManagedPetscMatrix from a set of parameters
      * \param[in]  params  The description of the Matrix
      */
    explicit ManagedPetscMatrix( MatrixParameters::shared_ptr params );

    /** \brief Destructor
      */
    virtual ~ManagedPetscMatrix();

    /** \brief Create a NativePetscMatrix with the same non-zero
      * structure
      * \param[in] m  The matrix to duplicate
      * \param[in] comm The communicator to duplicate on
      * \return A new matrix with the same non-zero structure
      */
    static Matrix::shared_ptr duplicateMat( Mat m, AMP_MPI comm );
    /** \brief Copy data from a PETSc Mat
      * \param[in] m  The matrix with the data
      */
    void copyFromMat( Mat m );

    shared_ptr cloneMatrix() const;
};
}
}


#endif
