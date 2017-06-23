#ifndef included_AMP_PetscMatrix
#define included_AMP_PetscMatrix

#include "vectors/petsc/PetscHelpers.h"


#ifdef MPICH_SKIP_MPICXX
#define _FIX_FOR_PETSC_MPICH_CXX
#undef MPICH_SKIP_MPICXX
#endif

#ifdef OMPI_SKIP_MPICXX
#define _FIX_FOR_PETSC_OMPI_CXX
#undef OMPI_SKIP_MPICXX
#endif

#include "petscmat.h"


#ifdef _FIX_FOR_PETSC_MPICH_CXX
#ifndef MPICH_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#endif
#endif

#ifdef _FIX_FOR_PETSC_OMPI_CXX
#ifndef OMPI_SKIP_MPICXX
#define OMPI_SKIP_MPICXX
#endif
#endif


#include "matrices/Matrix.h"

namespace AMP {
namespace LinearAlgebra {

/** \brief  Parameters to create a PetscMatrix
  */
typedef MatrixParameters PetscMatrixParameters;

/** \class  PetscMatrix
  * \brief  A matrix that can provide a PETSc Mat
  * \details  A PetscMatrix presents a Mat data structure.
  * Given an AMP::LinearAlgebra::Matrix, this class can create a Mat view
  * without copying the data.  As such, this class serves three
  * purposes:
  *  -# Provides a Mat for derived classes to use, fill, manage, etc.
  *  -# Provides an interface for accessing this Mat independent of base or derived classes
  *  -# Provides a static method for creating a Mat view of an AMP vector.
  */
class PetscMatrix : virtual public Matrix
{
protected:
    /** \brief Unused default constrcutor
      */
    PetscMatrix();

    /** \brief Unused copy constructor
      */
    PetscMatrix( const PetscMatrix &rhs );

    /** \brief Indicates if d_Mat was created internally
      */
    bool d_MatCreatedInternally;

    /** \brief  The Mat used by inherited classes
      */
    Mat d_Mat;

public:
    /** \brief Create a PetscMatrix from a set of parameters
      * \param[in] params  Description of the PetscMatrix
      */
    explicit PetscMatrix( MatrixParameters::shared_ptr params );

    /** \brief  Destructor
      */
    virtual ~PetscMatrix();

    /** \brief  Get the underlying PETSc Mat
      * \return a PETSc Mat view of this matrix
      */
    virtual Mat &getMat();

    /** \brief  Get the underlying PETSc Mat
      * \return a PETSc Mat view of this matrix
      */
    virtual Mat getMat() const;

    /** \brief  Create a view of an AMP::LinearAlgebra::Matrix that has
      * a PETSc Mat view
      * \return An AMP::LinearAlgebra::Matrix guaranteed to be a derived
      * class of PetscMatrix.
      */
    static shared_ptr createView( shared_ptr m );

};
}
}

#include "PetscMatrix.inline.h"

#endif
