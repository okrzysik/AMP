
#ifndef included_AMP_Petsc_Matrix
#define included_AMP_Petsc_Matrix

// AMP files
#include "matrices/petsc/PetscMatrix.h"

namespace AMP {
namespace LinearAlgebra {

/** \class NativePetscMatrix
  * \brief  This is a thin wrapper around PETSc Mat
  * \details  As opposed to ManagedPetscMatrix, this is a
  *    thin wrapper around a PETSc Mat.
  */
class NativePetscMatrix : public PetscMatrix {
protected:
    /** \brief Unused default constructor
      */
    NativePetscMatrix();

    virtual void multiply( shared_ptr other_op, shared_ptr &result );

public:
    /** \brief  Construct a matrix from a PETSc Mat.
      * \param[in] m  The Mat to wrap
      * \param[in] dele  Let this class deallocate the Mat
      */
    explicit NativePetscMatrix( Mat m, bool dele = false );


    /** \brief Destructor
      */
    virtual ~NativePetscMatrix();

    /** \brief Create a NativePetscMatrix with the same non-zero
      * structure
      * \param[in] m  The matrix to duplicate
      * \return A new matrix with the same non-zero structure
      */
    static Matrix::shared_ptr duplicateMat( Mat m );

    /** \brief Copy data from a PETSc Mat
    * \param[in] m  The matrix with the data
    */
    void copyFromMat( Mat m );

    virtual void mult( Vector::const_shared_ptr in, Vector::shared_ptr out );
    virtual void multTranspose( Vector::const_shared_ptr in, Vector::shared_ptr out );

    virtual shared_ptr cloneMatrix() const;


    virtual Vector::shared_ptr getRightVector() const;
    virtual Vector::shared_ptr getLeftVector() const;
    virtual Discretization::DOFManager::shared_ptr getRightDOFManager() const;
    virtual Discretization::DOFManager::shared_ptr getLeftDOFManager() const;

    virtual size_t numGlobalRows() const;
    virtual size_t numGlobalColumns() const;

    virtual void scale( double alpha );
    virtual void axpy( double alpha, const Matrix &x );

    virtual void
    addValuesByGlobalID( int num_rows, int num_cols, int *rows, int *cols, double *values );
    virtual void
    setValuesByGlobalID( int num_rows, int num_cols, int *rows, int *cols, double *values );
    virtual void
    getValuesByGlobalID( int num_rows, int num_cols, int *rows, int *cols, double *values ) const;
    virtual void
    getRowByGlobalID( int row, std::vector<unsigned int> &cols, std::vector<double> &values ) const;

    virtual void setScalar( double );
    virtual void setDiagonal( Vector::const_shared_ptr in );

    virtual void makeConsistent();
    virtual Vector::shared_ptr extractDiagonal( Vector::shared_ptr p = Vector::shared_ptr() ) const;
    virtual double L1Norm() const;
    virtual void setIdentity();
    virtual void zero();

private:
};
}
}

#include "NativePetscMatrix.inline.h"

#endif
