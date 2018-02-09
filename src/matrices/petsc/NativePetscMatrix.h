
#ifndef included_AMP_Petsc_Matrix
#define included_AMP_Petsc_Matrix

// AMP files
#include "AMP/matrices/petsc/PetscMatrix.h"

namespace AMP {
namespace LinearAlgebra {

/** \class NativePetscMatrix
 * \brief  This is a thin wrapper around PETSc Mat
 * \details  As opposed to ManagedPetscMatrix, this is a
 *    thin wrapper around a PETSc Mat.
 */
class NativePetscMatrix : public PetscMatrix
{
protected:
    /** \brief Unused default constructor
     */
    NativePetscMatrix();

    virtual void multiply( shared_ptr other_op, shared_ptr &result ) override;

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

    virtual void mult( Vector::const_shared_ptr in, Vector::shared_ptr out ) override;
    virtual void multTranspose( Vector::const_shared_ptr in, Vector::shared_ptr out ) override;

    virtual shared_ptr cloneMatrix() const override;


    virtual Vector::shared_ptr getRightVector() const override;
    virtual Vector::shared_ptr getLeftVector() const override;
    virtual Discretization::DOFManager::shared_ptr getRightDOFManager() const override;
    virtual Discretization::DOFManager::shared_ptr getLeftDOFManager() const override;

    virtual size_t numGlobalRows() const override;
    virtual size_t numGlobalColumns() const override;

    virtual void scale( double alpha ) override;
    virtual void axpy( double alpha, const Matrix &x ) override;

    virtual void addValuesByGlobalID(
        size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values ) override;
    virtual void setValuesByGlobalID(
        size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values ) override;
    virtual void getValuesByGlobalID( size_t num_rows,
                                      size_t num_cols,
                                      size_t *rows,
                                      size_t *cols,
                                      double *values ) const override;
    virtual void getRowByGlobalID( size_t row,
                                   std::vector<size_t> &cols,
                                   std::vector<double> &values ) const override;

    std::vector<size_t> getColumnIDs( size_t row ) const override;

    virtual void setScalar( double ) override;
    virtual void setDiagonal( Vector::const_shared_ptr in ) override;

    virtual void makeConsistent() override;
    virtual Vector::shared_ptr
    extractDiagonal( Vector::shared_ptr p = Vector::shared_ptr() ) const override;
    virtual double L1Norm() const override;
    virtual void setIdentity() override;
    virtual void zero() override;

private:
};
} // namespace LinearAlgebra
} // namespace AMP

#include "NativePetscMatrix.inline.h"

#endif
