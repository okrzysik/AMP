
#ifndef included_AMP_Petsc_Matrix
#define included_AMP_Petsc_Matrix

// AMP files
#include "AMP/matrices/petsc/PetscMatrix.h"

namespace AMP::LinearAlgebra {

/** \class NativePetscMatrix
 * \brief  This is a thin wrapper around PETSc Mat
 * \details  As opposed to ManagedPetscMatrix, this is a
 *    thin wrapper around a PETSc Mat.
 */
class NativePetscMatrix : public Matrix
{
protected:
    /** \brief Unused default constructor
     */
    NativePetscMatrix();

    void multiply( shared_ptr other_op, shared_ptr &result ) override;

public:
    /** \brief  Construct a matrix from a PETSc Mat.
     * \param[in] m  The Mat to wrap
     * \param[in] dele  Let this class deallocate the Mat
     */
    explicit NativePetscMatrix( Mat m, bool dele = false );


    /** \brief Destructor
     */
    virtual ~NativePetscMatrix();

    //! Return the type of the matrix
    virtual std::string type() const override { return "NativePetscMatrix"; }

    /** \brief Create a NativePetscMatrix with the same non-zero
     * structure
     * \param[in] m  The matrix to duplicate
     * \return A new matrix with the same non-zero structure
     */
    static std::shared_ptr<Matrix> duplicateMat( Mat m );

    /** \brief Copy data from a PETSc Mat
     * \param[in] m  The matrix with the data
     */
    void copyFromMat( Mat m );

    void mult( Vector::const_shared_ptr in, Vector::shared_ptr out ) override;
    void multTranspose( Vector::const_shared_ptr in, Vector::shared_ptr out ) override;

    shared_ptr cloneMatrix() const override;


    Vector::shared_ptr getRightVector() const override;
    Vector::shared_ptr getLeftVector() const override;
    std::shared_ptr<Discretization::DOFManager> getRightDOFManager() const override;
    std::shared_ptr<Discretization::DOFManager> getLeftDOFManager() const override;

    size_t numGlobalRows() const override;
    size_t numGlobalColumns() const override;

    void scale( double alpha ) override;
    void axpy( double alpha, const Matrix &x ) override;

    void addValuesByGlobalID(
        size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values ) override;
    void setValuesByGlobalID(
        size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values ) override;
    void getValuesByGlobalID( size_t num_rows,
                              size_t num_cols,
                              size_t *rows,
                              size_t *cols,
                              double *values ) const override;
    void getRowByGlobalID( size_t row,
                           std::vector<size_t> &cols,
                           std::vector<double> &values ) const override;

    std::vector<size_t> getColumnIDs( size_t row ) const override;

    void setScalar( double ) override;
    void setDiagonal( Vector::const_shared_ptr in ) override;

    void makeConsistent() override;
    Vector::shared_ptr
    extractDiagonal( Vector::shared_ptr p = Vector::shared_ptr() ) const override;
    double L1Norm() const override;
    void setIdentity() override;
    void zero() override;

private:
    Mat d_Mat;
    bool d_MatCreatedInternally;
};


} // namespace AMP::LinearAlgebra


#endif
