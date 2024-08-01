#ifndef included_AMP_Petsc_MatrixData
#define included_AMP_Petsc_MatrixData

// AMP files
#include "AMP/matrices/data/MatrixData.h"
#include "petscmat.h"

namespace AMP::LinearAlgebra {

/** \class NativePetscMatrixData
 * \brief  This is a thin wrapper around PETSc Mat
 */
class NativePetscMatrixData : public MatrixData
{
public:
    NativePetscMatrixData();

    explicit NativePetscMatrixData( std::shared_ptr<MatrixParametersBase> params );

    /** \brief  Construct a matrix from a PETSc Mat.
     * \param[in] m  The Mat to wrap
     * \param[in] dele  Let this class deallocate the Mat
     */
    explicit NativePetscMatrixData( Mat m, bool dele = false );

    /** \brief Destructor
     */
    virtual ~NativePetscMatrixData();

    /** \brief Create a NativePetscMatrixData with the same non-zero
     * structure
     * \param[in] m  The matrix to duplicate
     * \return A new matrix with the same non-zero structure
     */
    static std::shared_ptr<MatrixData> duplicateMat( Mat m );

    /** \brief Copy data from a PETSc Mat
     * \param[in] m  The matrix with the data
     */
    void copyFromMat( Mat m );

    //! Return the type of the matrix
    virtual std::string type() const override { return "NativePetscMatrixData"; }

    std::shared_ptr<MatrixData> cloneMatrixData() const override;

    std::shared_ptr<MatrixData> transpose() const override;

    void extractDiagonal( std::shared_ptr<Vector> buf ) const override;

    void addValuesByGlobalID( size_t num_rows,
                              size_t num_cols,
                              size_t *rows,
                              size_t *cols,
                              void *values,
                              const typeID &id ) override;
    void setValuesByGlobalID( size_t num_rows,
                              size_t num_cols,
                              size_t *rows,
                              size_t *cols,
                              void *values,
                              const typeID &id ) override;
    void getValuesByGlobalID( size_t num_rows,
                              size_t num_cols,
                              size_t *rows,
                              size_t *cols,
                              void *values,
                              const typeID &id ) const override;
    void getRowByGlobalID( size_t row,
                           std::vector<size_t> &cols,
                           std::vector<double> &values ) const override;

    std::vector<size_t> getColumnIDs( size_t row ) const override;

    void makeConsistent( AMP::LinearAlgebra::ScatterType t ) override;

    std::shared_ptr<Vector> getRightVector() const;
    std::shared_ptr<Vector> getLeftVector() const;
    std::shared_ptr<Discretization::DOFManager> getRightDOFManager() const override;
    std::shared_ptr<Discretization::DOFManager> getLeftDOFManager() const override;

    size_t numGlobalRows() const override;
    size_t numGlobalColumns() const override;

    Mat getMat() { return d_Mat; }

    void setMat( Mat mat, bool manage = true )
    {
        d_Mat                  = mat;
        d_MatCreatedInternally = manage;
    }

    AMP_MPI getComm() const override;

private:
    Mat d_Mat                   = nullptr;
    bool d_MatCreatedInternally = false;
};


} // namespace AMP::LinearAlgebra


#endif
