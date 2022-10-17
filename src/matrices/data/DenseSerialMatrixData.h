#ifndef included_AMP_DenseSerialMatrixData
#define included_AMP_DenseSerialMatrixData

#include "AMP/matrices/Matrix.h"
#include "AMP/vectors/Vector.h"
#include <memory>

class DenseSerialMatrix;

namespace AMP::LinearAlgebra {


/** \class DenseSerialMatrixData
 * \brief  An concrete class for dealing with dense serial matrices
 * \details  This is a concrete class that stores a dense local matrix.
 *    This is not a distributed matrix and requires that the comm is AMP_COMM_SELF.
 */
class DenseSerialMatrixData : public MatrixData
{
public:
    /** \brief Constructor
     * \param[in] params  Description of the matrix
     */
    explicit DenseSerialMatrixData( std::shared_ptr<MatrixParameters> params );

    DenseSerialMatrixData()                                = delete;
    DenseSerialMatrixData( const DenseSerialMatrixData & ) = delete;
    DenseSerialMatrixData &operator=( const DenseSerialMatrixData & ) = delete;

    /** \brief Destructor
     */
    virtual ~DenseSerialMatrixData();

    //! Return the type of the matrix
    std::string type() const override { return "DenseSerialMatrixData"; }

    /** \brief  Return a matrix with the same non-zero and distributed structure
     * \return  The new matrix
     */
    std::shared_ptr<MatrixData> cloneMatrixData() const override;

    /** \brief  Add values to those in the matrix
     * \param[in] num_rows The number of rows represented in values
     * \param[in] num_cols The number of cols represented in values
     * \param[in] rows  The row ids of values
     * \param[in] cols  The column ids of values
     * \param[in] values  The values to add to the matrix
     * \details  This method may fail if the matrix has not
     * allocated a particular(row,col) specified, depending
     * on the actual subclass of matrix used.
     */
    void addValuesByGlobalID(
        size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values ) override;

    /** \brief  Set values in the matrix
     * \param[in] num_rows The number of rows represented in values
     * \param[in] num_cols The number of cols represented in values
     * \param[in] rows  The row ids of values
     * \param[in] cols  The column ids of values
     * \param[in] values  The values to set to the matrix
     * \details  This method may fail if the matrix has not
     * allocated a particular(row,col) specified, depending
     * on the actual subclass of matrix used.
     */
    void setValuesByGlobalID(
        size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values ) override;

    /** \brief  Get values in the matrix
     * \param[in] num_rows The number of rows represented in values
     * \param[in] num_cols The number of cols represented in values
     * \param[in] rows  The row ids of values
     * \param[in] cols  The column ids of values
     * \param[in] values  The values to get from the matrix (row-major ordering)
     * \details  This method will return zero for any entries that
     *   have not been allocated or are not ghosts on the current processor.
     */
    void getValuesByGlobalID( size_t num_rows,
                              size_t num_cols,
                              size_t *rows,
                              size_t *cols,
                              double *values ) const override;


    /** \brief  Retrieve a row of the matrix in compressed format
     * \param[in]  row Which row
     * \param[out] cols  The column ids of the returned values
     * \param[out] values  The values in the row
     */
    void getRowByGlobalID( size_t row,
                           std::vector<size_t> &cols,
                           std::vector<double> &values ) const override;

    /** \brief  Given a row, retrieve the non-zero column indices of the matrix in compressed format
     * \param[in]  row Which row
     */
    std::vector<size_t> getColumnIDs( size_t row ) const override;

    /** \brief  Perform communication to ensure values in the
     * matrix are the same across cores.
     */
    void makeConsistent() override {}

    /** \brief Get the DOFManager associated with a right vector( For \f$\mathbf{y}^T\mathbf{Ax}\f$,
     * \f$\mathbf{x}\f$ is
     * a right vector )
     * \return  The DOFManager associated with a right vector
     */
    Discretization::DOFManager::shared_ptr getRightDOFManager() const override;

    /** \brief Get the DOFManager associated with a left vector( For \f$\mathbf{y}^T\mathbf{Ax}\f$,
     * \f$\mathbf{y}\f$ is
     * a left vector )
     * \return  The DOFManager associated with a left vector
     */
    Discretization::DOFManager::shared_ptr getLeftDOFManager() const override;

    /** \brief  Get the number of local rows in the matrix
     * \return  The number of local rows
     */
    size_t numLocalRows() const override { return d_rows; }

    /** \brief  Get the number of global rows in the matrix
     * \return  The number of global rows
     */
    size_t numGlobalRows() const override { return d_rows; }

    /** \brief  Get the number of local columns in the matrix
     * \return  The number of local columns
     */
    size_t numLocalColumns() const override { return d_cols; }

    /** \brief  Get the number of global columns in the matrix
     * \return  The number of global columns
     */
    size_t numGlobalColumns() const override { return d_cols; }

    std::shared_ptr<AMP::LinearAlgebra::Variable> getLeftVariable() { return d_VariableLeft; }
    std::shared_ptr<AMP::LinearAlgebra::Variable> getRightVariable() { return d_VariableRight; }

protected:
    // AMP variables and DOFManagers for the left and right vectors
    std::shared_ptr<AMP::LinearAlgebra::Variable> d_VariableLeft;
    std::shared_ptr<AMP::LinearAlgebra::Variable> d_VariableRight;
    std::shared_ptr<AMP::Discretization::DOFManager> d_DOFManagerLeft;
    std::shared_ptr<AMP::Discretization::DOFManager> d_DOFManagerRight;

    // Data for the matrix
    size_t d_rows;
    size_t d_cols;
    double *d_M;

    friend class DenseSerialMatrix;
};
} // namespace AMP::LinearAlgebra


#endif
