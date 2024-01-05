#ifndef included_AMP_CSRMatrixData_h
#define included_AMP_CSRMatrixData_h

#include "AMP/matrices/data/MatrixData.h"

#include <tuple>

namespace AMP::Discretization {
class DOFManager;
}

namespace AMP::LinearAlgebra {

template<typename Policy>
class CSRMatrixData : public MatrixData
{
public:
    using gidx_t   = typename Policy::gidx_t;
    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    /** \brief Constructor
     * \param[in] params  Description of the matrix
     */
    explicit CSRMatrixData( std::shared_ptr<MatrixParametersBase> params );

    //! Destructor
    virtual ~CSRMatrixData();

    //! Empty constructor
    CSRMatrixData();

    //! Copy constructor
    CSRMatrixData( const CSRMatrixData & ) = delete;

    //! Clone the data
    std::shared_ptr<MatrixData> cloneMatrixData() const override;

    //! Transpose
    std::shared_ptr<MatrixData> transpose() const override;

    //! Extract the diagonal vector
    void extractDiagonal( std::shared_ptr<Vector> buf ) const override;

    //! Return the type of the matrix
    std::string type() const override { return "CSRMatrixData"; }

    /** \brief  Retrieve a row of the matrix in compressed format
     * \param[in]  row Which row
     * \param[out] cols  The column ids of the returned values
     * \param[out] values  The values in the row
     */
    void getRowByGlobalID( size_t row,
                           std::vector<size_t> &cols,
                           std::vector<double> &values ) const override;

    /** \brief  Add values to those in the matrix
     * \param[in] num_rows The number of rows represented in values
     * \param[in] num_cols The number of cols represented in values
     * \param[in] rows  The row ids of values
     * \param[in] cols  The column ids of values
     * \param[in] values  The values to add to the matrix (row-major ordering)
     * \details  This method may fail if the matrix has not
     * allocated a particular (row,col) specified, depending
     * on the actual subclass of matrix used.
     */
    void addValuesByGlobalID( size_t num_rows,
                              size_t num_cols,
                              size_t *rows,
                              size_t *cols,
                              void *values,
                              const typeID &id ) override;

    /** \brief  Set values in the matrix
     * \param[in] num_rows The number of rows represented in values
     * \param[in] num_cols The number of cols represented in values
     * \param[in] rows  The row ids of values
     * \param[in] cols  The column ids of values
     * \param[in] values  The values to set to the matrix (row-major ordering)
     * \details  This method may fail if the matrix has not
     * allocated a particular (row,col) specified, depending
     * on the actual subclass of matrix used.
     */
    void setValuesByGlobalID( size_t num_rows,
                              size_t num_cols,
                              size_t *rows,
                              size_t *cols,
                              void *values,
                              const typeID &id ) override;

    /** \brief  Get values in the matrix
     * \param[in] num_rows The number of rows represented in values
     * \param[in] num_cols The number of cols represented in values
     * \param[in] rows  The row ids of values
     * \param[in] cols  The column ids of values
     * \param[out] values  The values to get from the matrix (row-major ordering)
     * \details  This method will return zero for any entries that
     *   have not been allocated or are not ghosts on the current processor.
     */
    void getValuesByGlobalID( size_t num_rows,
                              size_t num_cols,
                              size_t *rows,
                              size_t *cols,
                              void *values,
                              const typeID &id ) const override;

    /** \brief  Given a row, retrieve the non-zero column indices of the matrix in compressed format
     * \param[in]  row Which row
     */
    std::vector<size_t> getColumnIDs( size_t row ) const override;

    /** \brief  Perform communication to ensure values in the
     * matrix are the same across cores.
     */
    void makeConsistent() override;

    /** \brief Get the DOFManager associated with a right vector ( For
     * \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$\mathbf{x}\f$
     * is a right vector )
     * \return  The DOFManager associated with a right vector
     */
    std::shared_ptr<Discretization::DOFManager> getRightDOFManager() const override;

    /** \brief Get the DOFManager associated with a left vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$,
     * \f$\mathbf{y}\f$ is
     * a left vector )
     * \return  The DOFManager associated with a left vector
     */
    std::shared_ptr<Discretization::DOFManager> getLeftDOFManager() const override;

    /** \brief  Get the number of local rows in the matrix
     * \return  The number of local rows
     */
    size_t numLocalRows() const override;

    /** \brief  Get the number of global rows in the matrix
     * \return  The number of global rows
     */
    size_t numGlobalRows() const override;

    /** \brief  Get the number of local columns in the matrix
     * \return  The number of local columns
     */
    size_t numLocalColumns() const override;

    /** \brief  Get the number of global columns in the matrix
     * \return  The number of global columns
     */
    size_t numGlobalColumns() const override;

    /** \brief  Get the global id of the beginning row
     * \return  beginning global row id
     */
    size_t beginRow() const override;

    /** \brief  Get the global id of the ending row
     * \return  ending global row id
     */
    size_t endRow() const override;

    std::tuple<lidx_t *, gidx_t const *, scalar_t const *> getCSRData()
    {
        return std::make_tuple( d_nnz_per_row, d_cols, d_coeffs );
    }

    bool isSquare() const noexcept { return d_is_square; }

protected:
    bool d_is_square;
    gidx_t d_first_row;
    gidx_t d_last_row;
    gidx_t d_first_col;
    gidx_t d_last_col;
    lidx_t *d_nnz_per_row;
    gidx_t const *d_cols;
    scalar_t const *d_coeffs;

    std::shared_ptr<Discretization::DOFManager> d_leftDOFManager;
    std::shared_ptr<Discretization::DOFManager> d_rightDOFManager;
};

} // namespace AMP::LinearAlgebra

#endif
