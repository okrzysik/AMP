#ifndef included_AMP_CSRSerialMatrixData_h
#define included_AMP_CSRSerialMatrixData_h

#include <tuple>

#include "AMP/utils/enable_shared_from_this.h"
#include "AMP/utils/typeid.h"

namespace AMP::Discretization {
class DOFManager;
}

namespace AMP::LinearAlgebra {

// This intentionally does not inherit from matrix data
template<typename Policy>
class CSRSerialMatrixData : public AMP::enable_shared_from_this<CSRSerialMatrixData<Policy>>
{
public:
    using gidx_t   = typename Policy::gidx_t;
    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    /** \brief Constructor
     * \param[in] params  Description of the matrix
     */
    explicit CSRSerialMatrixData( std::shared_ptr<MatrixParametersBase> params,
				  bool is_diag );

    //! Destructor
    virtual ~CSRSerialMatrixData();

    //! Empty constructor
    CSRSerialMatrixData();

    //! Copy constructor
    CSRSerialMatrixData( const CSRSerialMatrixData & ) = delete;

    //! Clone the data
    std::shared_ptr<MatrixData> cloneMatrixData() const;

    //! Transpose
    std::shared_ptr<MatrixData> transpose() const;

    //! Extract the diagonal vector
    void extractDiagonal( std::shared_ptr<Vector> buf ) const;

    //! Return the type of the matrix
    std::string type() const { return "CSRSerialMatrixData"; }

    /** \brief  Retrieve a row of the matrix in compressed format
     * \param[in]  row Which row
     * \param[out] cols  The column ids of the returned values
     * \param[out] values  The values in the row
     */
    void getRowByGlobalID( size_t row,
                           std::vector<size_t> &cols,
                           std::vector<double> &values ) const;

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
                              const typeID &id );

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
                              const typeID &id );

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
                              const typeID &id ) const;

    /** \brief  Given a row, retrieve the non-zero column indices of the matrix in compressed format
     * \param[in]  row Which row
     */
    std::vector<size_t> getColumnIDs( size_t row ) const;

    /** \brief  Get the number of local rows in the matrix
     * \return  The number of local rows
     */
    size_t numLocalRows() const;

    /** \brief  Get the number of local columns in the matrix
     * \return  The number of local columns
     */
    size_t numLocalColumns() const;

    /** \brief  Get the global id of the beginning row
     * \return  beginning global row id
     */
    size_t beginRow() const;

    /** \brief  Get the global id of the ending row
     * \return  ending global row id
     */
    size_t endRow() const;

    size_t beginCol() const { return d_first_col; }

    std::tuple<lidx_t *, lidx_t const *, scalar_t const *> getCSRData()
    {
        return std::make_tuple( d_nnz_per_row, d_cols_loc, d_coeffs );
    }

    bool isDiag() const noexcept { return d_is_diag; }

    bool isEmpty() const { return d_is_empty; }

    auto numberOfNonZeros() const { return d_nnz; }

protected:
    bool d_is_diag     = true;
    bool d_is_empty    = false;
    gidx_t d_first_row = 0;
    gidx_t d_last_row  = 0;
    gidx_t d_first_col = 0;
    gidx_t d_last_col  = 0;

    lidx_t *d_nnz_per_row = nullptr;
    lidx_t *d_row_starts  = nullptr;
    lidx_t *d_cols_loc    = nullptr;
    gidx_t *d_cols        = nullptr;
    scalar_t *d_coeffs    = nullptr;

    lidx_t d_nnz = 0;

    AMP::Utilities::MemoryType d_memory_location = AMP::Utilities::MemoryType::host;

    std::shared_ptr<MatrixParametersBase> d_pParameters;

    bool d_manage_cols   = true;
    bool d_manage_nnz    = true;
    bool d_manage_coeffs = true;
};
  
} // namespace AMP::LinearAlgebra

#endif

#include "AMP/matrices/data/CSRSerialMatrixData.hpp"
