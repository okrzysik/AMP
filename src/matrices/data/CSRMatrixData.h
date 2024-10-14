#ifndef included_AMP_CSRMatrixData_h
#define included_AMP_CSRMatrixData_h

#include "AMP/matrices/data/CSRLocalMatrixData.h"
#include "AMP/matrices/data/MatrixData.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/memory.h"

#include <functional>
#include <map>
#include <tuple>

namespace AMP::Discretization {
class DOFManager;
}

namespace AMP::LinearAlgebra {

template<typename Policy,
         class Allocator      = AMP::HostAllocator<int>,
         class DiagMatrixData = CSRLocalMatrixData<Policy, Allocator>,
         class OffdMatrixData = CSRLocalMatrixData<Policy, Allocator>>
class CSRMatrixData : public MatrixData
{
public:
    using gidx_t   = typename Policy::gidx_t;
    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;
    using gidxAllocator_t =
        typename std::allocator_traits<Allocator>::template rebind_alloc<gidx_t>;
    using lidxAllocator_t =
        typename std::allocator_traits<Allocator>::template rebind_alloc<lidx_t>;
    using scalarAllocator_t =
        typename std::allocator_traits<Allocator>::template rebind_alloc<scalar_t>;

    // Memory location, set by examining type of Allocator
    const AMP::Utilities::MemoryType d_memory_location;

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
     * \param[in] values  The values to add to the matrix
     * \param[in] id   typeID of raw data
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
     * \param[in] values  The values to set to the matrix
     * \param[in] id   typeID of raw data
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
     * \param[in] values  The values to get from the matrix (row-major ordering)
     * \param[in] id   typeID of raw data
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
    void makeConsistent( AMP::LinearAlgebra::ScatterType t ) override;

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

    /** \brief  Get the global id of the beginning row (inclusive)
     * \return  beginning global row id
     */
    size_t beginRow() const override;

    /** \brief  Get the global id of the ending row (exclusive)
     * \return  ending global row id
     */
    size_t endRow() const override;

    /** \brief  Get the global id of the beginning column (inclusive)
     * \return  beginning global row id
     */
    size_t beginCol() const { return d_first_col; }

    std::shared_ptr<DiagMatrixData> getDiagMatrix() { return d_diag_matrix; }

    std::shared_ptr<DiagMatrixData> getOffdMatrix() { return d_offd_matrix; }

    lidx_t *getDiagRowStarts() { return d_diag_matrix->d_row_starts.get(); }

    lidx_t *getOffDiagRowStarts() { return d_offd_matrix->d_row_starts.get(); }

    bool isSquare() const noexcept { return d_is_square; }

    std::shared_ptr<AMP::LinearAlgebra::Variable> getLeftVariable()
    {
        return d_pParameters->d_VariableLeft;
    }
    std::shared_ptr<AMP::LinearAlgebra::Variable> getRightVariable()
    {
        return d_pParameters->d_VariableRight;
    }

    auto numberOfNonZeros() const { return d_nnz; }

    auto numberOfNonZerosDiag() const { return d_diag_matrix->d_nnz; }

    auto numberOfNonZerosOffDiag() const { return d_offd_matrix->d_nnz; }

    bool hasOffDiag() const { return !d_offd_matrix->d_is_empty; }

    auto getMemoryLocation() const { return d_memory_location; }

    void sortColumns( MatrixSortScheme sort_type )
    {
        d_diag_matrix->sortColumns( sort_type );
        d_offd_matrix->sortColumns( sort_type );
    }

    template<typename idx_t>
    idx_t *getOffDiagColumnMap() const
    {
        return d_offd_matrix->getColumnMap();
    }

    template<typename idx_t>
    void getOffDiagColumnMap( std::vector<idx_t> &colMap ) const
    {
        d_offd_matrix->getColumnMap( colMap );
    }

protected:
    bool d_is_square   = true;
    gidx_t d_first_row = 0;
    gidx_t d_last_row  = 0;
    gidx_t d_first_col = 0;
    gidx_t d_last_col  = 0;
    lidx_t d_nnz       = 0;

    gidxAllocator_t d_gidxAllocator;
    lidxAllocator_t d_lidxAllocator;
    scalarAllocator_t d_scalarAllocator;

    std::shared_ptr<DiagMatrixData> d_diag_matrix;
    std::shared_ptr<OffdMatrixData> d_offd_matrix;

    std::shared_ptr<Discretization::DOFManager> d_leftDOFManager;
    std::shared_ptr<Discretization::DOFManager> d_rightDOFManager;

    //!  \f$A_{i,j}\f$ storage of off core matrix data
    std::map<gidx_t, std::map<gidx_t, scalar_t>> d_other_data;

    //!  \f$A_{i,j}\f$ storage of off core matrix data to set
    std::map<gidx_t, std::map<gidx_t, scalar_t>> d_ghost_data;

    //!  Update matrix data off-core
    void setOtherData( std::map<gidx_t, std::map<gidx_t, scalar_t>> &,
                       AMP::LinearAlgebra::ScatterType );
};

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
static CSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData> const *
getCSRMatrixData( MatrixData const &A )
{
    auto ptr =
        dynamic_cast<CSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData> const *>(
            &A );
    AMP_INSIST( ptr, "dynamic cast from const MatrixData to const CSRMatrixData failed" );
    return ptr;
}

template<typename Policy, class Allocator, class DiagMatrixData, class OffdMatrixData>
static CSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData> *
getCSRMatrixData( MatrixData &A )
{
    auto ptr =
        dynamic_cast<CSRMatrixData<Policy, Allocator, DiagMatrixData, OffdMatrixData> *>( &A );
    AMP_INSIST( ptr, "dynamic cast from const MatrixData to const CSRMatrixData failed" );
    return ptr;
}

} // namespace AMP::LinearAlgebra

#endif
