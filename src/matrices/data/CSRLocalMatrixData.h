#ifndef included_AMP_CSRLocalMatrixData_h
#define included_AMP_CSRLocalMatrixData_h

#include "AMP/matrices/data/MatrixData.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/memory.h"

#include <algorithm>
#include <functional>
#include <map>
#include <numeric>
#include <tuple>

namespace AMP::Discretization {
class DOFManager;
}

namespace AMP::LinearAlgebra {

// Forward declare CSRMatrixData to make it a friend
template<typename P, class A, class DIAG, class OFFD>
class CSRMatrixData;

template<typename Policy, class Allocator>
class CSRLocalMatrixData :
    public AMP::enable_shared_from_this<CSRLocalMatrixData<Policy, Allocator>>
{
public:
    template<typename P, class A, class DIAG, class OFFD>
    friend class CSRMatrixData;

    using gidx_t   = typename Policy::gidx_t;
    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;
    using gidxAllocator_t =
        typename std::allocator_traits<Allocator>::template rebind_alloc<gidx_t>;
    using lidxAllocator_t =
        typename std::allocator_traits<Allocator>::template rebind_alloc<lidx_t>;
    using scalarAllocator_t =
        typename std::allocator_traits<Allocator>::template rebind_alloc<scalar_t>;

    /** \brief Constructor
     * \param[in] outer Containing CSRMatrixData object
     * \param[in] params Description of the matrix
     * \param[in] is_diag True if this is the diag block, influences which dofs are used/ignored
     */
    explicit CSRLocalMatrixData( std::shared_ptr<MatrixParametersBase> params,
                                 AMP::Utilities::MemoryType memory_location,
                                 typename Policy::gidx_t first_row,
                                 typename Policy::gidx_t last_row,
                                 typename Policy::gidx_t first_col,
                                 typename Policy::gidx_t last_col,
                                 bool is_diag );

    //! Destructor
    virtual ~CSRLocalMatrixData();

    std::tuple<lidx_t *, gidx_t *, lidx_t *, scalar_t *> getDataFields()
    {
        return std::make_tuple(
            d_row_starts.get(), d_cols.get(), d_cols_loc.get(), d_coeffs.get() );
    }

    bool isDiag() const { return d_is_diag; }

    lidx_t numberOfNonZeros() const { return d_nnz; }

    lidx_t numLocalRows() const { return d_last_row - d_first_row; }

    lidx_t numLocalColumns() const { return d_last_col - d_first_col; }

    lidx_t numUniqueColumns() const { return d_ncols_unq; }

    lidx_t beginRow() const { return d_first_row; }

    lidx_t endRow() const { return d_last_row; }

    lidx_t beginCol() const { return d_first_col; }

    lidx_t endCol() const { return d_last_col; }

    void globalToLocalColumns();

    gidx_t *getColumnMap() const
    {
        if ( d_is_diag ) {
            return nullptr;
        }
        return d_cols_unq.get();
    }

    template<typename idx_t>
    void getColumnMap( std::vector<idx_t> &colMap ) const
    {
        // Don't do anything if offd and empty
        if ( !d_is_diag && d_is_empty ) {
            return;
        }

        AMP_INSIST( d_memory_location < AMP::Utilities::MemoryType::device,
                    "Copies from device to host memory not implemented yet" );

        colMap.resize( d_is_diag ? ( d_last_col - d_first_col ) : d_ncols_unq );

        if ( d_is_diag ) {
            std::iota( colMap.begin(), colMap.end(), d_first_col );
        } else {
            if constexpr ( std::is_same_v<idx_t, gidx_t> ) {
                std::copy( d_cols_unq.get(), d_cols_unq.get() + d_ncols_unq, colMap.begin() );
            } else {
                std::transform( d_cols_unq.get(),
                                d_cols_unq.get() + d_ncols_unq,
                                colMap.begin(),
                                []( gidx_t c ) -> idx_t { return c; } );
            }
        }
    }

    void setNNZ( const std::vector<lidx_t> &nnz );

    void printStats( bool show_zeros ) const
    {
        std::cout << ( d_is_diag ? "  diag block:" : "  offd block:" ) << std::endl;
        if ( d_is_empty ) {
            std::cout << "    EMPTY" << std::endl;
            return;
        }
        std::cout << "    first | last row: " << d_first_row << " | " << d_last_row << std::endl;
        std::cout << "    first | last col: " << d_first_col << " | " << d_last_col << std::endl;
        auto tot_nnz     = d_row_starts[d_num_rows];
        scalar_t avg_nnz = static_cast<scalar_t>( tot_nnz ) / static_cast<scalar_t>( d_num_rows );
        std::cout << "    avg nnz per row: " << avg_nnz << std::endl;
        std::cout << "    tot nnz: " << tot_nnz << std::endl;
        std::cout << "    row 0: ";
        for ( auto n = d_row_starts[0]; n < d_row_starts[1]; ++n ) {
            if ( d_coeffs[n] != 0 || show_zeros ) {
                std::cout << "(" << d_cols[n] << "," << d_coeffs[n] << "), ";
            }
        }
        std::cout << "    row last: ";
        for ( auto n = d_row_starts[d_num_rows - 1]; n < d_row_starts[d_num_rows]; ++n ) {
            if ( d_coeffs[n] != 0 || show_zeros ) {
                std::cout << "(" << d_cols[n] << "," << d_coeffs[n] << "), ";
            }
        }
        std::cout << std::endl << std::endl;
    }

protected:
    std::shared_ptr<CSRLocalMatrixData> cloneMatrixData();

    void getRowByGlobalID( const size_t local_row,
                           std::vector<size_t> &cols,
                           std::vector<double> &values ) const;

    void getValuesByGlobalID( const size_t local_row,
                              const size_t col,
                              void *values,
                              const typeID &id ) const;

    void addValuesByGlobalID( const size_t num_cols,
                              const size_t rows,
                              const size_t *cols,
                              const scalar_t *vals,
                              const typeID &id );

    void setValuesByGlobalID( const size_t num_cols,
                              const size_t rows,
                              const size_t *cols,
                              const scalar_t *vals,
                              const typeID &id );

    std::vector<size_t> getColumnIDs( const size_t local_row ) const;

    void sortColumns();

    // Data members passed from outer CSRMatrixData object
    const AMP::Utilities::MemoryType d_memory_location;
    const gidx_t d_first_row;
    const gidx_t d_last_row;
    const gidx_t d_first_col;
    const gidx_t d_last_col;

    const bool d_is_diag;
    bool d_is_empty = true;

    std::shared_ptr<lidx_t[]> d_row_starts;
    std::shared_ptr<gidx_t[]> d_cols;
    std::shared_ptr<gidx_t[]> d_cols_unq;
    std::shared_ptr<lidx_t[]> d_cols_loc;
    std::shared_ptr<scalar_t[]> d_coeffs;

    lidx_t d_num_rows  = 0;
    lidx_t d_nnz       = 0;
    lidx_t d_ncols_unq = 0;

    gidxAllocator_t d_gidxAllocator;
    lidxAllocator_t d_lidxAllocator;
    scalarAllocator_t d_scalarAllocator;

    std::shared_ptr<MatrixParametersBase> d_pParameters;
};

} // namespace AMP::LinearAlgebra

#endif
