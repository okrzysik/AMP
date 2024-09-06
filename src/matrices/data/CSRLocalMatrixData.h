#ifndef included_AMP_CSRLocalMatrixData_h
#define included_AMP_CSRLocalMatrixData_h

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

    void findColumnMap();

    // Data members passed from outer CSRMatrixData object
    const AMP::Utilities::MemoryType d_memory_location;
    const gidx_t d_first_row;
    const gidx_t d_last_row;
    const gidx_t d_first_col;
    const gidx_t d_last_col;

    const bool d_is_diag;
    bool d_is_empty = true;

    std::shared_ptr<lidx_t[]> d_nnz_per_row;
    std::shared_ptr<lidx_t[]> d_row_starts;
    std::shared_ptr<gidx_t[]> d_cols;
    std::shared_ptr<gidx_t[]> d_cols_unq;
    std::shared_ptr<lidx_t[]> d_cols_loc;
    std::shared_ptr<scalar_t[]> d_coeffs;

    lidx_t d_num_rows  = 0;
    lidx_t d_nnz       = 0;
    lidx_t d_nnz_pad   = 0;
    lidx_t d_ncols_unq = 0;

    gidxAllocator_t gidxAllocator;
    lidxAllocator_t lidxAllocator;
    scalarAllocator_t scalarAllocator;

    std::shared_ptr<MatrixParametersBase> d_pParameters;
};

} // namespace AMP::LinearAlgebra

#endif
