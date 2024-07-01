#ifndef included_AMP_CSRMatrixData_h
#define included_AMP_CSRMatrixData_h

#include "AMP/matrices/data/MatrixData.h"

#include <map>
#include <tuple>

namespace AMP::Discretization {
class DOFManager;
}

namespace AMP::LinearAlgebra {

template<typename Policy, class Allocator = std::allocator<int>>
class CSRMatrixData : public MatrixData
{
public:
    using gidx_t   = typename Policy::gidx_t;
    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;
    using gidxAllocator_t   = typename std::allocator_traits<Allocator>::template rebind_alloc<gidx_t>;
    using lidxAllocator_t   = typename std::allocator_traits<Allocator>::template rebind_alloc<lidx_t>;
    using scalarAllocator_t = typename std::allocator_traits<Allocator>::template rebind_alloc<scalar_t>;


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

    /** \brief  Get the global id of the beginning row
     * \return  beginning global row id
     */
    size_t beginRow() const override;

    /** \brief  Get the global id of the ending row
     * \return  ending global row id
     */
    size_t endRow() const override;

    size_t beginCol() const { return d_first_col; }

    std::tuple<lidx_t *, gidx_t *, lidx_t *, scalar_t *> getCSRDiagData()
    {
        return std::make_tuple( d_diag_matrix->d_nnz_per_row,
                                d_diag_matrix->d_cols,
                                d_diag_matrix->d_cols_loc,
                                d_diag_matrix->d_coeffs );
    }

    std::tuple<lidx_t *, gidx_t *, lidx_t *, scalar_t *> getCSROffDiagData()
    {
        return std::make_tuple( d_off_diag_matrix->d_nnz_per_row,
                                d_off_diag_matrix->d_cols,
                                d_off_diag_matrix->d_cols_loc,
                                d_off_diag_matrix->d_coeffs );
    }

    lidx_t* getDiagRowStarts() {
      return d_diag_matrix->d_row_starts;
    }

    lidx_t* getOffDiagRowStarts() {
      return d_off_diag_matrix->d_row_starts;
    }

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

    auto numberOfNonZerosOffDiag() const { return d_off_diag_matrix->d_nnz; }

    bool hasOffDiag() const { return !d_off_diag_matrix->d_is_empty; }

    auto getMemoryLocation() const { return d_memory_location; }

    template<typename idx_t>
    void getOffDiagColumnMap( std::vector<idx_t> &colMap ) const
    {
        // Don't do anything if empty
        if ( d_off_diag_matrix->d_is_empty ) { return; }

	// Column maps formed lazily, ensure it exists
	d_off_diag_matrix->findColumnMap();

	if ( d_memory_location < AMP::Utilities::MemoryType::device ) {
	
	  // Resize and fill colMap
	  colMap.resize( d_off_diag_matrix->d_ncols_unq );

	  if constexpr ( std::is_same_v<idx_t, gidx_t> ) {
	    std::copy( d_off_diag_matrix->d_cols_unq,
		       d_off_diag_matrix->d_cols_unq + d_off_diag_matrix->d_ncols_unq,
		       colMap.begin() );
	  } else {
	    std::transform( d_off_diag_matrix->d_cols_unq,
			    d_off_diag_matrix->d_cols_unq + d_off_diag_matrix->d_ncols_unq,
			    colMap.begin(),
			    []( gidx_t c ) -> idx_t { return c; } );
	  }
	} else {
	  AMP_ERROR( "Copies from device to host memory not implemented yet" );
	}
    }

private:
    // Private internal data class for managing the non-zero structure of the matrix
    // One instance will be made for the diagonal block and another for the off-diagonal block
    class CSRSerialMatrixData : public AMP::enable_shared_from_this<CSRSerialMatrixData>
    {
        // The outer CSRMatrixData class should have direct access to the internals of this class
        friend class CSRMatrixData<Policy>;

    public:
        /** \brief Constructor
         * \param[in] params Description of the matrix
         * \param[in] is_diag True if this is the diag block, influences which dofs are used/ignored
         */
        explicit CSRSerialMatrixData( const CSRMatrixData<Policy> &outer,
                                      std::shared_ptr<MatrixParametersBase> params,
                                      bool is_diag );

        explicit CSRSerialMatrixData( const CSRMatrixData<Policy> &outer );

        //! Destructor
        virtual ~CSRSerialMatrixData();

      std::shared_ptr<CSRSerialMatrixData> cloneMatrixData( const CSRMatrixData<Policy,Allocator> &outer );

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

    protected:
        const CSRMatrixData<Policy> &d_outer; // reference to the containing CSRMatrixData object
        bool d_is_diag  = true;
        bool d_is_empty = false;

        lidx_t *d_nnz_per_row = nullptr;
        lidx_t *d_row_starts  = nullptr;
        gidx_t *d_cols        = nullptr;
        gidx_t *d_cols_unq    = nullptr;
        lidx_t *d_cols_loc    = nullptr;
        scalar_t *d_coeffs    = nullptr;

        lidx_t d_num_rows  = 0;
        lidx_t d_nnz       = 0;
        lidx_t d_nnz_pad   = 0;
        lidx_t d_ncols_unq = 0;

        AMP::Utilities::MemoryType d_memory_location = AMP::Utilities::MemoryType::host;
        gidxAllocator_t gidxAllocator;
        lidxAllocator_t lidxAllocator;
        scalarAllocator_t scalarAllocator;

        std::shared_ptr<MatrixParametersBase> d_pParameters;

        bool d_own_data = true;
    };

protected:
    bool d_is_square   = true;
    gidx_t d_first_row = 0;
    gidx_t d_last_row  = 0;
    gidx_t d_first_col = 0;
    gidx_t d_last_col  = 0;
    lidx_t d_nnz       = 0;

    AMP::Utilities::MemoryType d_memory_location = AMP::Utilities::MemoryType::host;
    gidxAllocator_t gidxAllocator;
    lidxAllocator_t lidxAllocator;
    scalarAllocator_t scalarAllocator;


    std::shared_ptr<CSRSerialMatrixData> d_diag_matrix     = nullptr;
    std::shared_ptr<CSRSerialMatrixData> d_off_diag_matrix = nullptr;

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

template<typename Policy>
static CSRMatrixData<Policy> const *getCSRMatrixData( MatrixData const &A )
{
    auto ptr = dynamic_cast<CSRMatrixData<Policy> const *>( &A );
    AMP_INSIST( ptr, "dynamic cast from const MatrixData to const CSRMatrixData failed" );
    return ptr;
}

template<typename Policy>
static CSRMatrixData<Policy> *getCSRMatrixData( MatrixData &A )
{
    auto ptr = dynamic_cast<CSRMatrixData<Policy> *>( &A );
    AMP_INSIST( ptr, "dynamic cast from const MatrixData to const CSRMatrixData failed" );
    return ptr;
}

} // namespace AMP::LinearAlgebra

#endif
