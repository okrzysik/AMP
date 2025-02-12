#ifndef included_AMP_Matrix_GetRowHelper
#define included_AMP_Matrix_GetRowHelper

#include "AMP/discretization/DOF_Manager.h"

#include <array>
#include <memory>
#include <vector>


namespace AMP::LinearAlgebra {


class GetRowHelper final
{
public:
    /** \brief Construct GetRowHelper
     * \details  This will construct the GetRowHelper from the left and right DOFManager
     * \param[in]  leftDOF      The left DOFManager
     * \param[in]  rightDOF     The right DOFManager
     */
    GetRowHelper( std::shared_ptr<const AMP::Discretization::DOFManager> leftDOF,
                  std::shared_ptr<const AMP::Discretization::DOFManager> rightDOF );


    //! Destructor
    ~GetRowHelper();


public: // Standard interfaces
    /** \brief  Get the number of non-zeros
     * \details  This will return the number of non-zeros for the row as [local,remote]
     * \param[in]  row          The row
     */
    std::array<size_t, 2> NNZ( size_t row ) const;


    /** \brief  Get the number of non-zeros
     * \details  This will return the number local of non-zeros
     * \param[in]  row          The row
     */
    size_t localNNZ( size_t row ) const;


    /** \brief  Get the number of non-zeros
     * \details  This will return the number remote of non-zeros
     * \param[in]  row          The row
     */
    size_t remoteNNZ( size_t row ) const;


    /** \brief  Get the row
     * \details  This will return the local and remote non-zero entries for the row
     * \param[in]  row          The row of interest
     * \param[out] local        The local non-zero entries (may be null)
     * \param[out] remote       The remote non-zero entries (may be null)
     */
    template<class BIGINT_TYPE>
    void getRow( BIGINT_TYPE row, BIGINT_TYPE *local, BIGINT_TYPE *remote ) const;

public: // Multi-row interfaces
    /** \brief  Get the number of non-zeros
     * \details  This will return the number remote of non-zeros for each row
     * \param[in]  start        The first row of interest
     * \param[in]  end          One past the last row of interest
     * \param[out] N_local        The number of local non-zero entries (may be null)
     * \param[out] N_remote       The number of remote non-zero entries (may be null)
     */
    template<class BIGINT_TYPE, class INT_TYPE>
    void NNZ( BIGINT_TYPE start, BIGINT_TYPE end, INT_TYPE *N_local, INT_TYPE *N_remote ) const;


    /** \brief  Get the number of non-zeros
     * \details  This will return the number remote of non-zeros for each row
     * \param[in]  N_rows       The number of rows of interest
     * \param[in]  rows         The rows of interest
     * \param[out] N_local        The number of local non-zero entries (may be null)
     * \param[out] N_remote       The number of remote non-zero entries (may be null)
     */
    template<class BIGINT_TYPE, class INT_TYPE>
    void NNZ( BIGINT_TYPE N_rows, BIGINT_TYPE *rows, INT_TYPE *N_local, INT_TYPE *N_remote ) const;


    /** \brief  Get the rows
     * \details  This will return the local and remote non-zero entries for the rows
     * \param[in]  start        The first row of interest
     * \param[in]  end          One past the last row of interest
     * \param[out] local        The local non-zero entries (may be null)
     * \param[out] remote       The remote non-zero entries (may be null)
     */
    template<class BIGINT_TYPE>
    void
    getRow( BIGINT_TYPE start, BIGINT_TYPE end, BIGINT_TYPE **local, BIGINT_TYPE **remote ) const;


    /** \brief  Get the rows
     * \details  This will return the local and remote non-zero entries for the rows
     * \param[in]  N_rows     The number of rows of interest
     * \param[in]  rows       The rows of interest
     * \param[out] local      The local non-zero entries (NOT null, but could be filled with nulls)
     * \param[out] remote     The remove non-zero entries (NOT null, but could be filled with nulls)
     */
    template<class BIGINT_TYPE>
    void getRow( BIGINT_TYPE N_rows,
                 BIGINT_TYPE *rows,
                 BIGINT_TYPE **local,
                 BIGINT_TYPE **remote ) const;


private: // Private routines
    void getRow2( size_t row, size_t *local, size_t *remote ) const;
    void reserve( size_t N );

private: // Member data
    std::shared_ptr<const AMP::Discretization::DOFManager> d_leftDOF;
    std::shared_ptr<const AMP::Discretization::DOFManager> d_rightDOF;
    std::array<size_t, 2> *d_NNZ = nullptr;
    size_t *d_local              = nullptr;
    size_t *d_remote             = nullptr;
    size_t *d_localOffset        = nullptr;
    size_t *d_remoteOffset       = nullptr;
    size_t d_size[2]             = { 0, 0 };
    size_t d_capacity[2]         = { 0, 0 };
};


} // namespace AMP::LinearAlgebra

#include "AMP/matrices/GetRowHelper.hpp"

#endif
