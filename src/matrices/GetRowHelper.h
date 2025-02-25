#ifndef included_AMP_Matrix_GetRowHelper
#define included_AMP_Matrix_GetRowHelper

#include "AMP/discretization/DOF_Manager.h"

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
                  std::shared_ptr<const AMP::Discretization::DOFManager> rightDOF )
        : d_leftDOF( leftDOF ),
          d_rightDOF( rightDOF ),
          d_beginCol( d_rightDOF->beginDOF() ),
          d_endCol( d_rightDOF->endDOF() )
    {
        AMP_ASSERT( d_leftDOF && d_rightDOF );
    }

    //! Destructor
    ~GetRowHelper() = default;

public:
    template<class BIGINT_TYPE, class INT_TYPE>
    void NNZ( BIGINT_TYPE row, INT_TYPE &N_local, INT_TYPE &N_remote );

    template<class BIGINT_TYPE>
    void getRow( BIGINT_TYPE row, BIGINT_TYPE *cols_local, BIGINT_TYPE *cols_remote );

private: // Member data
    const std::shared_ptr<const AMP::Discretization::DOFManager> d_leftDOF;
    const std::shared_ptr<const AMP::Discretization::DOFManager> d_rightDOF;
    const size_t d_beginCol;
    const size_t d_endCol;
    std::vector<size_t> d_rowDOFs;
};

} // namespace AMP::LinearAlgebra

#include "AMP/matrices/GetRowHelper.hpp"

#endif
