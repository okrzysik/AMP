#include "AMP/matrices/MatrixParameters.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Utilities.h"


namespace AMP::LinearAlgebra {


MatrixParameters::MatrixParameters( std::shared_ptr<AMP::Discretization::DOFManager> left,
                                    std::shared_ptr<AMP::Discretization::DOFManager> right,
                                    const AMP_MPI &comm )
    : d_comm( comm )

{
    AMP_ASSERT( left );
    AMP_ASSERT( right );
    AMP_ASSERT( !d_comm.isNull() );
    d_DOFManagerLeft  = left;
    d_DOFManagerRight = right;
}
} // namespace AMP::LinearAlgebra
