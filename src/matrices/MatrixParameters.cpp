#include "AMP/matrices/MatrixParameters.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Utilities.h"


namespace AMP {
namespace LinearAlgebra {


MatrixParameters::MatrixParameters( AMP::Discretization::DOFManager::shared_ptr left,
                                    AMP::Discretization::DOFManager::shared_ptr right,
                                    const AMP_MPI &comm )
    : d_comm( comm )

{
    AMP_ASSERT( left );
    AMP_ASSERT( right );
    d_DOFManagerLeft  = left;
    d_DOFManagerRight = right;
}
} // namespace LinearAlgebra
} // namespace AMP
