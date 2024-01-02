#include "AMP/matrices/data/CSRMatrixData.hpp"
#include "AMP/matrices/data/hypre/HypreCSRPolicy.h"

namespace AMP::LinearAlgebra {
template class CSRMatrixData<HypreCSRPolicy>;
}
