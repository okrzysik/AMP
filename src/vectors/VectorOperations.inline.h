#ifndef included_AMP_VectorOperations_inline
#define included_AMP_VectorOperations_inline


#include "vectors/VectorData.h"


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
* Get the comm                                                  *
****************************************************************/
inline bool VectorOperations::hasComm() const
{
    if ( d_VectorData == nullptr )
        return false;
    return d_VectorData->getCommunicationList() != nullptr;
}
inline const AMP_MPI& VectorOperations::getComm() const
{
    return d_VectorData->getCommunicationList()->getComm();
}


} // LinearAlgebra namespace
} // AMP namespace

#endif
