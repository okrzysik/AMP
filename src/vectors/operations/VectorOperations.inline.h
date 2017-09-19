#ifndef included_AMP_VectorOperations_inline
#define included_AMP_VectorOperations_inline


#include "vectors/data/VectorData.h"


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


/****************************************************************
* Access ghost values                                           *
****************************************************************/
inline bool VectorOperations::haGhosts() const
{
    if ( d_VectorData == nullptr )
        return false;
    return d_VectorData->d_Ghosts!=nullptr;
}
inline std::vector<double>& VectorOperations::getGhosts()
{
    return *(d_VectorData->d_Ghosts);
}


/****************************************************************
* Wrappers for shared_ptr                                       *
****************************************************************/
inline bool VectorOperations::equals( AMP::shared_ptr<const VectorOperations> x, double tol )
{
    return equals( *x, tol );
}
inline void VectorOperations::copy( AMP::shared_ptr<const VectorOperations> x )
{
    return copy( *x );
}
inline void VectorOperations::scale( double alpha, AMP::shared_ptr<const VectorOperations> x )
{
    return scale( alpha, *x );
}
inline void VectorOperations::add( AMP::shared_ptr<const VectorOperations> x, AMP::shared_ptr<const VectorOperations> y )
{
    return add( *x, *y );
}
void VectorOperations::addScalar( AMP::shared_ptr<const VectorOperations> x, double alpha )
{
    return addScalar( *x, alpha );
}
inline void VectorOperations::subtract( AMP::shared_ptr<const VectorOperations> x, AMP::shared_ptr<const VectorOperations> y )
{
    return subtract( *x, *y );
}
inline void VectorOperations::multiply( AMP::shared_ptr<const VectorOperations> x, AMP::shared_ptr<const VectorOperations> y )
{
    return multiply( *x, *y );
}
inline void VectorOperations::divide( AMP::shared_ptr<const VectorOperations> x, AMP::shared_ptr<const VectorOperations> y )
{
    return divide( *x, *y );
}
inline void VectorOperations::reciprocal( AMP::shared_ptr<const VectorOperations> x )
{
    return reciprocal( *x );
}
inline void VectorOperations::linearSum( double alpha,
    AMP::shared_ptr<const VectorOperations> x,
    double beta,
    AMP::shared_ptr<const VectorOperations> y )
{
    return linearSum( alpha, *x, beta, *y );
}
inline void VectorOperations::axpy( double alpha,
    AMP::shared_ptr<const VectorOperations> x,
    AMP::shared_ptr<const VectorOperations> y )
{
    return axpy( alpha, *x, *y );
}
inline void VectorOperations::axpby( double alpha, double beta, AMP::shared_ptr<const VectorOperations> x )
{
    return axpby( alpha, beta, *x );
}
inline void VectorOperations::abs( AMP::shared_ptr<const VectorOperations> x )
{
    return abs( *x );
}
inline double VectorOperations::dot( AMP::shared_ptr<const VectorOperations> x ) const
{
    return dot( *x );
}
inline double VectorOperations::minQuotient( AMP::shared_ptr<const VectorOperations> x, AMP::shared_ptr<const VectorOperations> y )
{
    return minQuotient( *x, *y );
}
inline double VectorOperations::wrmsNorm( 
    AMP::shared_ptr<const VectorOperations> x,
    AMP::shared_ptr<const VectorOperations> y )
{
    return wrmsNorm( *x, *y );
}
inline double VectorOperations::wrmsNormMask( AMP::shared_ptr<const VectorOperations> x,
   AMP::shared_ptr<const VectorOperations> y,
   AMP::shared_ptr<const VectorOperations> mask )
{
    return wrmsNormMask( *x, *y, *mask );
}


} // LinearAlgebra namespace
} // AMP namespace

#endif
