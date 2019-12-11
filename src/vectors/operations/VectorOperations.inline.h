#ifndef included_AMP_VectorOperations_inline
#define included_AMP_VectorOperations_inline


#include "AMP/vectors/data/VectorData.h"


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
inline const AMP_MPI &VectorOperations::getComm() const
{
    return d_VectorData->getCommunicationList()->getComm();
}


/****************************************************************
 * Access ghost values                                           *
 ****************************************************************/
inline bool VectorOperations::hasGhosts() const
{
    if ( d_VectorData == nullptr )
        return false;
    return d_VectorData->d_Ghosts != nullptr;
}
inline std::vector<double> &VectorOperations::getGhosts() { return *( d_VectorData->d_Ghosts ); }


/****************************************************************
 * Wrappers for shared_ptr                                       *
 ****************************************************************/
inline bool VectorOperations::equals( std::shared_ptr<const VectorOperations> x, double tol )
{
    return equals( *x, tol );
}
inline void VectorOperations::copy( std::shared_ptr<const VectorOperations> x )
{
    return copy( *x );
}
inline void VectorOperations::scale( double alpha, std::shared_ptr<const VectorOperations> x )
{
    return scale( alpha, *x );
}
inline void VectorOperations::add( std::shared_ptr<const VectorOperations> x,
                                   std::shared_ptr<const VectorOperations> y )
{
    return add( *x, *y );
}
void VectorOperations::addScalar( std::shared_ptr<const VectorOperations> x, double alpha )
{
    return addScalar( *x, alpha );
}
inline void VectorOperations::subtract( std::shared_ptr<const VectorOperations> x,
                                        std::shared_ptr<const VectorOperations> y )
{
    return subtract( *x, *y );
}
inline void VectorOperations::multiply( std::shared_ptr<const VectorOperations> x,
                                        std::shared_ptr<const VectorOperations> y )
{
    return multiply( *x, *y );
}
inline void VectorOperations::divide( std::shared_ptr<const VectorOperations> x,
                                      std::shared_ptr<const VectorOperations> y )
{
    return divide( *x, *y );
}
inline void VectorOperations::reciprocal( std::shared_ptr<const VectorOperations> x )
{
    return reciprocal( *x );
}
inline void VectorOperations::linearSum( double alpha,
                                         std::shared_ptr<const VectorOperations> x,
                                         double beta,
                                         std::shared_ptr<const VectorOperations> y )
{
    return linearSum( alpha, *x, beta, *y );
}
inline void VectorOperations::axpy( double alpha,
                                    std::shared_ptr<const VectorOperations> x,
                                    std::shared_ptr<const VectorOperations> y )
{
    return axpy( alpha, *x, *y );
}
inline void
VectorOperations::axpby( double alpha, double beta, std::shared_ptr<const VectorOperations> x )
{
    return axpby( alpha, beta, *x );
}
inline void VectorOperations::abs( std::shared_ptr<const VectorOperations> x ) { return abs( *x ); }
inline double VectorOperations::dot( std::shared_ptr<const VectorOperations> x ) const
{
    return dot( *x );
}
inline double VectorOperations::minQuotient( std::shared_ptr<const VectorOperations> x,
                                             std::shared_ptr<const VectorOperations> y )
{
    return minQuotient( *x, *y );
}
inline double VectorOperations::wrmsNorm( std::shared_ptr<const VectorOperations> x,
                                          std::shared_ptr<const VectorOperations> y )
{
    return wrmsNorm( *x, *y );
}
inline double VectorOperations::wrmsNormMask( std::shared_ptr<const VectorOperations> x,
                                              std::shared_ptr<const VectorOperations> y,
                                              std::shared_ptr<const VectorOperations> mask )
{
    return wrmsNormMask( *x, *y, *mask );
}


} // namespace LinearAlgebra
} // namespace AMP

#endif
