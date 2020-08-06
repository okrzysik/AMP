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
#if 0
inline bool VectorOperations::equals( std::shared_ptr<const VectorData> x, std::shared_ptr<const VectorData> y, double tol )
{
  return equals( *x, *y, tol );
}
inline void VectorOperations::copy( std::shared_ptr<const VectorData> x, std::shared_ptr<const VectorData> y )
{
  return copy( *x, *y );
}
inline void VectorOperations::scale( double alpha, std::shared_ptr<const VectorData> x, std::shared_ptr<const VectorData> y )
{
  return scale( alpha, *x, *y );
}
inline void VectorOperations::add( std::shared_ptr<const VectorData> x,
                                   std::shared_ptr<const VectorData> y,
				   std::shared_ptr<const VectorData> z )
{
  return add( *x, *y, *z );
}
void VectorOperations::addScalar( std::shared_ptr<const VectorData> x, double alpha, std::shared_ptr<const VectorData> y )
{
  return addScalar( *x, alpha, *y );
}
inline void VectorOperations::subtract( std::shared_ptr<const VectorData> x,
                                        std::shared_ptr<const VectorData> y,
					std::shared_ptr<const VectorData> z )
{
  return subtract( *x, *y, *z );
}
inline void VectorOperations::multiply( std::shared_ptr<const VectorData> x,
                                        std::shared_ptr<const VectorData> y,
					std::shared_ptr<const VectorData> z )
{
    return multiply( *x, *y, *z );
}
inline void VectorOperations::divide( std::shared_ptr<const VectorData> x,
                                      std::shared_ptr<const VectorData> y,
				      std::shared_ptr<const VectorData> z )
{
    return divide( *x, *y, *z );
}
inline void VectorOperations::reciprocal( std::shared_ptr<const VectorData> x, std::shared_ptr<const VectorData> y )
{
    return reciprocal( *x, *y );
}
inline void VectorOperations::linearSum( double alpha,
                                         std::shared_ptr<const VectorData> x,
                                         double beta,
                                         std::shared_ptr<const VectorData> y,
					 std::shared_ptr<const VectorData> z )
{
    return linearSum( alpha, *x, beta, *y, *z );
}
inline void VectorOperations::axpy( double alpha,
                                    std::shared_ptr<const VectorData> x,
                                    std::shared_ptr<const VectorData> y,
				    std::shared_ptr<const VectorData> z )
{
    return axpy( alpha, *x, *y, *z );
}
inline void
VectorOperations::axpby( double alpha, double beta, std::shared_ptr<const VectorData> x, std::shared_ptr<const VectorData> y )
{
    return axpby( alpha, beta, *x, *y );
}
 inline void VectorOperations::abs( std::shared_ptr<const VectorData> x, std::shared_ptr<const VectorData> y ) { return abs( *x, *y ); }
inline double VectorOperations::dot( std::shared_ptr<const VectorData> x, std::shared_ptr<const VectorData> y ) const
{
    return dot( *x, *y );
}
inline double VectorOperations::minQuotient( std::shared_ptr<const VectorData> x,
                                             std::shared_ptr<const VectorData> y )
{
    return minQuotient( *x, *y );
}
inline double VectorOperations::wrmsNorm( std::shared_ptr<const VectorData> x,
                                          std::shared_ptr<const VectorData> y )
{
    return wrmsNorm( *x, *y );
}
inline double VectorOperations::wrmsNormMask( std::shared_ptr<const VectorData> x,
					      std::shared_ptr<const VectorData> mask,
					      std::shared_ptr<const VectorData> y )
{
  return wrmsNormMask( *x, *mask, *y );
}
#else
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
#endif
 

} // namespace LinearAlgebra
} // namespace AMP

#endif
