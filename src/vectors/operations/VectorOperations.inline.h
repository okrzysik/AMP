#ifndef included_AMP_VectorOperations_inline
#define included_AMP_VectorOperations_inline


#include "AMP/vectors/data/VectorData.h"


namespace AMP {
namespace LinearAlgebra {

/****************************************************************
 * Wrappers for shared_ptr                                       *
 ****************************************************************/

inline bool VectorOperations::equals( std::shared_ptr<const VectorData> x, std::shared_ptr<const VectorData> y, double tol ) const
{
  return equals( *x, *y, tol );
}

inline void VectorOperations::zero( std::shared_ptr<VectorData> x )
{
  zero(*x);
}
inline void VectorOperations::setToScalar( double alpha, std::shared_ptr<VectorData> x )
{
  setToScalar(alpha, *x);
}

inline void VectorOperations::setRandomValues( std::shared_ptr<VectorData> x )
{
  setRandomValues(*x);
}
 
inline void VectorOperations::copy( std::shared_ptr<const VectorData> x, std::shared_ptr<VectorData> y )
{
  copy( *x, *y );
}
inline void VectorOperations::scale( double alpha, std::shared_ptr<VectorData> x )
{
  scale( alpha, *x );
}
inline void VectorOperations::scale( double alpha, std::shared_ptr<const VectorData> x, std::shared_ptr<VectorData> y )
{
  scale( alpha, *x, *y );
}
inline void VectorOperations::add( std::shared_ptr<const VectorData> x,
                                   std::shared_ptr<const VectorData> y,
				   std::shared_ptr<VectorData> z )
{
  add( *x, *y, *z );
}
void VectorOperations::addScalar( std::shared_ptr<const VectorData> x, double alpha, std::shared_ptr<VectorData> y )
{
  addScalar( *x, alpha, *y );
}
inline void VectorOperations::subtract( std::shared_ptr<const VectorData> x,
                                        std::shared_ptr<const VectorData> y,
					std::shared_ptr<VectorData> z )
{
  subtract( *x, *y, *z );
}
inline void VectorOperations::multiply( std::shared_ptr<const VectorData> x,
                                        std::shared_ptr<const VectorData> y,
					std::shared_ptr<VectorData> z )
{
  multiply( *x, *y, *z );
}
inline void VectorOperations::divide( std::shared_ptr<const VectorData> x,
                                      std::shared_ptr<const VectorData> y,
				      std::shared_ptr<VectorData> z )
{
  divide( *x, *y, *z );
}
inline void VectorOperations::reciprocal( std::shared_ptr<const VectorData> x, std::shared_ptr<VectorData> y )
{
  reciprocal( *x, *y );
}
inline void VectorOperations::linearSum( double alpha,
                                         std::shared_ptr<const VectorData> x,
                                         double beta,
                                         std::shared_ptr<const VectorData> y,
					 std::shared_ptr<VectorData> z )
{
  linearSum( alpha, *x, beta, *y, *z );
}
inline void VectorOperations::axpy( double alpha,
                                    std::shared_ptr<const VectorData> x,
                                    std::shared_ptr<const VectorData> y,
				    std::shared_ptr<VectorData> z )
{
  return axpy( alpha, *x, *y, *z );
}
inline void
VectorOperations::axpby( double alpha, double beta, std::shared_ptr<const VectorData> x, std::shared_ptr<VectorData> y )
{
  axpby( alpha, beta, *x, *y );
}
 inline void VectorOperations::abs( std::shared_ptr<const VectorData> x, std::shared_ptr<VectorData> y ) { return abs( *x, *y ); }
inline double VectorOperations::dot( std::shared_ptr<const VectorData> x, std::shared_ptr<const VectorData> y ) const
{
    return dot( *x, *y );
}

inline double VectorOperations::minQuotient( std::shared_ptr<const VectorData> x,
                                             std::shared_ptr<const VectorData> y ) const
{
    return minQuotient( *x, *y );
}
inline double VectorOperations::wrmsNorm( std::shared_ptr<const VectorData> x,
                                          std::shared_ptr<const VectorData> y ) const
{
    return wrmsNorm( *x, *y );
}
inline double VectorOperations::wrmsNormMask( std::shared_ptr<const VectorData> x,
                                              std::shared_ptr<const VectorData> mask,
                                              std::shared_ptr<const VectorData> y ) const
{
  return wrmsNormMask( *x, *mask, *y );
}
 
} // namespace LinearAlgebra
} // namespace AMP

#endif
