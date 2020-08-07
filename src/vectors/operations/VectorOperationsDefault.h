#ifndef included_AMP_VectorOperationsDefault
#define included_AMP_VectorOperationsDefault


#include "AMP/vectors/operations/VectorOperations.h"


namespace AMP {
namespace LinearAlgebra {


/**
 * \brief  A default set of vector operations
 * \details VectorOperationsDefault impliments a default set of
 *    vector operations on the CPU.
 */
template<typename TYPE = double>
class VectorOperationsDefault : virtual public VectorOperations
{
public:
    // Constructor
    VectorOperationsDefault() {}

    //! Destructor
    virtual ~VectorOperationsDefault() {}

    //! Clone the operations
    std::shared_ptr<VectorOperations> cloneOperations() const override;

public:
//  function that operate on VectorData 
    void copy( const VectorData &x, VectorData &z ) override;
    void zero( VectorData &z ) override;
    void setToScalar( double alpha, VectorData &z ) override;
    void setRandomValues( VectorData &x ) override;    
    void setRandomValues( RNG::shared_ptr rng, VectorData &x ) override;    
    void scale( double alpha, const VectorData &x, VectorData &y ) override;
    void scale( double alpha, VectorData &x ) override;
    void add( const VectorData &x, const VectorData &y, VectorData &z ) override;
    void subtract( const VectorData &x, const VectorData &y, VectorData &z ) override;
    void multiply( const VectorData &x, const VectorData &y, VectorData &z ) override;
    void divide( const VectorData &x, const VectorData &y, VectorData &z ) override;
    void reciprocal( const VectorData &x, VectorData &y ) override;
    void linearSum( double alpha,
			   const VectorData &x,
			   double beta,
			   const VectorData &y,
			   VectorData &z) override;
    void axpy( double alpha, const VectorData &x, const VectorData &y, VectorData &z ) override;
    void axpby( double alpha, double beta, const VectorData &x, VectorData &y ) override;
    void abs( const VectorData &x, VectorData &z ) override;
    void addScalar( const VectorData &x, double alpha_in, VectorData &y ) override;

    double localMin( const VectorData &x ) const override;
    double localMax( const VectorData &x ) const override;
    double localL1Norm( const VectorData &x ) const override;
    double localL2Norm( const VectorData &x  ) const override;
    double localMaxNorm( const VectorData &x ) const override;
    double localDot( const VectorData &x, const VectorData &y ) const override;
    double localMinQuotient( const VectorData &x, const VectorData &y ) const override;
    double localWrmsNorm( const VectorData &x, const VectorData &y ) const override;
    double localWrmsNormMask( const VectorData &x, const VectorData &mask, const VectorData &y ) const override;
    bool   localEquals( const VectorData &x, const VectorData &y, double tol = 0.000001 ) const override;

public: // Pull VectorOperations into the current scope
    using VectorOperations::setRandomValues;
    using VectorOperations::scale;
    using VectorOperations::add;
    using VectorOperations::subtract;
    using VectorOperations::multiply;
    using VectorOperations::divide;
    using VectorOperations::reciprocal;
    using VectorOperations::linearSum;
    using VectorOperations::axpy;
    using VectorOperations::axpby;
    using VectorOperations::abs;
    using VectorOperations::min;
    using VectorOperations::max;
    using VectorOperations::dot;
    using VectorOperations::L1Norm;
    using VectorOperations::L2Norm;
    using VectorOperations::maxNorm;
    using VectorOperations::localL1Norm;
    using VectorOperations::localL2Norm;
    using VectorOperations::localMaxNorm;
    using VectorOperations::addScalar;
    using VectorOperations::equals;
    using VectorOperations::minQuotient;
    using VectorOperations::wrmsNorm;
    using VectorOperations::wrmsNormMask;
};


} // namespace LinearAlgebra
} // namespace AMP


#endif
