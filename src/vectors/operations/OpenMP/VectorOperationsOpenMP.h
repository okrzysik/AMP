#ifndef included_AMP_VectorOperationsOpenMP
#define included_AMP_VectorOperationsOpenMP


#include "AMP/vectors/operations/VectorOperations.h"


namespace AMP {
namespace LinearAlgebra {


/**
 * \brief  A default set of vector operations
 * \details VectorOperationsOpenMP impliments a default set of
 *    vector operations on the CPU.
 */
template<typename TYPE = double>
class VectorOperationsOpenMP : virtual public VectorOperations
{
public:
    // Constructor
    VectorOperationsOpenMP() {}

    //! Destructor
    virtual ~VectorOperationsOpenMP() {}

    //! Clone the operations
    std::shared_ptr<VectorOperations> cloneOperations() const override;

public:
//  function that operate on VectorData 
    void zero( VectorData &z ) override;
    void setToScalar( double alpha, VectorData &z ) override;
    void setRandomValues( VectorData &x ) override;    
    void setRandomValues( RNG::shared_ptr rng, VectorData &x ) override;    
    void copy( const VectorData &x, VectorData &z ) override;
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
    using VectorOperations::abs;
    using VectorOperations::add;
    using VectorOperations::addScalar;
    using VectorOperations::axpby;
    using VectorOperations::axpy;
    using VectorOperations::divide;
    using VectorOperations::dot;
    using VectorOperations::equals;
    using VectorOperations::linearSum;
    using VectorOperations::minQuotient;
    using VectorOperations::multiply;
    using VectorOperations::reciprocal;
    using VectorOperations::scale;
    using VectorOperations::setRandomValues;
    using VectorOperations::subtract;
    using VectorOperations::wrmsNorm;
    using VectorOperations::wrmsNormMask;
};


} // namespace LinearAlgebra
} // namespace AMP


#endif
