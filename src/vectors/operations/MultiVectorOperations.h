#ifndef included_AMP_MultiVectorOperations
#define included_AMP_MultiVectorOperations


#include "AMP/vectors/operations/VectorOperations.h"


namespace AMP {
namespace LinearAlgebra {

class MultiVectorData;

/**
 * \brief  A set of vector operations for multivectors
 * \details MultiVectorOperations impliments a default set of
 *    vector operations for multivectors.
 */
class MultiVectorOperations : public VectorOperations
{
public:
    // Constructor
    MultiVectorOperations() :VectorOperations() {}

    //! Destructor
    virtual ~MultiVectorOperations() {}

    //! Clone the operations
    std::shared_ptr<VectorOperations> cloneOperations() const override;

private:
    //  static function that operate on VectorData
    static VectorData *getVectorDataComponent( VectorData &x, size_t i );
    static const VectorData *getVectorDataComponent( const VectorData &x, size_t i );
    static MultiVectorData *getMultiVectorData( VectorData &x );
    static const MultiVectorData *getMultiVectorData( const VectorData &x );

public:
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
                    VectorData &z ) override;
    void axpy( double alpha, const VectorData &x, const VectorData &y, VectorData &z ) override;
    void axpby( double alpha, double beta, const VectorData &x, VectorData &y ) override;
    void abs( const VectorData &x, VectorData &z ) override;
    void addScalar( const VectorData &x, double alpha_in, VectorData &y ) override;

    double localMin( const VectorData &x ) const override;
    double localMax( const VectorData &x ) const override;
    double localL1Norm( const VectorData &x ) const override;
    double localL2Norm( const VectorData &x ) const override;
    double localMaxNorm( const VectorData &x ) const override;
    double localDot( const VectorData &x, const VectorData &y ) const override;
    double localMinQuotient( const VectorData &x, const VectorData &y ) const override;
    double localWrmsNorm( const VectorData &x, const VectorData &y ) const override;
    double localWrmsNormMask( const VectorData &x,
                              const VectorData &mask,
                              const VectorData &y ) const override;
    bool
    localEquals( const VectorData &x, const VectorData &y, double tol = 0.000001 ) const override;

    void updateVectorOperations( std::vector<std::shared_ptr<VectorOperations>> &ops );
protected:
    // Internal data
    std::vector<std::shared_ptr<VectorOperations>> d_operations;


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
