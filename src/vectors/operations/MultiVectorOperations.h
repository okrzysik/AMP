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
    MultiVectorOperations() : VectorOperations() {}

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
    std::string VectorOpName() const override { return "MultiVectorOperations"; }
    void zero( VectorData &z ) override;
    void setToScalar( const Scalar &alpha, VectorData &z ) override;
    void setRandomValues( VectorData &x ) override;
    void setRandomValues( std::shared_ptr<RNG> rng, VectorData &x ) override;
    void copy( const VectorData &x, VectorData &z ) override;
    void scale( const Scalar &alpha, const VectorData &x, VectorData &y ) override;
    void scale( const Scalar &alpha, VectorData &x ) override;
    void add( const VectorData &x, const VectorData &y, VectorData &z ) override;
    void subtract( const VectorData &x, const VectorData &y, VectorData &z ) override;
    void multiply( const VectorData &x, const VectorData &y, VectorData &z ) override;
    void divide( const VectorData &x, const VectorData &y, VectorData &z ) override;
    void reciprocal( const VectorData &x, VectorData &y ) override;
    void linearSum( const Scalar &alpha,
                    const VectorData &x,
                    const Scalar &beta,
                    const VectorData &y,
                    VectorData &z ) override;
    void
    axpy( const Scalar &alpha, const VectorData &x, const VectorData &y, VectorData &z ) override;
    void
    axpby( const Scalar &alpha, const Scalar &beta, const VectorData &x, VectorData &y ) override;
    void abs( const VectorData &x, VectorData &z ) override;
    void addScalar( const VectorData &x, const Scalar &alpha_in, VectorData &y ) override;

    Scalar localMin( const VectorData &x ) const override;
    Scalar localMax( const VectorData &x ) const override;
    Scalar localL1Norm( const VectorData &x ) const override;
    Scalar localL2Norm( const VectorData &x ) const override;
    Scalar localMaxNorm( const VectorData &x ) const override;
    Scalar localDot( const VectorData &x, const VectorData &y ) const override;
    Scalar localMinQuotient( const VectorData &x, const VectorData &y ) const override;
    Scalar localWrmsNorm( const VectorData &x, const VectorData &y ) const override;
    Scalar localWrmsNormMask( const VectorData &x,
                              const VectorData &mask,
                              const VectorData &y ) const override;
    bool localEquals( const VectorData &x,
                      const VectorData &y,
                      const Scalar &tol = 1e-6 ) const override;

    void resetVectorOperations( std::vector<std::shared_ptr<VectorOperations>> &ops );

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
