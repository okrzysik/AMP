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
class VectorOperationsDefault : public VectorOperations
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
    std::string VectorOpName() const override { return "VectorOperationsDefault"; }
    void copy( const VectorData &x, VectorData &z ) override;
    void zero( VectorData &z ) override;
    void setToScalar( const Scalar &alpha, VectorData &z ) override;
    void setRandomValues( VectorData &x ) override;
    void setRandomValues( std::shared_ptr<RNG> rng, VectorData &x ) override;
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
                      const Scalar &tol = 0.000001 ) const override;

public: // Pull VectorOperations into the current scope
    using VectorOperations::abs;
    using VectorOperations::add;
    using VectorOperations::addScalar;
    using VectorOperations::axpby;
    using VectorOperations::axpy;
    using VectorOperations::divide;
    using VectorOperations::dot;
    using VectorOperations::equals;
    using VectorOperations::L1Norm;
    using VectorOperations::L2Norm;
    using VectorOperations::linearSum;
    using VectorOperations::localL1Norm;
    using VectorOperations::localL2Norm;
    using VectorOperations::localMaxNorm;
    using VectorOperations::max;
    using VectorOperations::maxNorm;
    using VectorOperations::min;
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
