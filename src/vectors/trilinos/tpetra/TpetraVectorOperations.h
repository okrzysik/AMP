#ifndef included_AMP_NativeTpetraVectorOperations_H_
#define included_AMP_NativeTpetraVectorOperations_H_

#include "AMP/vectors/Scalar.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/operations/VectorOperations.h"

#include <Teuchos_OrdinalTraits.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_Vector.hpp>

namespace AMP::LinearAlgebra {


/** \class TpetraVectorOperations
 * \brief An AMP Vector that uses Tpetra for parallel data management, linear algebra,
 * etc.
 * \details  This is an AMP wrapper to Tpetra.
 * This class is used when Tpetra is chosen as the default linear algebra engine.
 * This class is not to be used directly, just through base class interfaces.
 * \see TpetraVector
 */
template<typename ST = double,
         typename LO = int32_t,
         typename GO = int64_t,
         typename NT = Tpetra::Vector<>::node_type>
class TpetraVectorOperations : public VectorOperations
{
public:
    TpetraVectorOperations()          = default;
    virtual ~TpetraVectorOperations() = default;
    //  function that operate on VectorData
    void setToScalar( const Scalar &alpha, VectorData &z ) override;
    void setRandomValues( VectorData &x ) override;
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

    Scalar min( const VectorData &x ) const override;
    Scalar max( const VectorData &x ) const override;
    Scalar L1Norm( const VectorData &x ) const override;
    Scalar L2Norm( const VectorData &x ) const override;
    Scalar maxNorm( const VectorData &x ) const override;
    Scalar dot( const VectorData &x, const VectorData &y ) const override;
};


} // namespace AMP::LinearAlgebra

#endif
