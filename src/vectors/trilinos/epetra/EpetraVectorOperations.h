#ifndef included_AMP_EpetraVectorOperations
#define included_AMP_EpetraVectorOperations


#include "AMP/vectors/operations/VectorOperationsDefault.h"

#include <Epetra_Vector.h>


namespace AMP {
namespace LinearAlgebra {


/** \class EpetraVectorOperations
 * \brief A linear algebra engine that uses Epetra
 * \details  Use the Epetra implementation of the L1 BLAS routines.  Unlike other
 * libraries, it is very difficult to separate the data from the engine.  For this
 * reason, the EpetraVectorEngine contains the Epetra_Vector to operate on.
 */
class EpetraVectorOperations : virtual public VectorOperationsDefault<double>
{
public:
    // Constructor
    EpetraVectorOperations() {}

    // virtual void addScalar ( const VectorOperations & , double );
    virtual void setToScalar( double alpha ) override;
    virtual void scale( double alpha, const VectorOperations &x ) override;
    virtual void scale( double alpha ) override;
    virtual void add( const VectorOperations &x, const VectorOperations &y ) override;
    virtual void subtract( const VectorOperations &x, const VectorOperations &y ) override;
    virtual void multiply( const VectorOperations &x, const VectorOperations &y ) override;
    virtual void divide( const VectorOperations &x, const VectorOperations &y ) override;
    virtual void reciprocal( const VectorOperations &x ) override;
    virtual void linearSum( double alpha,
                            const VectorOperations &x,
                            double beta,
                            const VectorOperations &y ) override;
    virtual void
    axpy( double alpha, const VectorOperations &x, const VectorOperations &y ) override;
    virtual void axpby( double alpha, double beta, const VectorOperations &x ) override;
    virtual void abs( const VectorOperations &x ) override;
    virtual double min( void ) const override;
    virtual double max( void ) const override;
    virtual void setRandomValues( void ) override;
    virtual double L1Norm( void ) const override;
    virtual double L2Norm( void ) const override;
    virtual double maxNorm( void ) const override;
    virtual double dot( const VectorOperations &x ) const override;
    virtual double localMin( void ) const override;
    virtual double localMax( void ) const override;
    virtual double localL1Norm( void ) const override;
    virtual double localL2Norm( void ) const override;
    virtual double localMaxNorm() const override;
    virtual double localDot( const VectorOperations &x ) const override;


public: // Pull VectorOperations into the current scope
    using VectorOperationsDefault::abs;
    using VectorOperationsDefault::add;
    using VectorOperationsDefault::axpby;
    using VectorOperationsDefault::axpy;
    using VectorOperationsDefault::divide;
    using VectorOperationsDefault::dot;
    using VectorOperationsDefault::linearSum;
    using VectorOperationsDefault::minQuotient;
    using VectorOperationsDefault::multiply;
    using VectorOperationsDefault::reciprocal;
    using VectorOperationsDefault::scale;
    using VectorOperationsDefault::setRandomValues;
    using VectorOperationsDefault::subtract;
    using VectorOperationsDefault::wrmsNorm;
    using VectorOperationsDefault::wrmsNormMask;

protected:
    Epetra_Vector &getEpetra_Vector();
    const Epetra_Vector &getEpetra_Vector() const;
};


} // namespace LinearAlgebra
} // namespace AMP


#endif
