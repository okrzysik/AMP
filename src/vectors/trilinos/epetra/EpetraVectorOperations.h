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
    void setToScalar( double alpha ) override;
    void scale( double alpha, const VectorOperations &x ) override;
    void scale( double alpha ) override;
    void add( const VectorOperations &x, const VectorOperations &y ) override;
    void subtract( const VectorOperations &x, const VectorOperations &y ) override;
    void multiply( const VectorOperations &x, const VectorOperations &y ) override;
    void divide( const VectorOperations &x, const VectorOperations &y ) override;
    void reciprocal( const VectorOperations &x ) override;
    void linearSum( double alpha,
                            const VectorOperations &x,
                            double beta,
                            const VectorOperations &y ) override;
    void
    axpy( double alpha, const VectorOperations &x, const VectorOperations &y ) override;
    void axpby( double alpha, double beta, const VectorOperations &x ) override;
    void abs( const VectorOperations &x ) override;
    double min( void ) const override;
    double max( void ) const override;
    void setRandomValues( void ) override;
    double L1Norm( void ) const override;
    double L2Norm( void ) const override;
    double maxNorm( void ) const override;
    double dot( const VectorOperations &x ) const override;
    double localMin( void ) const override;
    double localMax( void ) const override;
    double localL1Norm( void ) const override;
    double localL2Norm( void ) const override;
    double localMaxNorm() const override;
    double localDot( const VectorOperations &x ) const override;


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
