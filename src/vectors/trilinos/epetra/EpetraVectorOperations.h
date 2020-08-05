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
    void axpy( double alpha, const VectorOperations &x, const VectorOperations &y ) override;
    void axpby( double alpha, double beta, const VectorOperations &x ) override;
    void abs( const VectorOperations &x ) override;
    double min( void ) const override;
    double max( void ) const override;
    void setRandomValues( void ) override;
    double L1Norm( void ) const override;
    double L2Norm( void ) const override;
    double maxNorm( void ) const override;
    double dot( const VectorOperations &x ) const override;

 protected:
    void setToScalar( double alpha, VectorData &z );
    void setRandomValues( VectorData &x );    
    void scale( double alpha, const VectorData &x, VectorData &y );
    void scale( double alpha, VectorData &x );
    void add( const VectorData &x, const VectorData &y, VectorData &z );
    void subtract( const VectorData &x, const VectorData &y, VectorData &z );
    void multiply( const VectorData &x, const VectorData &y, VectorData &z );
    void divide( const VectorData &x, const VectorData &y, VectorData &z );
    void reciprocal( const VectorData &x, VectorData &y );
    void linearSum( double alpha,
			   const VectorData &x,
			   double beta,
			   const VectorData &y,
			   VectorData &z);
    void axpy( double alpha, const VectorData &x, const VectorData &y, VectorData &z );
    void axpby( double alpha, double beta, const VectorData &x, VectorData &y );
    void abs( const VectorData &x, VectorData &z );

    double min( const VectorData &x ) const;
    double max( const VectorData &x ) const;
    double dot( const VectorData &x, const VectorData &y ) const;
    double L1Norm( const VectorData &x ) const;
    double L2Norm( const VectorData &x ) const;
    double maxNorm( const VectorData &x ) const;
    //    double localL1Norm( const VectorData &x ) const;
    //    double localL2Norm( const VectorData &x  ) const;
    //    double localMaxNorm( const VectorData &x ) const;

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

    using VectorOperationsDefault::localDot;
    using VectorOperationsDefault::localL1Norm;
    using VectorOperationsDefault::localL2Norm;
    using VectorOperationsDefault::localMax;
    using VectorOperationsDefault::localMaxNorm;
    using VectorOperationsDefault::localMin;

protected:
    Epetra_Vector &getEpetra_Vector();
    const Epetra_Vector &getEpetra_Vector() const;
};


} // namespace LinearAlgebra
} // namespace AMP


#endif
