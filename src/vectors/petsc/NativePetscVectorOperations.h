#ifndef included_AMP_NativePetscVectorOperations
#define included_AMP_NativePetscVectorOperations

#include "petsc.h"
#include "petsc/private/vecimpl.h"
#include "petscvec.h"

#include "AMP/vectors/operations/VectorOperationsDefault.h"
#include "AMP/vectors/petsc/NativePetscVector.h"

namespace AMP {
namespace LinearAlgebra {


/*! \struct Vec
    \brief PETSc vector
*/

/** \class NativePetscVector
 * \brief An AMP Vector that uses PETSc for parallel data management, linear algebra,
 * etc.
 * \details  This is an AMP wrapper to PETSc.  This is different from ManagedPetscVector
 * in that this class does not replace calls to Vec*.  Rather, it wraps these calls.
 * This class is used when PETSc is chosen as the default linear algebra engine.
 *
 * This class is not to be used directly, just through base class interfaces.
 * \see PetscVector
 * \see ManagedPetscVector
 */
class NativePetscVectorOperations : public VectorOperationsDefault<double>
{
private:
    void resetArray();
    void resetArray() const;

    //**********************************************************************
    // functions that operate on VectorData
    static Vec getPetscVec( VectorData &x );
    static Vec getPetscVec( const VectorData &x );
    static Vec getConstPetscVec( const VectorData &x );

    static NativePetscVector *getNativeVec( VectorData &vx );
    static const NativePetscVector *getNativeVec( const VectorData &vx );

public:
    void copy( const VectorData &x, VectorData &z ) override;
    void setToScalar( double alpha, VectorData &z ) override;
    void setRandomValues( VectorData &x ) override;
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

    double min( const VectorData &x ) const override;
    double max( const VectorData &x ) const override;
    double dot( const VectorData &x, const VectorData &y ) const override;
    double L1Norm( const VectorData &x ) const override;
    double L2Norm( const VectorData &x ) const override;
    double maxNorm( const VectorData &x ) const override;
    double localL1Norm( const VectorData &x ) const override;
    double localL2Norm( const VectorData &x ) const override;
    double localMaxNorm( const VectorData &x ) const override;

    void axpbypcz( double alpha,
                   const VectorData &x,
                   double beta,
                   const VectorData &y,
                   double gamma,
                   VectorData &z );
    //**********************************************************************

private:
public: // Pull VectorOperations into the current scope
    using VectorOperationsDefault::abs;
    using VectorOperationsDefault::add;
    using VectorOperationsDefault::axpby;
    using VectorOperationsDefault::axpy;
    using VectorOperationsDefault::divide;
    using VectorOperationsDefault::dot;
    using VectorOperationsDefault::L1Norm;
    using VectorOperationsDefault::L2Norm;
    using VectorOperationsDefault::linearSum;
    using VectorOperationsDefault::max;
    using VectorOperationsDefault::maxNorm;
    using VectorOperationsDefault::min;
    using VectorOperationsDefault::minQuotient;
    using VectorOperationsDefault::multiply;
    using VectorOperationsDefault::reciprocal;
    using VectorOperationsDefault::scale;
    using VectorOperationsDefault::setRandomValues;
    using VectorOperationsDefault::subtract;
    using VectorOperationsDefault::wrmsNorm;
    using VectorOperationsDefault::wrmsNormMask;
};


} // namespace LinearAlgebra
} // namespace AMP

#endif
