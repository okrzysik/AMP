#include "AMP/vectors/operations/ManagedVectorOperations.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/ManagedVector.h"
#include "AMP/vectors/data/ManagedVectorData.h"
#include "AMP/vectors/data/VectorData.h"

#include <iostream>
#include <stdexcept>
#include <string>
#include <typeinfo>


namespace AMP {
namespace LinearAlgebra {

static inline const ManagedVectorData *getManagedVector( const VectorData &x )
{
    auto y = dynamic_cast<const ManagedVectorData *>( &x );
    return y;
}
static inline ManagedVectorData *getManagedVector( VectorData &x )
{
    auto y = dynamic_cast<ManagedVectorData *>( &x );
    AMP_INSIST( y != nullptr, "x is not a ManagedVectorData" );
    return y;
}

//**********************************************************************
// Functions that operate on VectorData objects

void ManagedVectorOperations::copy( const VectorData &src, VectorData &dst )
{
    auto dst_managed = getManagedVector( dst );
    AMP_ASSERT( dst_managed );
    std::shared_ptr<Vector> vec1;
    std::shared_ptr<const Vector> vec2;
    auto src_managed = getManagedVector( src );
    if ( src_managed ) {
        // We are dealing with two managed vectors, check if they both have data engines
        if ( dst_managed->getVectorEngine() )
            vec1 = std::dynamic_pointer_cast<Vector>( dst_managed->getVectorEngine() );
        if ( src_managed->getVectorEngine() != nullptr )
            vec2 = std::dynamic_pointer_cast<const Vector>( src_managed->getVectorEngine() );
    }
    // Perform the copy
    if ( vec1 != nullptr && vec2 != nullptr ) {
        // We have two data engines, perform the copy between them
        vec1->copy( *vec2 );
    } else {
        // Default, general case
        VectorOperationsDefault::copy( src, dst );
    }
    dst_managed->dataChanged();
    dst_managed->makeConsistent( VectorData::ScatterType::CONSISTENT_SET );
}

void ManagedVectorOperations::setToScalar( const Scalar &alpha, VectorData &x )
{
    auto xm = getManagedVector( x );
    xm->getVectorEngine()->setToScalar( alpha );
    xm->dataChanged();
    xm->makeConsistent( VectorData::ScatterType::CONSISTENT_SET );
}

void ManagedVectorOperations::setRandomValues( VectorData &x )
{
    auto xm = getManagedVector( x );
    xm->getVectorEngine()->setRandomValues();
    xm->dataChanged();
    xm->makeConsistent( VectorData::ScatterType::CONSISTENT_SET );
}

void ManagedVectorOperations::scale( const Scalar &alpha, const VectorData &x, VectorData &y )
{
    auto x2 = getManagedVector( x );
    auto y2 = getManagedVector( y );
    if ( x2 != nullptr ) {
        y2->getVectorEngine()->scale( alpha, *x2->getVectorEngine() );
    } else {
        VectorOperationsDefault::scale( alpha, x, y );
    }
    y2->dataChanged();
}

void ManagedVectorOperations::scale( const Scalar &alpha, VectorData &x )
{
    auto y = getManagedVector( x );
    y->getVectorEngine()->scale( alpha );
    y->dataChanged();
}

void ManagedVectorOperations::add( const VectorData &x, const VectorData &y, VectorData &z )
{
    auto x2 = getManagedVector( x );
    auto y2 = getManagedVector( y );
    auto z2 = getManagedVector( z );
    if ( x2 != nullptr && y2 != nullptr ) {
        z2->getVectorEngine()->add( *x2->getVectorEngine(), *y2->getVectorEngine() );
    } else {
        VectorOperationsDefault::add( x, y, z );
    }
    z2->dataChanged();
}

void ManagedVectorOperations::subtract( const VectorData &x, const VectorData &y, VectorData &z )
{
    auto x2 = getManagedVector( x );
    auto y2 = getManagedVector( y );
    auto z2 = getManagedVector( z );
    if ( x2 != nullptr && y2 != nullptr ) {
        z2->getVectorEngine()->subtract( *x2->getVectorEngine(), *y2->getVectorEngine() );
    } else {
        VectorOperationsDefault::subtract( x, y, z );
    }
    z2->dataChanged();
}

void ManagedVectorOperations::multiply( const VectorData &x, const VectorData &y, VectorData &z )
{
    auto x2 = getManagedVector( x );
    auto y2 = getManagedVector( y );
    auto z2 = getManagedVector( z );
    if ( x2 != nullptr && y2 != nullptr ) {
        z2->getVectorEngine()->multiply( *x2->getVectorEngine(), *y2->getVectorEngine() );
    } else {
        VectorOperationsDefault::multiply( x, y, z );
    }
    z2->dataChanged();
}

void ManagedVectorOperations::divide( const VectorData &x, const VectorData &y, VectorData &z )
{
    auto x2 = getManagedVector( x );
    auto y2 = getManagedVector( y );
    auto z2 = getManagedVector( z );
    if ( x2 != nullptr && y2 != nullptr ) {
        z2->getVectorEngine()->divide( *x2->getVectorEngine(), *y2->getVectorEngine() );
    } else {
        VectorOperationsDefault::divide( x, y, z );
    }
    z2->dataChanged();
}

void ManagedVectorOperations::reciprocal( const VectorData &x, VectorData &y )
{
    auto x2 = getManagedVector( x );
    auto y2 = getManagedVector( y );
    if ( x2 != nullptr ) {
        y2->getVectorEngine()->reciprocal( *x2->getVectorEngine() );
    } else {
        VectorOperationsDefault::reciprocal( x, y );
    }
    y2->dataChanged();
}

void ManagedVectorOperations::linearSum( const Scalar &alpha,
                                         const VectorData &x,
                                         const Scalar &beta,
                                         const VectorData &y,
                                         VectorData &z )
{
    auto x2 = getManagedVector( x );
    auto y2 = getManagedVector( y );
    auto z2 = getManagedVector( z );
    if ( x2 != nullptr && y2 != nullptr ) {
        z2->getVectorEngine()->linearSum(
            alpha, *x2->getVectorEngine(), beta, *y2->getVectorEngine() );
    } else {
        VectorOperationsDefault::linearSum( alpha, x, beta, y, z );
    }
    z2->dataChanged();
}

void ManagedVectorOperations::axpy( const Scalar &alpha,
                                    const VectorData &x,
                                    const VectorData &y,
                                    VectorData &z )
{
    linearSum( alpha, x, 1.0, y, z );
}

void ManagedVectorOperations::axpby( const Scalar &alpha,
                                     const Scalar &beta,
                                     const VectorData &x,
                                     VectorData &z )
{
    linearSum( alpha, x, beta, z, z );
}

void ManagedVectorOperations::abs( const VectorData &x, VectorData &y )
{
    auto x2 = getManagedVector( x );
    auto y2 = getManagedVector( y );
    if ( x2 != nullptr ) {
        y2->getVectorEngine()->abs( *x2->getVectorEngine() );
    } else {
        VectorOperationsDefault::abs( x, y );
    }
    y2->dataChanged();
}

Scalar ManagedVectorOperations::min( const VectorData &x ) const
{
    auto x2 = getManagedVector( x );
    return x2->getVectorEngine()->min();
}

Scalar ManagedVectorOperations::max( const VectorData &x ) const
{
    auto x2 = getManagedVector( x );
    return x2->getVectorEngine()->max();
}

Scalar ManagedVectorOperations::dot( const VectorData &x, const VectorData &y ) const
{
    auto x2 = getManagedVector( x );
    auto y2 = getManagedVector( y );
    if ( ( x2 != nullptr ) && ( y2 != nullptr ) ) {
        return y2->getVectorEngine()->dot( *x2->getVectorEngine() );
    }
    return VectorOperationsDefault::dot( x, y );
}

Scalar ManagedVectorOperations::L1Norm( const VectorData &x ) const
{
    auto x2 = getManagedVector( x );
    return x2->getVectorEngine()->L1Norm();
}

Scalar ManagedVectorOperations::L2Norm( const VectorData &x ) const
{
    auto x2 = getManagedVector( x );
    return x2->getVectorEngine()->L2Norm();
}

Scalar ManagedVectorOperations::maxNorm( const VectorData &x ) const
{
    auto x2 = getManagedVector( x );
    return x2->getVectorEngine()->maxNorm();
}

} // namespace LinearAlgebra
} // namespace AMP
