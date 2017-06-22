#include "utils/Utilities.h"
#include "vectors/VectorEngine.h"
#include "vectors/petsc/ManagedPetscVector.h"
#include "vectors/petsc/PetscVector.h"
#ifdef USE_EXT_TRILINOS
#include "vectors/trilinos/epetra/EpetraVectorEngine.h"
#endif


#include "petscsys.h"


// Macro to cast the petsc vector y to an AMP ManagedPetscVector x
#define PETSC_RECAST( x, y )                    \
    AMP::LinearAlgebra::ManagedPetscVector *x = \
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( y->data )


// Overridden Petsc functions!
PetscErrorCode
_AMP_setvaluesblocked( Vec, PetscInt, const PetscInt[], const PetscScalar[], InsertMode )
{
    AMP_ERROR( "31 Not implemented" );
    return 0;
}
PetscErrorCode _AMP_view( Vec, PetscViewer )
{
    AMP_ERROR( "32 Not implemented" );
    return 0;
}
PetscErrorCode _AMP_placearray( Vec, const PetscScalar * )
{
    AMP_ERROR( "33 Not implemented" );
    return 0;
} /* place data array */
PetscErrorCode _AMP_replacearray( Vec, const PetscScalar * )
{
    AMP_ERROR( "34 Not implemented" );
    return 0;
} /* replace data array */
PetscErrorCode _AMP_viewnative( Vec, PetscViewer )
{
    AMP_ERROR( "42 Not implemented" );
    return 0;
}
PetscErrorCode _AMP_setlocaltoglobalmapping( Vec, ISLocalToGlobalMapping )
{
    AMP_ERROR( "44 Not implemented" );
    return 0;
}
PetscErrorCode
_AMP_setvalueslocal( Vec, PetscInt, const PetscInt *, const PetscScalar *, InsertMode )
{
    AMP_ERROR( "45 Not implemented" );
    return 0;
}

PetscErrorCode _AMP_getvalues( Vec, PetscInt, const PetscInt[], PetscScalar[] )
{
    AMP_ERROR( "52 Not implemented" );
    return 0;
}


PetscErrorCode _AMP_assemblybegin( Vec ) { return 0; }
PetscErrorCode _AMP_assemblyend( Vec ) { return 0; }

#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
PetscErrorCode _AMP_setoption( Vec, VecOption, PetscTruth ) { return 0; }
PetscErrorCode _AMP_load( PetscViewer, const VecType, Vec * )
{
    AMP_ERROR( "48 Not implemented" );
    return 0;
}
PetscErrorCode _AMP_loadintovector( PetscViewer, Vec )
{
    AMP_ERROR( "40 Not implemented" );
    return 0;
}
PetscErrorCode _AMP_loadintovectornative( PetscViewer, Vec )
{
    AMP_ERROR( "41 Not implemented" );
    return 0;
}
#elif PETSC_VERSION_GE(3,2,0)
PetscErrorCode _AMP_setoption( Vec, VecOption, PetscBool ) { return 0; }
PetscErrorCode _AMP_load( Vec, PetscViewer )
{
    AMP_ERROR( "48 Not implemented" );
    return 0;
}
#else
#error Not programmed for this version of petsc
#endif


// Inserts or adds values into certain locations of a vector.
PetscErrorCode
_AMP_setvalues( Vec px, PetscInt ni, const PetscInt ix[], const PetscScalar y[], InsertMode iora )
{
    PETSC_RECAST( x, px );
    auto indices = new size_t[ni];
    auto vals    = new double[ni];
    for ( PetscInt i = 0; i < ni; i++ ) {
        indices[i] = static_cast<size_t>( ix[i] );
        vals[i]    = static_cast<double>( y[i] );
    }
    if ( iora == INSERT_VALUES ) {
        x->setValuesByGlobalID( ni, indices, vals );
    } else if ( iora == ADD_VALUES ) {
        x->addValuesByGlobalID( ni, indices, vals );
    } else {
        AMP_ERROR( "Invalid option for InsertMode" );
    }
    delete[] indices;
    delete[] vals;
    return 0;
}

#if PETSC_VERSION_(3,7,5)
PetscErrorCode _AMP_shift( Vec px, PetscScalar s)
{
    PETSC_RECAST( x, px );

    auto cur = x->VectorData::begin();
    auto end = x->VectorData::end();
    while ( cur != end ) {
        *cur = s+ ( *cur );
        cur++;
    }
    return 0;
}
#elif PETSC_VERSION_LT(3,7,5)
// This function makes no sense wrt the PETSc interface VecShift( Vec, PetscScalar );
PetscErrorCode _AMP_shift( Vec )
{
    AMP_ERROR( "This function cannot be implemented as designed" );
    return 0;
}
#else
#error Not programmed for this version of petsc
#endif

PetscErrorCode
_AMP_axpbypcz( Vec c, PetscScalar alpha, PetscScalar beta, PetscScalar gamma, Vec a, Vec b )
{
    PETSC_RECAST( x, a );
    PETSC_RECAST( y, b );
    PETSC_RECAST( z, c );

    if ( z->isAnAliasOf( *x ) ) {
        z->linearSum( alpha + gamma, *x, beta, *y );
    } else if ( z->isAnAliasOf( *y ) ) {
        z->linearSum( alpha, *x, beta + gamma, *y );
    } else {
        z->linearSum( alpha, *x, gamma, *z );
        z->linearSum( beta, *y, 1., *z );
    }

    PetscObjectStateIncrease( reinterpret_cast<::PetscObject>( c ) );
    return 0;
}

PetscErrorCode _AMP_max( Vec a, PetscInt *p, PetscReal *ans )
{
    if ( ( p != nullptr ) && ( p != PETSC_NULL ) ) {
        AMP_ERROR( "Cannot find position for max" );
    }
    PETSC_RECAST( x, a );
    *ans = x->max();
    return 0;
}

PetscErrorCode _AMP_min( Vec a, PetscInt *p, PetscReal *ans )
{
    if ( ( p != nullptr ) && ( p != PETSC_NULL ) ) {
        AMP_ERROR( "Cannot find position for max" );
    }
    PETSC_RECAST( x, a );
    *ans = x->min();
    return 0;
}

PetscErrorCode _AMP_aypx( Vec b, PetscScalar alpha, Vec a )
{
    PETSC_RECAST( x, a );
    PETSC_RECAST( y, b );
    y->linearSum( alpha, *y, 1., *x );
    return 0;
} /* y = x + alpha * y */

PetscErrorCode _AMP_dot_local( Vec a, Vec b, PetscScalar *ans )
{
    PETSC_RECAST( x, a );
    PETSC_RECAST( y, b );
    *ans = x->localDot( *y );
    return 0;
}

PetscErrorCode _AMP_tdot_local( Vec a, Vec b, PetscScalar *ans )
{
    return _AMP_dot_local( a, b, ans );
}

PetscErrorCode _AMP_mdot_local( Vec a, PetscInt num, const Vec array[], PetscScalar *ans )
{
    for ( PetscInt i = 0; i != num; i++ )
        _AMP_dot_local( a, array[i], ans + i );
    return 0;
}

PetscErrorCode _AMP_mtdot_local( Vec a, PetscInt num, const Vec array[], PetscScalar *ans )
{
    for ( PetscInt i = 0; i != num; i++ )
        _AMP_dot_local( a, array[i], ans + i );
    return 0;
}

PetscErrorCode _AMP_exp( Vec a )
{
    PETSC_RECAST( x, a );

    auto cur = x->VectorData::begin();
    auto end = x->VectorData::end();
    while ( cur != end ) {
        *cur = exp( *cur );
        cur++;
    }
    return 0;
}

PetscErrorCode _AMP_log( Vec a )
{
    PETSC_RECAST( x, a );

    auto cur = x->VectorData::begin();
    auto end = x->VectorData::end();
    while ( cur != end ) {
        *cur = log( *cur );
        cur++;
    }
    return 0;
}

PetscErrorCode _AMP_pointwisemin( Vec a, Vec b, Vec c )
{
    PETSC_RECAST( x, a );
    PETSC_RECAST( y, b );
    PETSC_RECAST( z, c );

    AMP_INSIST( x->getLocalSize() == y->getLocalSize(), "Incompatible vectors" );
    AMP_INSIST( x->getLocalSize() == z->getLocalSize(), "Incompatible vectors" );

    auto xi = x->VectorData::begin();
    auto xe = x->VectorData::end();
    auto yi = y->VectorData::constBegin();
    auto zi = z->VectorData::constBegin();
    while ( xi != xe ) {
        *xi = std::min( *yi, *zi );
        xi++;
        yi++;
        zi++;
    }
    return 0;
}

PetscErrorCode _AMP_pointwisemax( Vec a, Vec b, Vec c )
{
    PETSC_RECAST( x, a );
    PETSC_RECAST( y, b );
    PETSC_RECAST( z, c );

    AMP_INSIST( x->getLocalSize() == y->getLocalSize(), "Incompatible vectors" );
    AMP_INSIST( x->getLocalSize() == z->getLocalSize(), "Incompatible vectors" );

    auto xi = x->VectorData::begin();
    auto xe = x->VectorData::end();
    auto yi = y->VectorData::constBegin();
    auto zi = z->VectorData::constBegin();
    while ( xi != xe ) {
        *xi = std::max( *yi, *zi );
        xi++;
        yi++;
        zi++;
    }
    return 0;
}

PetscErrorCode _AMP_pointwisemaxabs( Vec a, Vec b, Vec c )
{
    PETSC_RECAST( x, a );
    PETSC_RECAST( y, b );
    PETSC_RECAST( z, c );

    AMP_INSIST( x->getLocalSize() == y->getLocalSize(), "Incompatible vectors" );
    AMP_INSIST( x->getLocalSize() == z->getLocalSize(), "Incompatible vectors" );

    auto xi = x->VectorData::begin();
    auto yi = y->VectorData::begin();
    auto zi = z->VectorData::begin();
    auto xe = x->VectorData::end();
    while ( xi != xe ) {
        *xi = std::max( fabs( *yi ), fabs( *zi ) );
        xi++;
        yi++;
        zi++;
    }
    return 0;
}

PetscErrorCode _AMP_pointwisemult( Vec a, Vec b, Vec c )
{
    PETSC_RECAST( x, a );
    PETSC_RECAST( y, b );
    PETSC_RECAST( z, c );
    x->multiply( *y, *z );
    return 0;
}

PetscErrorCode _AMP_pointwisedivide( Vec a, Vec b, Vec c )
{
    PETSC_RECAST( x, a );
    PETSC_RECAST( y, b );
    PETSC_RECAST( z, c );
    x->divide( *y, *z );
    return 0;
}

PetscErrorCode _AMP_sqrt( Vec a )
{
    PETSC_RECAST( x, a );
    auto cur = x->VectorData::begin();
    auto end = x->VectorData::end();
    while ( cur != end ) {
        *cur = sqrt( fabs( *cur ) );
        cur++;
    }
    return 0;
}

PetscErrorCode _AMP_setrandom( Vec a, PetscRandom )
{
    PETSC_RECAST( x, a );
    x->setRandomValues();
    return 0;
} /* set y[j] = random numbers */

PetscErrorCode _AMP_conjugate( Vec )
{
    return 0; // Not dealing with complex right now
}

PetscErrorCode _AMP_axpby( Vec b, PetscScalar alpha, PetscScalar beta, Vec a )
{
    PETSC_RECAST( x, a );
    PETSC_RECAST( y, b );
    y->axpby( alpha, beta, *x );
    return 0;
}

PetscErrorCode _AMP_swap( Vec a, Vec b )
{
    PETSC_RECAST( x, a );
    PETSC_RECAST( y, b );
    x->swapVectors( *y );
    return 0;
}

PetscErrorCode _AMP_getsize( Vec a, PetscInt *ans )
{
    PETSC_RECAST( x, a );
    size_t size = x->getGlobalSize();
    if ( sizeof( PetscInt ) < 8 ) {
        AMP_ASSERT( size < 0x80000000 );
    }
    *ans = (PetscInt) size;
    return 0;
}

PetscErrorCode _AMP_maxpointwisedivide( Vec a, Vec b, PetscReal *res )
{
    PETSC_RECAST( x, a );
    PETSC_RECAST( y, b );

    auto cur_x = x->VectorData::constBegin();
    auto cur_y = y->VectorData::constBegin();
    auto end_x = x->VectorData::constEnd();
    double local_res = 0.0;
    while ( cur_x != end_x ) {
        if ( *cur_y == 0.0 ) {
            local_res = std::max( local_res, fabs( *cur_x ) );
        } else {
            local_res = std::max( local_res, fabs( ( *cur_x ) / ( *cur_y ) ) );
        }
        cur_x++;
        cur_y++;
    }

    *res = x->getComm().maxReduce( local_res );

    return 0;
}

PetscErrorCode _AMP_scale( Vec a, PetscScalar alpha )
{
    PETSC_RECAST( x, a );
    x->scale( alpha );
    return 0;
}

PetscErrorCode _AMP_copy( Vec in, Vec out )
{
    PETSC_RECAST( x, in );
    PETSC_RECAST( y, out );
    y->copyVector( x->shared_from_this() );
    return 0;
}

PetscErrorCode _AMP_maxpy( Vec v, PetscInt num, const PetscScalar *alpha, Vec *vecs )
{
    for ( int i = 0; i != num; i++ )
        VecAXPY( v, alpha[i], vecs[i] );
    return 0;
}

PetscErrorCode _AMP_dot( Vec a, Vec b, PetscScalar *ans )
{
    PETSC_RECAST( x, a );
    PETSC_RECAST( y, b );
    *ans = x->dot( *y );
    return 0;
}

PetscErrorCode _AMP_mdot( Vec v, PetscInt num, const Vec vec[], PetscScalar *ans )
{
    for ( PetscInt i = 0; i != num; i++ )
        VecDot( v, vec[i], ans + i );
    return 0;
}

PetscErrorCode _AMP_tdot( Vec a, Vec b, PetscScalar *ans )
{
    PETSC_RECAST( x, a );
    PETSC_RECAST( y, b );
    *ans = x->dot( *y );
    return 0;
}

PetscErrorCode _AMP_mtdot( Vec v, PetscInt num, const Vec vec[], PetscScalar *ans )
{
    for ( PetscInt i = 0; i != num; i++ )
        VecTDot( v, vec[i], ans + i );
    return 0;
}


PetscErrorCode _AMP_destroyvecs( PetscInt num, Vec vecArray[] )
{
    for ( PetscInt i = 0; i != num; i++ )
        PETSC::vecDestroy( &vecArray[i] );
    delete[] vecArray;
    return 0;
}

PetscErrorCode _AMP_axpy( Vec out, PetscScalar alpha, Vec in )
{
    AMP::LinearAlgebra::ManagedPetscVector *x =
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( in->data );
    AMP::LinearAlgebra::ManagedPetscVector *y =
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( out->data );
    y->axpy( alpha, *x, *y );
    return 0;
}


PetscErrorCode _AMP_waxpy( Vec w, PetscScalar alpha, Vec x, Vec y )
{
    AMP::LinearAlgebra::ManagedPetscVector *xIn =
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( x->data );
    AMP::LinearAlgebra::ManagedPetscVector *yIn =
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( y->data );
    AMP::LinearAlgebra::ManagedPetscVector *wOut =
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( w->data );

    AMP_INSIST( ( wOut != xIn ) && ( wOut != yIn ),
                "ERROR: _AMP_waxpy: w cannot be the same as x or y" );

    wOut->axpy( alpha, *xIn, *yIn );
    return 0;
}


PetscErrorCode _AMP_norm_local( Vec in, NormType type, PetscReal *ans )
{
    if ( type == NORM_1 )
        *ans =
            reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( in->data )->localL1Norm();
    else if ( type == NORM_2 )
        *ans =
            reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( in->data )->localL2Norm();
    else if ( type == NORM_INFINITY )
        *ans =
            reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( in->data )->localMaxNorm();
    else if ( type == NORM_1_AND_2 ) {
        *ans =
            reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( in->data )->localL1Norm();
        *( ans + 1 ) =
            reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( in->data )->localL2Norm();
    } else
        AMP_ERROR( "Unknown norm type" );
    if ( type != NORM_1_AND_2 ) {
        PetscObjectComposedDataSetReal(
            reinterpret_cast<::PetscObject>( in ), NormIds[type], ans[0] );
    }
    PetscObjectStateIncrease( reinterpret_cast<::PetscObject>( in ) );
    return 0;
}

PetscErrorCode _AMP_norm( Vec in, NormType type, PetscReal *ans )
{
    if ( type == NORM_1 )
        *ans = reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( in->data )->L1Norm();
    else if ( type == NORM_2 )
        *ans = reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( in->data )->L2Norm();
    else if ( type == NORM_INFINITY )
        *ans = reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( in->data )->maxNorm();
    else if ( type == NORM_1_AND_2 ) {
        *ans = reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( in->data )->L1Norm();
        *( ans + 1 ) =
            reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( in->data )->L2Norm();
    } else
        AMP_ERROR( "Unknown norm type" );
    if ( type != NORM_1_AND_2 ) {
        PetscObjectComposedDataSetReal(
            reinterpret_cast<::PetscObject>( in ), NormIds[type], ans[0] );
    }
    PetscObjectStateIncrease( reinterpret_cast<::PetscObject>( in ) );
    return 0;
}

bool _Verify_Memory( AMP::LinearAlgebra::Vector *p1, AMP::LinearAlgebra::Vector *p2 )
{
    for ( size_t i = 0; i != p1->numberOfDataBlocks(); i++ ) {
        if ( p1->getRawDataBlock<double>() == p2->getRawDataBlock<double>() )
            return false;
    }
    return true;
}

PetscErrorCode _AMP_duplicate( Vec in, Vec *out )
{

    AMP::LinearAlgebra::ManagedPetscVector *p =
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( in->data );
    AMP::LinearAlgebra::ManagedPetscVector *dup = p->petscDuplicate();
    AMP_ASSERT( _Verify_Memory( p, dup ) );
    *out = dup->getVec();
    return 0;
}

PetscErrorCode _AMP_duplicatevecs( Vec v, PetscInt num, Vec **vecArray )
{
    auto tvecArray = new Vec[num];
    for ( PetscInt i = 0; i != num; i++ )
        VecDuplicate( v, tvecArray + i );
    *vecArray = tvecArray;
    return 0;
}


PetscErrorCode _AMP_restorearray( Vec y, PetscScalar ** )
{
    PetscObjectStateIncrease( reinterpret_cast<::PetscObject>( y ) );
    return 0;
}

PetscErrorCode _AMP_getarray( Vec in, PetscScalar **out )
{
    AMP::LinearAlgebra::ManagedPetscVector *p =
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( in->data );
    *out = p->getRawDataBlock<PetscScalar>();
    return 0;
}

PetscErrorCode _AMP_getlocalsize( Vec in, PetscInt *out )
{
    AMP::LinearAlgebra::ManagedPetscVector *p =
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( in->data );
    *out = (PetscInt) p->getLocalSize();
    return 0;
}

PetscErrorCode _AMP_setfromoptions( Vec )
{
    AMP_ERROR( "This should not be thrown" );
    return 0;
}

PetscErrorCode _AMP_reciprocal( Vec v )
{
    AMP::LinearAlgebra::ManagedPetscVector *p =
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( v->data );
    p->reciprocal( *p );
    return 0;
}

PetscErrorCode _AMP_abs( Vec v )
{
    AMP::LinearAlgebra::ManagedPetscVector *p =
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( v->data );
    p->abs( *p );
    PetscObjectStateIncrease( reinterpret_cast<::PetscObject>( v ) );
    return 0;
}

PetscErrorCode _AMP_resetarray( Vec ) { return 0; }
PetscErrorCode _AMP_destroy( Vec v )
{
    AMP::LinearAlgebra::ManagedPetscVector *p =
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( v->data );
    if ( p->constructedWithPetscDuplicate() ) {
        delete p;
    }
    return 0;
}

PetscErrorCode _AMP_create( Vec ) { return 0; }


PetscErrorCode _AMP_set( Vec x, PetscScalar alpha )
{
    AMP::LinearAlgebra::ManagedPetscVector *p =
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( x->data );
    p->setToScalar( alpha );
    // petsc calls object state increase for this function
    return 0;
}

void reset_vec_ops( Vec t )
{

    // Then, replace the functions
    t->ops->duplicate       = _AMP_duplicate;
    t->ops->duplicatevecs   = _AMP_duplicatevecs;
    t->ops->destroyvecs     = _AMP_destroyvecs;
    t->ops->dot             = _AMP_dot;
    t->ops->mdot            = _AMP_mdot;
    t->ops->norm            = _AMP_norm;
    t->ops->tdot            = _AMP_tdot;
    t->ops->mtdot           = _AMP_mtdot;
    t->ops->scale           = _AMP_scale;
    t->ops->copy            = _AMP_copy;
    t->ops->set             = _AMP_set;
    t->ops->swap            = _AMP_swap;
    t->ops->axpy            = _AMP_axpy;
    t->ops->axpby           = _AMP_axpby;
    t->ops->maxpy           = _AMP_maxpy;
    t->ops->aypx            = _AMP_aypx;
    t->ops->waxpy           = _AMP_waxpy;
    t->ops->axpbypcz        = _AMP_axpbypcz;
    t->ops->pointwisemult   = _AMP_pointwisemult;
    t->ops->pointwisedivide = _AMP_pointwisedivide;
    t->ops->setvalues       = _AMP_setvalues;
    t->ops->assemblybegin   = _AMP_assemblybegin;
    t->ops->assemblyend     = _AMP_assemblyend;
    t->ops->getarray        = _AMP_getarray;
    t->ops->getsize         = _AMP_getsize;
    t->ops->getlocalsize    = _AMP_getlocalsize;
    t->ops->restorearray    = _AMP_restorearray;
    t->ops->max             = _AMP_max;
    t->ops->min             = _AMP_min;
    t->ops->setrandom       = _AMP_setrandom;


    t->ops->setoption = _AMP_setoption;


    t->ops->setvaluesblocked        = _AMP_setvaluesblocked;
    t->ops->destroy                 = _AMP_destroy;
    t->ops->view                    = _AMP_view;
    t->ops->placearray              = _AMP_placearray;
    t->ops->replacearray            = _AMP_replacearray;
    t->ops->dot_local               = _AMP_dot_local;
    t->ops->tdot_local              = _AMP_tdot_local;
    t->ops->norm_local              = _AMP_norm_local;
    t->ops->mdot_local              = _AMP_mdot_local;
    t->ops->mtdot_local             = _AMP_mtdot_local;
    t->ops->reciprocal              = _AMP_reciprocal;
    t->ops->conjugate               = _AMP_conjugate;
    t->ops->setlocaltoglobalmapping = _AMP_setlocaltoglobalmapping;
    t->ops->setvalueslocal          = _AMP_setvalueslocal;
    t->ops->resetarray              = _AMP_resetarray;
    t->ops->maxpointwisedivide      = _AMP_maxpointwisedivide;
    t->ops->load                    = _AMP_load;
    t->ops->pointwisemax            = _AMP_pointwisemax;
    t->ops->pointwisemaxabs         = _AMP_pointwisemaxabs;
    t->ops->pointwisemin            = _AMP_pointwisemin;
    t->ops->getvalues               = _AMP_getvalues;
    t->ops->sqrt                    = _AMP_sqrt;
    t->ops->abs                     = _AMP_abs;
    t->ops->exp                     = _AMP_exp;
    t->ops->log                     = _AMP_log;
    t->ops->shift                   = _AMP_shift;
    t->ops->create                  = _AMP_create;
#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
    t->ops->loadintovector       = _AMP_loadintovector;
    t->ops->loadintovectornative = _AMP_loadintovectornative;
    t->ops->viewnative           = _AMP_viewnative;
#endif

    /*** The following functions do not need to be overridden
      t->ops->setfromoptions = _AMP_setfromoptions;
     ***/
}


namespace AMP {
namespace LinearAlgebra {


void ManagedPetscVector::initPetsc()
{
    AMP_MPI comm =
        AMP::dynamic_pointer_cast<ManagedVectorParameters>( getParameters() )->d_Engine->getComm();
    VecCreate( comm.getCommunicator(), &d_petscVec );

    d_petscVec->data        = this;
    d_petscVec->petscnative = PETSC_FALSE;

    reset_vec_ops( d_petscVec );

#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
    PetscMapInitialize( comm.getCommunicator(), d_petscVec->map );
    PetscMapSetBlockSize( d_petscVec->map, 1 );
    PetscMapSetSize( d_petscVec->map, this->getGlobalSize() );
    PetscMapSetLocalSize( d_petscVec->map, this->getLocalSize() );
    d_petscVec->map->rstart = static_cast<PetscInt>( this->getDOFManager()->beginDOF() );
    d_petscVec->map->rend   = static_cast<PetscInt>( this->getDOFManager()->endDOF() );
#elif PETSC_VERSION_GE(3,2,0)
    PetscLayoutSetBlockSize( d_petscVec->map, 1 );
    PetscLayoutSetSize( d_petscVec->map, this->getGlobalSize() );
    PetscLayoutSetLocalSize( d_petscVec->map, this->getLocalSize() );
    PetscLayoutSetUp( d_petscVec->map );
#else
#error Not programmed for this version yet
#endif

    d_bMadeWithPetscDuplicate = false;

    const std::string my_name = "AMPManagedPetscVectorReal";
    int ierr                  = 0;

    if ( ( (PetscObject) d_petscVec )->type_name ) {
        ierr = PetscFree( ( (PetscObject) d_petscVec )->type_name );
    }

    ierr =
        PetscObjectChangeTypeName( reinterpret_cast<PetscObject>( d_petscVec ), my_name.c_str() );
    AMP_INSIST( ierr == 0, "PetscObjectChangeTypeName returned non-zero error code" );

    VecSetFromOptions( d_petscVec );
}


ManagedPetscVector::ManagedPetscVector( VectorParameters::shared_ptr params )
    : ManagedVector( params ), PetscVector()
{
    initPetsc();
    registerListener( this );
}


ManagedPetscVector::ManagedPetscVector( Vector::shared_ptr alias )
    : ManagedVector( alias ), PetscVector()
{
    initPetsc();
    auto firer = dynamic_pointer_cast<DataChangeFirer>( alias );
    if ( firer )
        firer->registerListener( this );
}


ManagedPetscVector::~ManagedPetscVector()
{
    int refct = ( ( (PetscObject) d_petscVec )->refct );
    if ( !d_bMadeWithPetscDuplicate ) {
        if ( refct > 1 )
            AMP_ERROR( "Deleting a vector still held by PETSc" );
        PETSC::vecDestroy( &d_petscVec );
    }
}


bool ManagedPetscVector::petscHoldsView() const
{
    int refct = ( ( (PetscObject) d_petscVec )->refct );
    if ( !d_bMadeWithPetscDuplicate && refct > 1 )
        return true;
    return false;
}


ManagedPetscVector *ManagedPetscVector::petscDuplicate()
{
    ManagedPetscVector *pAns = rawClone();
    pAns->setVariable( getVariable() );
    pAns->d_bMadeWithPetscDuplicate = true;
    return pAns;
}


void ManagedPetscVector::copyFromPetscVec( Vector &dest, Vec source )
{
    auto params = AMP::dynamic_pointer_cast<ManagedVectorParameters>(
        dynamic_cast<ManagedVector*>(&dest)->getParameters() );
    if ( !params )
        throw( "Incompatible vector types" );

    if ( sizeof( PetscInt ) < 8 )
        AMP_INSIST( dest.getGlobalSize() < 0x80000000,
                    "PETsc is compiled with 32-bit integers and "
                    "we are trying to use a vector with more "
                    "than 2^31 elements" );

    const VectorEngineParameters &eparams = *( params->d_Engine->getEngineParameters() );
    auto ids                              = new PetscInt[eparams.getLocalSize()];
    PetscInt begin                        = (PetscInt) eparams.beginDOF();
    PetscInt end                          = (PetscInt) eparams.endDOF();
    for ( PetscInt i   = begin; i < end; i++ )
        ids[i - begin] = i;
    VecGetValues( source, dest.getLocalSize(), ids, dest.getRawDataBlock<double>() );
    delete[] ids;
}


AMP::shared_ptr<AMP::LinearAlgebra::Vector> ManagedPetscVector::createFromPetscVec( Vec source,
                                                                                    AMP_MPI &comm )
{
#ifdef USE_EXT_TRILINOS
    PetscInt local_size, global_size, local_start, local_end;
    VecGetLocalSize( source, &local_size );
    VecGetSize( source, &global_size );
    VecGetOwnershipRange( source, &local_start, &local_end );
    AMP::shared_ptr<ManagedVectorParameters> t( new ManagedPetscVectorParameters() );
    VectorEngineParameters::shared_ptr ve_params(
        new EpetraVectorEngineParameters( local_size, global_size, comm ) );
    t->d_Engine = VectorEngine::shared_ptr( new EpetraVectorEngine(
        ve_params, VectorEngine::BufferPtr( new VectorEngine::Buffer( local_size ) ) ) );
    ManagedPetscVector *pRetVal_t =
        new ManagedPetscVector( AMP::dynamic_pointer_cast<VectorParameters>( t ) );
    Vector::shared_ptr pRetVal( pRetVal_t );
    return pRetVal;
#else
    AMP_ERROR( "General case not programmed yet" );
    NULL_USE( source );
    NULL_USE( comm );
    return AMP::shared_ptr<AMP::LinearAlgebra::Vector>();
#endif
}


void ManagedPetscVector::swapVectors( Vector &other )
{
    auto tmp = dynamic_cast<ManagedPetscVector*>( &other );
    AMP_ASSERT( tmp != nullptr );
    ParentVector::swapVectors( *tmp );
}

}
}
