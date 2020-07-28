#include "AMP/vectors/sundials/ManagedSundialsVector.h"
#include "AMP/vectors/data/VectorDataCPU.h"

#include "AMP/utils/UtilityMacros.h"


namespace AMP {
namespace LinearAlgebra {

/**
 * This macro extracts the content field of an N_Vector and casts it as a pointer to a
 * SundialsVector.
 */

//#define AMPVEC_CAST(v) (static_cast<ManagedSundialsVector*>(v->content))


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
ManagedSundialsVector::ManagedSundialsVector( VectorParameters::shared_ptr params )
    : ManagedVector( params ), SundialsVector()
{
    // Create N_Vector
    d_n_vector = (N_Vector) malloc( sizeof *d_n_vector );
    AMP_ASSERT( d_n_vector != nullptr );
    // Attach the content and the ops fields
    d_n_vector->content = this;
    d_n_vector->ops     = ManagedSundialsVector::createNVectorOps();
}
ManagedSundialsVector::ManagedSundialsVector( shared_ptr alias )
    : ManagedVector( alias ), SundialsVector()
{
    // Create N_Vector
    d_n_vector = (N_Vector) malloc( sizeof *d_n_vector );
    AMP_ASSERT( d_n_vector != nullptr );
    // Attach the content and the ops fields
    d_n_vector->content = this;
    d_n_vector->ops     = ManagedSundialsVector::createNVectorOps();
}


/************************************************************************
 * Destructor for SundialsVector                                         *
 * Frees the memory allocated for the member N_Vector and its ops field  *
 ************************************************************************/
ManagedSundialsVector::~ManagedSundialsVector()
{
    if ( d_n_vector ) {
        if ( d_n_vector->ops ) {
            free( d_n_vector->ops );
            d_n_vector->ops = nullptr;
        }
        free( d_n_vector );
        d_n_vector = nullptr;
    }
}


/************************************************************************
 * Clone the vector                                                      *
 ************************************************************************/
Vector::shared_ptr ManagedSundialsVector::cloneVector( const Variable::shared_ptr var ) const
{
    Vector::shared_ptr retVal( rawClone() );
    retVal->setVariable( var );
    return retVal;
}
ManagedSundialsVector *ManagedSundialsVector::rawClone() const
{
    auto p      = std::make_shared<ManagedSundialsVectorParameters>();
    p->d_Buffer = d_vBuffer->cloneData();
    if ( !p->d_Buffer ) {
        auto vec    = std::dynamic_pointer_cast<Vector>( d_Engine );
        auto vec2   = vec->cloneVector();
        p->d_Buffer = std::dynamic_pointer_cast<VectorData>( vec2 );
        p->d_Engine = std::dynamic_pointer_cast<VectorOperations>( vec2 );
    } else {
        p->d_Engine = cloneVectorEngine( p->d_Buffer );
    }
    p->d_CommList   = getCommunicationList();
    p->d_DOFManager = getDOFManager();
    auto retVal     = new ManagedSundialsVector( p );
    return retVal;
}


/**
 * Creates ops, the structure of vector operations which gets attached to the member N_Vector
 * Functions with no_impl in their names are not implemented at the moment
 */
N_Vector_Ops ManagedSundialsVector::createNVectorOps()
{
    N_Vector_Ops ops;
    ops               = (N_Vector_Ops) malloc( sizeof( struct _generic_N_Vector_Ops ) );
    ops->nvclone      = cloneVector_AMP;
    ops->nvcloneempty = cloneempty_no_impl;
    ops->nvdestroy    = freeVectorComponents_AMP;
    // ops->nvspace         = space_no_impl;
    ops->nvspace           = nullptr;
    ops->nvgetarraypointer = getarraypointer_no_impl;
    ops->nvsetarraypointer = setarraypointer_no_impl;
    ops->nvlinearsum       = linearSum_AMP;
    ops->nvconst           = setToScalar_AMP;
    ops->nvprod            = multiply_AMP;
    ops->nvdiv             = divide_AMP;
    ops->nvscale           = scale_AMP;
    ops->nvabs             = abs_AMP;
    ops->nvinv             = reciprocal_AMP;
    ops->nvaddconst        = addScalar_AMP;
    ops->nvdotprod         = dot_AMP;
    ops->nvmaxnorm         = maxNorm_AMP;
    ops->nvwrmsnorm        = WRMSNorm_AMP;
    ops->nvwrmsnormmask    = WRMSNormMask_AMP;
    ops->nvmin             = min_AMP;
    ops->nvwl2norm         = wl2norm_no_impl;
    ops->nvl1norm          = L1Norm_AMP;
    ops->nvcompare         = compare_no_impl;
    ops->nvinvtest         = invtest_no_impl;
    ops->nvconstrmask      = constrmask_no_impl;
    ops->nvminquotient     = minquotient_AMP;
    return ops;
}


N_Vector ManagedSundialsVector::cloneVector_AMP( N_Vector n_vector )
{
    /**
     * Extracts the content filed of n_vector
     */
    auto *srcAMPVector = static_cast<ManagedSundialsVector *>( n_vector->content );
    ManagedSundialsVector *newSundialsVector = srcAMPVector->rawClone();

    newSundialsVector->setVariable( srcAMPVector->getVariable() );
    return ( newSundialsVector->getNVector() );
}

N_Vector ManagedSundialsVector::cloneempty_no_impl( N_Vector )
{
    AMP_ERROR( "nvcloneempty not implemented" );
    return N_Vector();
}

void ManagedSundialsVector::freeVectorComponents_AMP( N_Vector v )
{

    auto *ptr = static_cast<ManagedSundialsVector *>( v->content );
    delete ptr;
}

realtype *ManagedSundialsVector::getarraypointer_no_impl( N_Vector )
{
    AMP_ERROR( "nvgetarraypointer not implemented" );
    return nullptr;
}

void ManagedSundialsVector::setarraypointer_no_impl( realtype *, N_Vector )
{
    AMP_ERROR( "nvsetarraypointer not implemented" );
}

/**
 * Set z = a*x + b*y
 */
void ManagedSundialsVector::linearSum_AMP(
    realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z )
{
    Vector *px = static_cast<ManagedSundialsVector *>( x->content );
    Vector *py = static_cast<ManagedSundialsVector *>( y->content );

    ( static_cast<ManagedSundialsVector *>( z->content ) )->linearSum( a, *px, b, *py );
}

/**
 * Set each entry of z to c
 */

void ManagedSundialsVector::setToScalar_AMP( realtype c, N_Vector z )
{
    ( static_cast<ManagedSundialsVector *>( z->content ) )->setToScalar( c );
}

/**
 * Set z_i = x_i * y_i
 */

void ManagedSundialsVector::multiply_AMP( N_Vector x, N_Vector y, N_Vector z )
{
    Vector *px = static_cast<ManagedSundialsVector *>( x->content );
    Vector *py = static_cast<ManagedSundialsVector *>( y->content );
    Vector *pz = static_cast<ManagedSundialsVector *>( z->content );

    pz->multiply( *px, *py );
}

/**
 * Set z_i = x_i / y_i
 */

void ManagedSundialsVector::divide_AMP( N_Vector x, N_Vector y, N_Vector z )
{
    Vector *px = static_cast<ManagedSundialsVector *>( x->content );
    Vector *py = static_cast<ManagedSundialsVector *>( y->content );
    Vector *pz = static_cast<ManagedSundialsVector *>( z->content );

    pz->divide( *px, *py );
}

/**
 * Set z_i = c * x_i
 */
void ManagedSundialsVector::scale_AMP( realtype c, N_Vector x, N_Vector z )
{
    Vector *px = static_cast<ManagedSundialsVector *>( x->content );

    ( static_cast<ManagedSundialsVector *>( z->content ) )->scale( c, *px );
}

/**
 * Set z_i = | x_i |
 */

void ManagedSundialsVector::abs_AMP( N_Vector x, N_Vector z )
{
    Vector *px = static_cast<ManagedSundialsVector *>( x->content );

    ( static_cast<ManagedSundialsVector *>( z->content ) )->abs( *px );
}

/**
 * Set z_i = 1.0 / x_i
 */
void ManagedSundialsVector::reciprocal_AMP( N_Vector x, N_Vector z )
{
    Vector *px = static_cast<ManagedSundialsVector *>( x->content );

    ( static_cast<ManagedSundialsVector *>( z->content ) )->reciprocal( *px );
}

/**
 * Set z_i = x_i + b
 */

void ManagedSundialsVector::addScalar_AMP( N_Vector x, realtype b, N_Vector z )
{
    Vector *px           = static_cast<ManagedSundialsVector *>( x->content );
    Vector::shared_ptr t = px->cloneVector();
    t->setToScalar( b );

    ( static_cast<ManagedSundialsVector *>( z->content ) )->add( *px, *t );
}

/**
 * Returns the dot product of x and y
 */
realtype ManagedSundialsVector::dot_AMP( N_Vector x, N_Vector y )
{
    double dotprod;

    Vector *px = static_cast<ManagedSundialsVector *>( x->content );

    dotprod = ( static_cast<ManagedSundialsVector *>( y->content ) )->dot( *px );
    return ( dotprod );
}

/**
 * Returns the maximum norm of the input vector
 */
realtype ManagedSundialsVector::maxNorm_AMP( N_Vector x )
{
    double maxnorm;

    maxnorm = ( static_cast<ManagedSundialsVector *>( x->content ) )->maxNorm();

    return ( maxnorm );
}

/**
 * Returns the minimum entry of the input vector
 */
realtype ManagedSundialsVector::min_AMP( N_Vector x )
{
    double min;

    min = ( static_cast<ManagedSundialsVector *>( x->content ) )->min();

    return ( min );
}

/**
 * Returns the L1 norm of the input vector
 */
realtype ManagedSundialsVector::L1Norm_AMP( N_Vector x )
{
    double l1norm;
    l1norm = ( static_cast<ManagedSundialsVector *>( x->content ) )->L1Norm();
    return ( l1norm );
}

/**
 * This function is supposed to return the weighted RMS-norm of the input vector x, with w a
 * weighting vector.
 * Proper implementation needs to be provided.
 */

realtype ManagedSundialsVector::WRMSNorm_AMP( N_Vector x, N_Vector w )
{
    Vector *px = static_cast<ManagedSundialsVector *>( x->content );
    Vector *pw = static_cast<ManagedSundialsVector *>( w->content );
    // is this OK?
    return Vector::wrmsNorm( *px, *pw );
}


realtype ManagedSundialsVector::WRMSNormMask_AMP( N_Vector x, N_Vector w, N_Vector id )
{
    Vector *px = static_cast<ManagedSundialsVector *>( x->content );
    Vector *pw = static_cast<ManagedSundialsVector *>( w->content );
    Vector *pm = static_cast<ManagedSundialsVector *>( id->content );

    // is this OK?
    return Vector::wrmsNormMask( *px, *pw, *pm );
}


realtype ManagedSundialsVector::wl2norm_no_impl( N_Vector, N_Vector )
{
    AMP_ERROR( "nvwl2norm not implemented" );
    return 0;
}

void ManagedSundialsVector::compare_no_impl( realtype, N_Vector, N_Vector )
{
    AMP_ERROR( "nvcompare not implemented" );
}

booleantype ManagedSundialsVector::invtest_no_impl( N_Vector, N_Vector )
{
    AMP_ERROR( "nvinvtest not implemented" );
    return false;
}

booleantype ManagedSundialsVector::constrmask_no_impl( N_Vector, N_Vector, N_Vector )
{
    AMP_ERROR( "nvconstrmask not implemented" );
    return false;
}

realtype ManagedSundialsVector::minquotient_AMP( N_Vector x, N_Vector w )
{
    Vector *px = static_cast<ManagedSundialsVector *>( x->content );
    Vector *pw = static_cast<ManagedSundialsVector *>( w->content );

    // is this OK?
    return ( px->minQuotient( *px, *pw ) );
}

ManagedVector *ManagedSundialsVector::getNewRawPtr() const
{
    return new ManagedSundialsVector(
        std::dynamic_pointer_cast<VectorParameters>( d_pParameters ) );
}

std::string ManagedSundialsVector::type() const
{
    std::string retVal = "Managed SUNDIALS Vector";
    retVal += ManagedVector::type();
    return retVal;
}

void ManagedSundialsVector::assemble() {}

} // namespace LinearAlgebra
} // namespace AMP
