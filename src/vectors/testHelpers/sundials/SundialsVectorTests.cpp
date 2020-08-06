#ifdef USE_EXT_SUNDIALS

#include "AMP/vectors/testHelpers/sundials/SundialsVectorTests.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/sundials/ManagedSundialsVector.h"
#include "AMP/vectors/sundials/SundialsVector.h"


namespace AMP {
namespace LinearAlgebra {


static inline N_Vector getVec( AMP::LinearAlgebra::Vector::shared_ptr vector )
{
    auto sundials = std::dynamic_pointer_cast<AMP::LinearAlgebra::SundialsVector>( vector );
    AMP_ASSERT( sundials != nullptr );
    return sundials->getNVector();
}


void SundialsVectorTests::testSundialsVector( AMP::UnitTest *ut )
{
    CloneSundialsVector( ut );
    LinearSumSundialsVector( ut );
    ConstSundialsVector( ut );
    ProdSundialsVector( ut );
    DivSundialsVector( ut );
    ScaleSundialsVector( ut );
    AbsSundialsVector( ut );
    InvSundialsVector( ut );
    AddConstSundialsVector( ut );
    DotProdSundialsVector( ut );
    MaxNormSundialsVector( ut );
    WRMSNormSundialsVector( ut );
    L1NormSundialsVector( ut );
}


void SundialsVectorTests::CloneSundialsVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getVector() );
    N_Vector vec_a                                     = getVec( vectora );
    N_Vector vec_b                                     = N_VClone( vec_a );
    AMP::LinearAlgebra::ManagedSundialsVector *vectorb = getVector( vec_b );
    bool pass                                          = true;
    for ( size_t i = 0; i != vectorb->numberOfDataBlocks(); i++ ) {
        if ( vectorb->getRawDataBlock<double>( i ) == vectora->getRawDataBlock<double>( i ) )
            pass = false;
    }
    if ( pass )
        utils->passes( "Clone created" );
    else
        utils->failure( "Failed to create clone" );
    N_VDestroy( vec_b );
    utils->passes( "N_VDestroy returned" );
}


void SundialsVectorTests::LinearSumSundialsVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorc( d_factory->getVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectord( d_factory->getVector() );

    N_Vector vec_a = getVec( vectora );
    N_Vector vec_b = getVec( vectorb );
    N_Vector vec_c = getVec( vectorc );

    vectora->setRandomValues();
    vectorb->setRandomValues();
    N_VLinearSum( .2, vec_a, .5, vec_b, vec_c );
    vectord->linearSum( .2, vectora, .5, vectorb, vectord );
    vectord->subtract( vectord, vectorc, vectord );
    if ( vectord->maxNorm() < 0.000001 )
        utils->passes( "random linear sum" );
    else
        utils->failure( "random linear sum" );
}


void SundialsVectorTests::ConstSundialsVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getVector() );
    N_Vector vec_a = getVec( vectora );
    N_VConst( 0., vec_a );
    double maxNorm = vectora->maxNorm();
    if ( maxNorm > 0 )
        utils->failure( "Nonzero inf norm" );
    else
        utils->passes( "Set vector to 0" );

    N_VConst( 1., vec_a );
    maxNorm       = vectora->maxNorm();
    double L1Norm = vectora->L1Norm();
    if ( ( maxNorm == 1. ) && ( L1Norm == (double) vectora->getGlobalSize() ) )
        utils->passes( "Set vector to 1" );
    else
        utils->failure( "Failed to set to 1" );
}


void SundialsVectorTests::ProdSundialsVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorc( d_factory->getVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectord( d_factory->getVector() );

    N_Vector vec_a = getVec( vectora );
    N_Vector vec_b = getVec( vectorb );
    N_Vector vec_c = getVec( vectorc );

    vectora->setRandomValues();
    vectorb->setRandomValues();
    N_VProd( vec_a, vec_b, vec_c );
    vectord->multiply( vectora, vectorb, vectord );
    vectord->subtract( vectorc, vectord, vectord );
    double norm = vectord->maxNorm();
    if ( norm < 0.000001 )
        utils->passes( "Products match" );
    else
        utils->failure( "Products are mis-matched" );
}


void SundialsVectorTests::DivSundialsVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorc( d_factory->getVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectord( d_factory->getVector() );

    N_Vector vec_a = getVec( vectora );
    N_Vector vec_b = getVec( vectorb );
    N_Vector vec_c = getVec( vectorc );

    vectora->setRandomValues();
    vectorb->setRandomValues();
    N_VDiv( vec_a, vec_b, vec_c );
    vectord->divide( vectora, vectorb, vectord );
    vectord->subtract( vectorc, vectord, vectord );
    if ( vectord->maxNorm() < 0.000001 )
        utils->passes( "Quotients match" );
    else
        utils->failure( "Quotients are mis-matched" );
}


void SundialsVectorTests::ScaleSundialsVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorc( d_factory->getVector() );

    N_Vector vec_a = getVec( vectora );
    N_Vector vec_b = getVec( vectorb );

    vectora->setRandomValues();
    N_VScale( 2.0, vec_a, vec_b );
    vectorc->scale( 2.0, vectora, vectorc );
    vectorc->subtract( vectorc, vectorb, vectorc );
    if ( vectorc->maxNorm() < 0.000001 )
        utils->passes( "Scalings match" );
    else
        utils->failure( "Scalings are mis-matched" );
}


void SundialsVectorTests::AbsSundialsVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorc( d_factory->getVector() );

    N_Vector vec_a = getVec( vectora );
    N_Vector vec_b = getVec( vectorb );

    vectora->setRandomValues();
    vectorc->setToScalar( 0.5 );
    vectora->subtract( vectora, vectorc, vectora );
    N_VAbs( vec_a, vec_b );
    vectorc->abs( vectora, vectorc );
    vectorc->subtract( vectorc, vectorb, vectorc );
    if ( vectorc->maxNorm() < 0.000001 )
        utils->passes( "Values match" );
    else
        utils->failure( "Values are mis-matched" );
}


void SundialsVectorTests::InvSundialsVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorc( d_factory->getVector() );

    N_Vector vec_a = getVec( vectora );
    N_Vector vec_b = getVec( vectorb );

    vectora->setRandomValues();
    N_VInv( vec_a, vec_b );
    vectorc->reciprocal( vectora, vectorc );
    vectorc->subtract( vectorc, vectorb, vectorc );
    if ( vectorc->maxNorm() < 0.000001 )
        utils->passes( "Scalings match" );
    else
        utils->failure( "Scalings are mis-matched" );
}


void SundialsVectorTests::AddConstSundialsVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorc( d_factory->getVector() );

    N_Vector vec_a = getVec( vectora );
    N_Vector vec_b = getVec( vectorb );

    vectora->setRandomValues();
    N_VAddConst( vec_a, .3, vec_b );
    vectorc->addScalar( vectora, .3, vectorc );
    vectorc->subtract( vectorb, vectorc, vectorc );
    double norm = vectorc->maxNorm();
    if ( norm < 0.00000001 )
        utils->passes( "N_VAddConst" );
    else
        utils->failure( "N_VAddConst" );
}


void SundialsVectorTests::DotProdSundialsVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getVector() );

    N_Vector vec_a = getVec( vectora );
    N_Vector vec_b = getVec( vectorb );

    vectora->setRandomValues();
    vectorb->setRandomValues();
    double d1 = N_VDotProd( vec_a, vec_b );
    double d2 = vectora->dot( vectorb, vectora );
    if ( fabs( d1 - d2 ) < 0.00000001 )
        utils->passes( "N_VDotProd" );
    else
        utils->failure( "N_VDotProd" );
}


void SundialsVectorTests::MaxNormSundialsVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getVector() );

    N_Vector vec_a = getVec( vectora );

    vectora->setRandomValues();

    double d1 = N_VMaxNorm( vec_a );
    double d2 = vectora->maxNorm();
    if ( fabs( d1 - d2 ) < 0.00000001 )
        utils->passes( "N_VMaxNorm" );
    else
        utils->failure( "N_VMaxNorm" );
}


void SundialsVectorTests::WRMSNormSundialsVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorc( d_factory->getVector() );
    if ( !vectorc )
        utils->failure( "N_VWrmsNorm" );

    N_Vector vec_a = getVec( vectora );
    N_Vector vec_b = getVec( vectorb );

    vectora->setRandomValues();
    vectorb->setRandomValues();

    double d1 = N_VWrmsNorm( vec_a, vec_b );
    double d2 = AMP::LinearAlgebra::Vector::wrmsNorm( vectora, vectorb );
    if ( fabs( d1 - d2 ) < 0.00000001 )
        utils->passes( "N_VWrmsNorm" );
    else
        utils->failure( "N_VWrmsNorm" );
}


void SundialsVectorTests::MinSundialsVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getVector() );

    N_Vector vec_a = getVec( vectora );

    vectora->setRandomValues();

    double d1 = N_VMin( vec_a );
    double d2 = vectora->min();
    if ( fabs( d1 - d2 ) < 0.00000001 )
        utils->passes( "N_VMin" );
    else
        utils->failure( "N_VMin" );
}


void SundialsVectorTests::L1NormSundialsVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getVector() );

    N_Vector vec_a = getVec( vectora );

    vectora->setRandomValues();

    double d1 = N_VL1Norm( vec_a );
    double d2 = vectora->L1Norm();
    if ( fabs( d1 - d2 ) < 0.00000001 )
        utils->passes( "N_VL1Norm" );
    else
        utils->failure( "N_VL1Norm" );
}


void SundialsVectorTests::MinQuotientSundialsVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getVector() );

    N_Vector vec_a = getVec( vectora );
    N_Vector vec_b = getVec( vectorb );

    vectora->setRandomValues();
    vectorb->setRandomValues();

    double d1 = N_VMinQuotient( vec_a, vec_b );
    double d2 = AMP::LinearAlgebra::Vector::minQuotient( vectora, vectorb );
    if ( fabs( d1 - d2 ) < 0.00000001 )
        utils->passes( "N_VMinQuotient" );
    else
        utils->failure( "N_VMinQuotient" );
}
} // namespace LinearAlgebra
} // namespace AMP

/// \endcond
#endif
