#ifdef USE_EXT_PETSC

#include "vectors/MultiVector.h"
#include "vectors/petsc/PetscHelpers.h"
#include "vectors/testHelpers/PetscVectorTests.h"
#include "utils/UnitTest.h"

#include "string"
#include <algorithm>


namespace AMP {
namespace LinearAlgebra {


void checkPetscError( AMP::UnitTest *utils, PetscErrorCode i )
{
    if ( i ) {
        char *ans;
        PetscErrorMessage( i, PETSC_NULL, &ans );
        utils->failure( ans );
        delete[] ans;
    }
}


void PetscVectorTests::testPetscVector( AMP::UnitTest *ut )
{
    InstantiatePetscVectors( ut );
    DuplicatePetscVector( ut );
    StaticCopyPetscVector( ut );
    StaticDuplicatePetscVector( ut );
    CopyPetscVector( ut );
    VerifyDotPetscVector( ut );
    VerifyNormsPetscVector( ut );
    VerifyScalePetscVector( ut );
    VerifyAXPYPetscVector( ut );
    VerifyAbsPetscVector( ut );
    VerifyMaxPointwiseDividePetscVector( ut );
    VerifyGetSizePetscVector( ut );
    VerifySwapPetscVector( ut );
    VerifyAXPBYPetscVector( ut );
    VerifySetPetscVector( ut );
    VerifySetRandomPetscVector( ut );
    VerifySqrtPetscVector( ut );
    VerifyPointwiseMultPetscVector( ut );
    VerifyPointwiseDividePetscVector( ut );
    VerifyPointwiseMaxPetscVector( ut );
    VerifyPointwiseMinPetscVector( ut );
    VerifyPointwiseMaxAbsPetscVector( ut );
    VerifyLogPetscVector( ut );
    VerifyExpPetscVector( ut );
    VerifyAYPXPetscVector( ut );
    VerifyAXPBYPCZPetscVector( ut );
}

void PetscVectorTests::InstantiatePetscVectors( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr native_vector( d_factory->getNativeVector() );
    if ( native_vector )
        utils->passes( "native created" );
    else
        utils->failure( "native created" );
    AMP::LinearAlgebra::Vector::shared_ptr managed_vector = d_factory->getManagedVector();
    if ( managed_vector )
        utils->passes( "managed created" );
    else
        utils->failure( "managed created" );
    utils->passes( "native destroyed" );
}


void PetscVectorTests::Bug_612( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getManagedVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getManagedVector() );

    Vec veca = vectora->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    Vec vecb = vectorb->castTo<AMP::LinearAlgebra::PetscVector>().getVec();

    vectora->setToScalar( 5.0 );
    vectorb->setToScalar( 0.0 );

    double result;
    VecMaxPointwiseDivide( veca, vecb, &result );
    if ( result == 5.0 )
        utils->passes( "Correct computation 1" );
    else
        utils->failure( "Incorrect computation 1" );

    vectorb->setToScalar( 5.0 );
    vectorb->getRawDataBlock<double>( 0 )[0] = 0.0;

    VecMaxPointwiseDivide( veca, vecb, &result );
    if ( result == 5.0 )
        utils->passes( "Correct computation 2" );
    else
        utils->failure( "Incorrect computation 2" );

    vectorb->setToScalar( 0.5 );
    vectorb->getRawDataBlock<double>( 0 )[0] = 0.0;

    VecMaxPointwiseDivide( veca, vecb, &result );
    if ( result == 10.0 )
        utils->passes( "Correct computation 3" );
    else
        utils->failure( "Incorrect computation 3" );
}


void PetscVectorTests::DuplicatePetscVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getManagedVector() );
    if ( vectora->isA<AMP::LinearAlgebra::NativePetscVector>() )
        return;

    vectora->setVariable( AMP::LinearAlgebra::Variable::shared_ptr(
        new AMP::LinearAlgebra::Variable( "dummy_variable" ) ) );
    Vec petsc_vec = vectora->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    Vec another_vec;
    checkPetscError( utils, VecDuplicate( petsc_vec, &another_vec ) );
    AMP::LinearAlgebra::ManagedPetscVector *dup =
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( another_vec->data );
    utils->passes( "managed duplicated" );
    if ( ( dup->getGlobalSize() == vectora->getGlobalSize() ) &&
         ( dup->getLocalSize() == vectora->getLocalSize() ) )
        utils->passes( "Allocated sizes are the same" );
    else
        utils->failure( "Allocated sizes are different" );
    if ( *( vectora->getVariable() ) == *( dup->getVariable() ) )
        utils->passes( "Associated variables are the same" );
    else
        utils->passes( "Associated variables are different" );
    checkPetscError( utils, PETSC::vecDestroy( &another_vec ) );
    utils->passes( "managed duplicated destroyed" );

    if ( vectora->isA<AMP::LinearAlgebra::MultiVector>() ) {
        AMP::LinearAlgebra::Vector::shared_ptr b =
            AMP::LinearAlgebra::PetscVector::view( vectora );
        bool passed = true;
        for ( size_t i = 0; i != b->numberOfDataBlocks(); i++ ) {
            if ( b->getRawDataBlock<double>( i ) == vectora->getRawDataBlock<double>( i ) ) {
                passed = false;
            }
        }
        if ( passed )
            utils->passes( "VecDuplicate of a multivector allocates new space" );
        else
            utils->passes( "VecDuplicate of a multivector does not allocate new space" );
    }
}


void PetscVectorTests::StaticDuplicatePetscVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vector_dup = d_factory->getManagedVector();
    vector_dup->setToScalar( 1. );
    double t = vector_dup->L1Norm();
    if ( vector_dup->getGlobalSize() == vectora->getGlobalSize() )
        utils->passes( "global size equality" );
    else
        utils->failure( "global size equality" );

    if ( vector_dup->getLocalSize() == vectora->getLocalSize() )
        utils->passes( "local size equality" );
    else
        utils->failure( "local size equality" );

    if ( fabs( t - (double) vector_dup->getGlobalSize() ) < 0.0000001 )
        utils->passes( "trivial set data" );
    else
        utils->failure( "trivial set data" );
}


void PetscVectorTests::StaticCopyPetscVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getNativeVector() );
    vectora->setToScalar( 1. );
    AMP::LinearAlgebra::Vector::shared_ptr vector_dup = d_factory->getManagedVector();
    vector_dup->copyVector( vectora );
    double t = vector_dup->L1Norm();
    if ( vector_dup->getGlobalSize() == vectora->getGlobalSize() )
        utils->passes( "global size equality" );
    else
        utils->failure( "global size equality" );

    if ( vector_dup->getLocalSize() == vectora->getLocalSize() )
        utils->passes( "local size equality" );
    else
        utils->failure( "local size equality" );

    if ( fabs( t - (double) vector_dup->getGlobalSize() ) < 0.0000001 )
        utils->passes( "trivial copy data" );
    else
        utils->failure( "trivial copy data" );
}


void PetscVectorTests::CopyPetscVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorc( d_factory->getNativeVector() );
    vectora->setRandomValues();
    vectorb->copyVector( vectora );
    vectorc->subtract( vectora, vectorb );
    if ( vectorc->maxNorm() < 0.0000001 )
        utils->passes( "native petsc copy" );
    else
        utils->failure( "native petsc copy" );
}


void PetscVectorTests::VerifyPointwiseMaxAbsPetscVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorc( d_factory->getNativeVector() );

    vectora->setRandomValues();
    vectorb->setToScalar( .65 );

    AMP::LinearAlgebra::Vector::shared_ptr vectord = d_factory->getManagedVector();
    AMP::LinearAlgebra::Vector::shared_ptr vectore = d_factory->getManagedVector();
    AMP::LinearAlgebra::Vector::shared_ptr vectorf = d_factory->getManagedVector();

    vectord->copyVector( vectora );
    vectore->setToScalar( .65 );

    Vec veca, vecb, vecc, vecd, vece, vecf;

    veca = vectora->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecb = vectorb->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecc = vectorc->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecd = vectord->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vece = vectore->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecf = vectorf->castTo<AMP::LinearAlgebra::PetscVector>().getVec();

    checkPetscError( utils, VecPointwiseMaxAbs( vecc, veca, vecb ) );
    checkPetscError( utils, VecPointwiseMaxAbs( vecf, vecd, vece ) );

    if ( vectorc->equals( vectorf ) )
        utils->passes( "VecPointwiseMaxAbs test" );
    else
        utils->failure( "VecPointwiseMaxAbs test" );
}


void PetscVectorTests::VerifyPointwiseMaxPetscVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorc( d_factory->getNativeVector() );

    vectora->setRandomValues();
    vectorb->setToScalar( .35 );

    AMP::LinearAlgebra::Vector::shared_ptr vectord = d_factory->getManagedVector();
    AMP::LinearAlgebra::Vector::shared_ptr vectore = d_factory->getManagedVector();
    AMP::LinearAlgebra::Vector::shared_ptr vectorf = d_factory->getManagedVector();

    vectord->copyVector( vectora );
    vectore->setToScalar( .35 );

    Vec veca, vecb, vecc, vecd, vece, vecf;

    veca = vectora->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecb = vectorb->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecc = vectorc->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecd = vectord->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vece = vectore->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecf = vectorf->castTo<AMP::LinearAlgebra::PetscVector>().getVec();

    checkPetscError( utils, VecPointwiseMax( vecc, veca, vecb ) );
    checkPetscError( utils, VecPointwiseMax( vecf, vecd, vece ) );

    if ( vectorc->equals( vectorf ) )
        utils->passes( "VecPointwiseMax test" );
    else
        utils->failure( "VecPointwiseMax test" );
}


void PetscVectorTests::VerifyPointwiseMinPetscVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorc( d_factory->getNativeVector() );

    vectora->setRandomValues();
    vectorb->setToScalar( .35 );

    AMP::LinearAlgebra::Vector::shared_ptr vectord = d_factory->getManagedVector();
    AMP::LinearAlgebra::Vector::shared_ptr vectore = d_factory->getManagedVector();
    AMP::LinearAlgebra::Vector::shared_ptr vectorf = d_factory->getManagedVector();

    vectord->copyVector( vectora );
    vectore->setToScalar( .35 );

    Vec veca, vecb, vecc, vecd, vece, vecf;

    veca = vectora->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecb = vectorb->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecc = vectorc->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecd = vectord->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vece = vectore->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecf = vectorf->castTo<AMP::LinearAlgebra::PetscVector>().getVec();

    checkPetscError( utils, VecPointwiseMin( vecc, veca, vecb ) );
    checkPetscError( utils, VecPointwiseMin( vecf, vecd, vece ) );

    if ( vectorc->equals( vectorf ) )
        utils->passes( "VecPointwiseMin test" );
    else
        utils->failure( "VecPointwiseMin test" );
}


void PetscVectorTests::VerifyAXPBYPCZPetscVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorc( d_factory->getNativeVector() );

    vectora->setRandomValues();
    vectorb->setToScalar( -.5 );
    vectorc->setToScalar( 3.45678 );

    AMP::LinearAlgebra::Vector::shared_ptr vectord = d_factory->getManagedVector();
    AMP::LinearAlgebra::Vector::shared_ptr vectore = d_factory->getManagedVector();
    AMP::LinearAlgebra::Vector::shared_ptr vectorf = d_factory->getManagedVector();

    vectord->copyVector( vectora );
    vectore->setToScalar( -.5 );
    vectorf->setToScalar( 3.45678 );

    Vec veca, vecb, vecc, vecd, vece, vecf;

    veca = vectora->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecb = vectorb->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecc = vectorc->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecd = vectord->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vece = vectore->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecf = vectorf->castTo<AMP::LinearAlgebra::PetscVector>().getVec();

    checkPetscError( utils,
                                     VecAXPBYPCZ( veca, 3.14159, 1.414, 2.1727, vecb, vecc ) );
    checkPetscError( utils,
                                     VecAXPBYPCZ( vecd, 3.14159, 1.414, 2.1727, vece, vecf ) );

    if ( vectora->equals( vectord ) )
        utils->passes( "VecAXPBYPCZ test" );
    else
        utils->failure( "VecAXPBYPCZ test" );
}


void PetscVectorTests::VerifyAYPXPetscVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getNativeVector() );

    vectora->setRandomValues();
    vectorb->setToScalar( -.5 );

    AMP::LinearAlgebra::Vector::shared_ptr vectorc = d_factory->getManagedVector();
    AMP::LinearAlgebra::Vector::shared_ptr vectord = d_factory->getManagedVector();

    vectorc->copyVector( vectora );
    vectord->copyVector( vectorb );

    Vec veca, vecb, vecc, vecd;

    veca = vectora->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecb = vectorb->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecc = vectorc->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecd = vectord->castTo<AMP::LinearAlgebra::PetscVector>().getVec();

    checkPetscError( utils, VecAYPX( veca, 2, vecb ) );
    checkPetscError( utils, VecAYPX( vecc, 2, vecd ) );

    if ( vectora->equals( vectorc ) )
        utils->passes( "VecAYPX test" );
    else
        utils->failure( "VecAYPX test" );
}


void PetscVectorTests::VerifyExpPetscVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getNativeVector() );

    vectora->setRandomValues();
    vectora->abs( vectora );

    AMP::LinearAlgebra::Vector::shared_ptr vectorb = d_factory->getManagedVector();
    vectorb->copyVector( vectora );

    Vec veca, vecb;

    veca = vectora->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecb = vectorb->castTo<AMP::LinearAlgebra::PetscVector>().getVec();

    checkPetscError( utils, VecExp( veca ) );
    checkPetscError( utils, VecExp( vecb ) );

    if ( vectora->equals( vectorb ) )
        utils->passes( "VecExp test" );
    else
        utils->failure( "VecExp test" );
}


void PetscVectorTests::VerifyLogPetscVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getNativeVector() );

    vectora->setRandomValues();
    vectora->abs( vectora );

    AMP::LinearAlgebra::Vector::shared_ptr vectorb = d_factory->getManagedVector();
    vectorb->copyVector( vectora );

    Vec veca, vecb;

    veca = vectora->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecb = vectorb->castTo<AMP::LinearAlgebra::PetscVector>().getVec();

    checkPetscError( utils, VecLog( veca ) );
    checkPetscError( utils, VecLog( vecb ) );

    if ( vectora->equals( vectorb ) )
        utils->passes( "VecLog test" );
    else
        utils->failure( "VecLog test" );
}


void PetscVectorTests::VerifyNormsPetscVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getNativeVector() );
    Vec veca;
    veca = vectora->castTo<AMP::LinearAlgebra::PetscVector>().getVec();

    vectora->setRandomValues();
    double l1norm_a1, l1norm_a2;
    double l2norm_a1, l2norm_a2;
    double infnorm_a1, infnorm_a2;
    checkPetscError( utils, VecNorm( veca, NORM_1, &l1norm_a1 ) );
    checkPetscError( utils, VecNorm( veca, NORM_2, &l2norm_a1 ) );
    checkPetscError( utils, VecNorm( veca, NORM_INFINITY, &infnorm_a1 ) );
    l1norm_a2  = vectora->L1Norm();
    l2norm_a2  = vectora->L2Norm();
    infnorm_a2 = vectora->maxNorm();
    if ( l1norm_a1 == l1norm_a2 ) // These should be identical, since same method called
        utils->passes( "l1 norm: native norm equals interface norm for native vector" );
    else
        utils->failure(
            "l1 norm: native norm does not equal interface norm for native vector" );
    if ( l2norm_a1 == l2norm_a2 ) // These should be identical, since same method called
        utils->passes( "l2 norm: native norm equals interface norm for native vector" );
    else
        utils->failure(
            "l2 norm: native norm does not equal interface norm for native vector" );
    if ( infnorm_a1 == infnorm_a2 ) // These should be identical, since same method called
        utils->passes( "inf norm: native norm equals interface norm for native vector" );
    else
        utils->failure(
            "inf norm: native norm does not equal interface norm for native vector" );

    AMP::LinearAlgebra::Vector::shared_ptr vectorc = d_factory->getManagedVector();
    vectorc->copyVector( vectora );


    Vec vecc;
    vecc = vectorc->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    double l1norm_c1, l1norm_c2;
    double l2norm_c1, l2norm_c2;
    double infnorm_c1, infnorm_c2;
    checkPetscError( utils, VecNorm( vecc, NORM_1, &l1norm_c1 ) );
    checkPetscError( utils, VecNorm( vecc, NORM_2, &l2norm_c1 ) );
    checkPetscError( utils, VecNorm( vecc, NORM_INFINITY, &infnorm_c1 ) );
    l1norm_c2  = vectorc->L1Norm();
    l2norm_c2  = vectorc->L2Norm();
    infnorm_c2 = vectorc->maxNorm();
    if ( l1norm_c1 == l1norm_c2 ) // These should be identical, since same method called
        utils->passes( "l1 norm: native norm equals interface norm for managed vector" );
    else
        utils->failure(
            "l1 norm: native norm does not equal interface norm for managed vector" );
    if ( l2norm_c1 == l2norm_c2 ) // These should be identical, since same method called
        utils->passes( "l2 norm: native norm equals interface norm for managed vector" );
    else
        utils->failure(
            "l2 norm: native norm does not equal interface norm for managed vector" );
    if ( infnorm_c1 == infnorm_c2 ) // These should be identical, since same method called
        utils->passes( "inf norm: native norm equals interface norm for managed vector" );
    else
        utils->failure(
            "inf norm: native norm does not equal interface norm for managed vector" );
    if ( fabs( l1norm_a1 - l1norm_c1 ) < 0.0000001 )
        utils->passes( "l1 norms equal for managed and native petsc vectors" );
    else
        utils->passes( "l1 norms not equal for managed and native petsc vectors" );
    if ( fabs( l2norm_a1 - l2norm_c1 ) < 0.0000001 )
        utils->passes( "l2 norms equal for managed and native petsc vectors" );
    else
        utils->passes( "l2 norms not equal for managed and native petsc vectors" );
    if ( fabs( infnorm_a1 - infnorm_c1 ) < 0.0000001 )
        utils->passes( "max norms equal for managed and native petsc vectors" );
    else
        utils->passes( "max norms not equal for managed and native petsc vectors" );
    checkPetscError( utils, VecNormBegin( vecc, NORM_1, &l1norm_c1 ) );
    checkPetscError( utils, VecNormBegin( vecc, NORM_2, &l2norm_c1 ) );
    checkPetscError( utils, VecNormBegin( vecc, NORM_INFINITY, &infnorm_c1 ) );
    checkPetscError( utils, VecNormEnd( vecc, NORM_1, &l1norm_c1 ) );
    checkPetscError( utils, VecNormEnd( vecc, NORM_2, &l2norm_c1 ) );
    checkPetscError( utils, VecNormEnd( vecc, NORM_INFINITY, &infnorm_c1 ) );
    l1norm_c2  = vectorc->L1Norm();
    l2norm_c2  = vectorc->L2Norm();
    infnorm_c2 = vectorc->maxNorm();
    if ( fabs( l1norm_c1 - l1norm_c2 ) <
         0.00001 ) // These should be identical, since same method called
        utils->passes( "l1 norm: native norm equals interface norm for managed vector" );
    else
        utils->failure(
            "l1 norm: native norm does not equal interface norm for managed vector" );
    if ( fabs( l2norm_c1 - l2norm_c2 ) <
         0.00001 ) // These should be identical, since same method called
        utils->passes( "l2 norm: native norm equals interface norm for managed vector" );
    else
        utils->failure(
            "l2 norm: native norm does not equal interface norm for managed vector" );
    if ( fabs( infnorm_c1 - infnorm_c2 ) <
         0.00001 ) // These should be identical, since same method called
        utils->passes( "inf norm: native norm equals interface norm for managed vector" );
    else
        utils->failure(
            "inf norm: native norm does not equal interface norm for managed vector" );
    if ( fabs( l1norm_a1 - l1norm_c1 ) < 0.0000001 )
        utils->passes( "l1 norms equal for managed and native petsc vectors" );
    else
        utils->passes( "l1 norms not equal for managed and native petsc vectors" );
    if ( fabs( l2norm_a1 - l2norm_c1 ) < 0.0000001 )
        utils->passes( "l2 norms equal for managed and native petsc vectors" );
    else
        utils->passes( "l2 norms not equal for managed and native petsc vectors" );
    if ( fabs( infnorm_a1 - infnorm_c1 ) < 0.0000001 )
        utils->passes( "max norms equal for managed and native petsc vectors" );
    else
        utils->passes( "max norms not equal for managed and native petsc vectors" );
}


void PetscVectorTests::VerifyAXPBYPetscVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getNativeVector() );

    vectora->setRandomValues();
    vectorb->setToScalar( 2.0 );

    AMP::LinearAlgebra::Vector::shared_ptr vectorc = d_factory->getManagedVector();
    vectorc->copyVector( vectora );
    AMP::LinearAlgebra::Vector::shared_ptr vectord = d_factory->getManagedVector();
    vectord->copyVector( vectorb );


    Vec veca, vecb, vecc, vecd;
    veca = vectora->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecb = vectorb->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecc = vectorc->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecd = vectord->castTo<AMP::LinearAlgebra::PetscVector>().getVec();

    checkPetscError( utils, VecAXPBY( veca, 1.234, 2.345, vecb ) );
    checkPetscError( utils, VecAXPBY( vecc, 1.234, 2.345, vecd ) );

    if ( vectora->L2Norm() != 0 )
        utils->passes( "Non-trivial axpby computed" );
    else
        utils->failure( "Trivial axpby computed" );
    if ( vectora->equals( vectorc ) )
        utils->passes( "Native axpby matches managed axpby" );
    else
        utils->failure( "Native axpby does not match managed axpby" );

    vectora->axpby( 1.234, 2.345, vectorb );
    vectorc->axpby( 1.234, 2.345, vectord );

    if ( vectora->L2Norm() != 0 )
        utils->passes( "Non-trivial axpby computed" );
    else
        utils->failure( "Trivial axpby computed" );
    if ( vectora->equals( vectorc ) )
        utils->passes( "Native axpby matches managed axpby" );
    else
        utils->failure( "Native axpby does not match managed axpby" );
}


void PetscVectorTests::VerifySwapPetscVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorc( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectord( d_factory->getNativeVector() );
    Vec veca, vecb;

    veca = vectora->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecb = vectorb->castTo<AMP::LinearAlgebra::PetscVector>().getVec();

    vectora->setRandomValues();
    vectorb->setToScalar( 99. );
    vectorc->copyVector( vectora );
    vectord->copyVector( vectorb );
    checkPetscError( utils, VecSwap( veca, vecb ) );
    vectorc->subtract( vectorc, vectorb );
    vectord->subtract( vectord, vectora );
    if ( ( vectorc->L1Norm() < 0.0000001 ) && ( vectord->L1Norm() < 0.000001 ) )
        utils->passes( "Swap vectors native interface works with native vectors" );
    else
        utils->failure( "Swap vectors native interface fails with native vectors" );
    vectorc->copyVector( vectora );
    vectord->copyVector( vectorb );
    vectora->swapVectors( vectorb );
    vectorc->subtract( vectorc, vectorb );
    vectord->subtract( vectord, vectora );
    if ( ( vectorc->L1Norm() < 0.0000001 ) && ( vectord->L1Norm() < 0.000001 ) )
        utils->passes( "Swap vectors AMP interface works with native vectors" );
    else
        utils->failure( "Swap vectors AMP interface fails with native vectors" );

    AMP::LinearAlgebra::Vector::shared_ptr vectore = d_factory->getManagedVector();
    AMP::LinearAlgebra::Vector::shared_ptr vectorf = d_factory->getManagedVector();
    AMP::LinearAlgebra::Vector::shared_ptr vectorg = d_factory->getManagedVector();
    AMP::LinearAlgebra::Vector::shared_ptr vectorh = d_factory->getManagedVector();

    Vec vece, vecf;
    vece = vectore->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecf = vectorf->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vectore->setRandomValues();
    vectorf->setToScalar( 99. );
    vectorg->copyVector( vectore );
    vectorh->copyVector( vectorf );
    checkPetscError( utils, VecSwap( vece, vecf ) );
    vectorg->subtract( vectorg, vectorf );
    vectorh->subtract( vectorh, vectore );
    if ( ( vectorg->L1Norm() < 0.0000001 ) && ( vectorh->L1Norm() < 0.000001 ) )
        utils->passes( "Swap vectors native interface works with managed vectors" );
    else
        utils->failure( "Swap vectors native interface fails with managed vectors" );

    vectorg->copyVector( vectore );
    vectorh->subtract( vectore, vectorg );
    vectorh->copyVector( vectorf );
    vectore->swapVectors( vectorf );
    vectorg->subtract( vectorg, vectorf );
    vectorh->subtract( vectorh, vectore );
    if ( ( vectorg->L1Norm() < 0.0000001 ) && ( vectorh->L1Norm() < 0.000001 ) )
        utils->passes( "Swap vectors managed interface works with managed vectors" );
    else
        utils->failure( "Swap vectors managed interface fails with managed vectors" );
}


void PetscVectorTests::VerifyGetSizePetscVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb = d_factory->getManagedVector();

    Vec veca, vecb;
    veca = vectora->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecb = vectorb->castTo<AMP::LinearAlgebra::PetscVector>().getVec();

    int sizea1, sizea2, sizeb1, sizeb2;
    sizea1 = vectora->getGlobalSize();
    checkPetscError( utils, VecGetSize( veca, &sizea2 ) );
    sizeb1 = vectorb->getGlobalSize();
    checkPetscError( utils, VecGetSize( vecb, &sizeb2 ) );

    if ( sizea1 == sizea2 )
        utils->passes( "Native PETSc: Native interface matches AMP interface" );
    else
        utils->failure( "Native PETSc: Native interface does not match AMP interface" );
    if ( sizeb1 == sizeb2 )
        utils->passes( "Managed PETSc: Native interface matches AMP interface" );
    else
        utils->failure( "Managed PETSc: Native interface does not match AMP interface" );
    if ( sizea1 == sizeb1 )
        utils->passes( "Managed PETSc matches native PETSc" );
}


void PetscVectorTests::VerifyMaxPointwiseDividePetscVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getNativeVector() );

    vectora->setRandomValues();
    vectorb->setToScalar( 3.14159 );

    AMP::LinearAlgebra::Vector::shared_ptr vectorc = d_factory->getManagedVector();
    vectorc->copyVector( vectora );
    AMP::LinearAlgebra::Vector::shared_ptr vectord = d_factory->getManagedVector();
    vectord->copyVector( vectorb );

    Vec veca, vecb, vecc, vecd;
    veca = vectora->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecb = vectorb->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecc = vectorc->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecd = vectord->castTo<AMP::LinearAlgebra::PetscVector>().getVec();

    double ans1, ans2;

    checkPetscError( utils, VecMaxPointwiseDivide( veca, vecb, &ans1 ) );
    checkPetscError( utils, VecMaxPointwiseDivide( vecc, vecd, &ans2 ) );

    if ( fabs( ans1 - ans2 ) < 0.000001 )
        utils->passes( "VecMaxPointwiseDivide working for both vector types" );
    else
        utils->failure( "VecMaxPointwiseDivide not working" );
}


void PetscVectorTests::VerifyAbsPetscVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getNativeVector() );
    Vec veca, vecb;
    veca = vectora->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecb = vectorb->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    if ( !veca || !vecb )
        utils->failure( "PETSC abs create" );

    vectora->setRandomValues();
    vectora->addScalar( vectora, 1. );
    vectorb->copyVector( vectora );
    vectora->scale( -1. );
    checkPetscError( utils, VecAbs( veca ) );
    vectorb->subtract( vectora, vectorb );
    if ( vectorb->L1Norm() < 0.000001 )
        utils->passes( "native interface on native petsc abs works" );
    else
        utils->failure( "native interface on native petsc abs doesn't work" );

    AMP::LinearAlgebra::Vector::shared_ptr vectorc = d_factory->getManagedVector();
    AMP::LinearAlgebra::Vector::shared_ptr vectord = d_factory->getManagedVector();

    Vec vecc, vecd;
    vecc = vectorc->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecd = vectord->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    if ( !vecc || !vecd )
        utils->failure( "PETSC abs create" );

    vectorc->setRandomValues();
    vectorc->addScalar( vectorc, 1. );
    vectord->copyVector( vectorc );
    vectorc->scale( -1. );
    checkPetscError( utils, VecAbs( vecc ) );
    vectord->subtract( vectorc, vectord );
    if ( vectord->L1Norm() < 0.000001 )
        utils->passes( "managed interface on native petsc abs works" );
    else
        utils->failure( "managed interface on native petsc abs doesn't work" );
}


void PetscVectorTests::VerifyPointwiseMultPetscVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorc( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectord( d_factory->getNativeVector() );
    if ( !vectora || !vectorb || !vectorc || !vectord )
        utils->failure( "PointwiseMult create" );

    vectora->setRandomValues();
    vectorb->setToScalar( 4.567 );
    vectorc->copyVector( vectora );
    vectord->copyVector( vectorb );

    Vec veca, vecb;
    veca = vectora->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecb = vectorb->castTo<AMP::LinearAlgebra::PetscVector>().getVec();

    checkPetscError( utils, VecPointwiseMult( veca, veca, vecb ) );
    vectorc->multiply( vectorc, vectord );
    vectorc->subtract( vectorc, vectora );
    if ( vectorc->L1Norm() < 0.000001 )
        utils->passes( "managed interface for native vector" );
    else
        utils->failure( "managed interface for native vector" );


    AMP::LinearAlgebra::Vector::shared_ptr vectore = d_factory->getManagedVector();
    AMP::LinearAlgebra::Vector::shared_ptr vectorf = d_factory->getManagedVector();
    AMP::LinearAlgebra::Vector::shared_ptr vectorg = d_factory->getManagedVector();
    AMP::LinearAlgebra::Vector::shared_ptr vectorh = d_factory->getManagedVector();

    vectore->setRandomValues();
    vectorf->setToScalar( 4.567 );
    vectorg->copyVector( vectore );
    vectorh->copyVector( vectorf );

    Vec vece, vecf;
    vece = vectore->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecf = vectorf->castTo<AMP::LinearAlgebra::PetscVector>().getVec();

    checkPetscError( utils, VecPointwiseMult( vece, vece, vecf ) );
    vectorg->multiply( vectorg, vectorh );
    vectorg->subtract( vectorg, vectore );

    if ( vectorg->L1Norm() < 0.000001 )
        utils->passes( "native interface for managed vector" );
    else
        utils->failure( "native interface for managed vector" );
}


void PetscVectorTests::VerifyPointwiseDividePetscVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorc( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectord( d_factory->getNativeVector() );
    vectora->setRandomValues();
    vectorb->setToScalar( 4.567 );
    vectorc->copyVector( vectora );
    vectord->copyVector( vectorb );

    Vec veca, vecb;
    veca = vectora->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecb = vectorb->castTo<AMP::LinearAlgebra::PetscVector>().getVec();

    checkPetscError( utils, VecPointwiseDivide( veca, veca, vecb ) );
    vectorc->divide( vectorc, vectord );
    vectorc->subtract( vectorc, vectora );
    if ( vectorc->L1Norm() < 0.000001 )
        utils->passes( "managed interface for native vector" );
    else
        utils->failure( "managed interface for native vector" );

    AMP::LinearAlgebra::Vector::shared_ptr vectore = d_factory->getManagedVector();
    AMP::LinearAlgebra::Vector::shared_ptr vectorf = d_factory->getManagedVector();
    AMP::LinearAlgebra::Vector::shared_ptr vectorg = d_factory->getManagedVector();
    AMP::LinearAlgebra::Vector::shared_ptr vectorh = d_factory->getManagedVector();

    vectore->setRandomValues();
    vectorf->setToScalar( 4.567 );
    vectorg->copyVector( vectore );
    vectorh->copyVector( vectorf );

    Vec vece, vecf;
    vece = vectore->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecf = vectorf->castTo<AMP::LinearAlgebra::PetscVector>().getVec();

    checkPetscError( utils, VecPointwiseDivide( vece, vece, vecf ) );
    vectorg->divide( vectorg, vectorh );
    vectorg->subtract( vectorg, vectore );

    if ( vectorg->L1Norm() < 0.000001 )
        utils->passes( "native interface for managed vector" );
    else
        utils->failure( "native interface for managed vector" );
}


void PetscVectorTests::VerifySqrtPetscVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getNativeVector() );

    vectora->setRandomValues();

    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getManagedVector() );
    vectorb->copyVector( vectora );

    Vec veca, vecb;
    veca = vectora->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecb = vectorb->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
    checkPetscError( utils, VecSqrt( veca ) );
    checkPetscError( utils, VecSqrt( vecb ) );
#elif PETSC_VERSION_GE(3,2,0)
    checkPetscError( utils, VecSqrtAbs( veca ) );
    checkPetscError( utils, VecSqrtAbs( vecb ) );
#else
#error Not programmed for this version yet
#endif
    bool equal = vectora->equals( vectorb );
    if ( equal )
        utils->passes( "Vector square root passes" );
    else
        utils->failure( "Vector square root fails" );
}


void PetscVectorTests::VerifySetRandomPetscVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getManagedVector() );

    Vec veca, vecb;
    veca = vectora->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecb = vectorb->castTo<AMP::LinearAlgebra::PetscVector>().getVec();

    vectora->setToScalar( 10.0 );
    vectorb->setToScalar( 10.0 );

    checkPetscError( utils, VecSetRandom( veca, PETSC_NULL ) );
    checkPetscError( utils, VecSetRandom( vecb, PETSC_NULL ) );

    if ( vectora->maxNorm() < 1. )
        utils->passes( "VecSetRandom passes for native petsc" );
    else
        utils->failure( "VecSetRandom fails for native petsc" );
    if ( vectorb->maxNorm() < 1. )
        utils->passes( "VecSetRandom passes for managed petsc" );
    else
        utils->failure( "VecSet fails for managed petsc" );

    vectora->setToScalar( 5.0 );
    vectorb->setToScalar( 6.0 );
    vectora->setRandomValues();
    vectorb->setRandomValues();

    if ( vectora->maxNorm() < 1. )
        utils->passes( "setToScalar passes for native petsc" );
    else
        utils->failure( "setToScalar fails for native petsc" );

    if ( vectorb->maxNorm() < 1. )
        utils->passes( "setToScalar passes for managed petsc" );
    else
        utils->failure( "setToScalar fails for managed petsc" );
}


void PetscVectorTests::VerifySetPetscVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getManagedVector() );

    Vec veca, vecb;
    veca = vectora->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecb = vectorb->castTo<AMP::LinearAlgebra::PetscVector>().getVec();

    checkPetscError( utils, VecSet( veca, 2.0 ) );
    checkPetscError( utils, VecSet( vecb, 3.0 ) );

    if ( ( fabs( vectora->L1Norm() - vectora->getGlobalSize() * 2.0 ) < 0.000001 ) &&
         ( fabs( vectora->maxNorm() - 2.0 ) < 0.000001 ) )
        utils->passes( "VecSet passes for native petsc" );
    else
        utils->failure( "VecSet fails for native petsc" );
    if ( ( fabs( vectorb->L1Norm() - vectorb->getGlobalSize() * 3.0 ) < 0.000001 ) &&
         ( fabs( vectorb->maxNorm() - 3.0 ) < 0.000001 ) )
        utils->passes( "VecSet passes for managed petsc" );
    else
        utils->failure( "VecSet fails for managed petsc" );
    vectora->setToScalar( 5.0 );
    vectorb->setToScalar( 6.0 );
    if ( ( fabs( vectora->L1Norm() - vectora->getGlobalSize() * 5.0 ) < 0.000001 ) &&
         ( fabs( vectora->maxNorm() - 5.0 ) < 0.000001 ) )
        utils->passes( "setToScalar passes for native petsc" );
    else
        utils->failure( "setToScalar fails for native petsc" );
    if ( ( fabs( vectorb->L1Norm() - vectorb->getGlobalSize() * 6.0 ) < 0.000001 ) &&
         ( fabs( vectorb->maxNorm() - 6.0 ) < 0.000001 ) )
        utils->passes( "setToScalar passes for managed petsc" );
    else
        utils->failure( "setToScalar fails for managed petsc" );
}


void PetscVectorTests::VerifyAXPYPetscVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectora2( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectora_orig( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb2( d_factory->getNativeVector() );
    Vec veca, vecb, veca2, veca_orig;
    veca      = vectora->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    veca2     = vectora2->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    veca_orig = vectora_orig->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecb      = vectorb->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    if ( !veca || !vecb || !veca || !veca_orig )
        utils->failure( "PETSc AXPY create" );

    vectora->setRandomValues();
    vectorb->setRandomValues();
    vectora_orig->copyVector( vectora );
    vectora2->copyVector( vectora );
    vectorb2->copyVector( vectorb );
    checkPetscError( utils, VecAXPY( veca, 1.23456, vecb ) );
    vectora2->axpy( 1.23456, vectorb2, vectora2 );
#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
    PetscTruth ans;
#elif PETSC_VERSION_GE(3,2,0)
    PetscBool ans;
#else
#error Not programmed for this version yet
#endif
    checkPetscError( utils, VecEqual( veca, veca2, &ans ) );
    if ( ans == PETSC_TRUE )
        utils->passes( "native interface on native petsc axpy works" );
    else
        utils->failure( "native interface on native petsc axpy doesn't work" );

    AMP::LinearAlgebra::Vector::shared_ptr vectorc  = d_factory->getManagedVector();
    AMP::LinearAlgebra::Vector::shared_ptr vectorc2 = d_factory->getManagedVector();
    AMP::LinearAlgebra::Vector::shared_ptr vectord  = d_factory->getManagedVector();
    AMP::LinearAlgebra::Vector::shared_ptr vectord2 = d_factory->getManagedVector();

    vectorc->copyVector( vectora_orig );
    vectorc2->copyVector( vectora_orig );
    vectord->copyVector( vectorb );
    vectord2->copyVector( vectorb2 );

    Vec vecc, vecd, vecc2, vecd2;
    vecc  = vectorc->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecc2 = vectorc2->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecd  = vectord->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecd2 = vectord2->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    checkPetscError( utils, VecAXPY( vecc, 1.23456, vecd ) );
    vectorc2->axpy( 1.23456, vectord2, vectorc2 );
    if ( !vecc || !vecd || !vecc2 || !vecd2 )
        utils->failure( "PETSC AXPY create" );

    if ( fabs( vectorc->L1Norm() - vectorc2->L1Norm() ) < 0.000001 )
        utils->passes( "managed interface passes l1 norm test of axpy" );
    else
        utils->failure( "managed interface fails l1 norm test of axpy" );
    if ( vectorc->L2Norm() == vectorc2->L2Norm() )
        utils->passes( "managed interface passes l2 norm test of axpy" );
    else
        utils->failure( "managed interface fails l2 norm test of axpy" );
    if ( vectorc->maxNorm() == vectorc2->maxNorm() )
        utils->passes( "managed interface passes inf norm test of axpy" );
    else
        utils->failure( "managed interface fails inf norm test of axpy" );

    if ( fabs( vectorc->L1Norm() - vectora->L1Norm() ) < 0.000001 )
        utils->passes( "managed and native L1 norms the same" );
    else
        utils->failure( "managed and native L1 norms different" );
    if ( fabs( vectorc->L2Norm() - vectora->L2Norm() ) < 0.000001 )
        utils->passes( "managed and native L2 norms the same" );
    else
        utils->failure( "managed and native L2 norms different" );
    if ( fabs( vectorc->maxNorm() - vectora->maxNorm() ) < 0.000001 )
        utils->passes( "managed and native inf norms the same" );
    else
        utils->failure( "managed and native inf norms different" );
}


void PetscVectorTests::VerifyScalePetscVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectora2( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getNativeVector() );
    Vec veca, vecb, veca2;
    veca  = vectora->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    veca2 = vectora2->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecb  = vectorb->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    if ( !veca || !veca2 || !vecb )
        utils->failure( "PETSc scale create" );

    vectora->setRandomValues();
    vectorb->copyVector( vectora );
    vectora2->copyVector( vectora );
    double norm1, norm2;
    checkPetscError( utils, VecScale( veca, 1.23456 ) );
    norm1 = vectora->L2Norm();
    norm2 = 1.23456 * vectorb->L2Norm();
    if ( fabs( norm1 - norm2 ) < 0.000001 )
        utils->passes( "native interface on native petsc scaling works" );
    else
        utils->failure( "native interface on native petsc scaling doesn't work" );
    vectora->scale( 1. / 1.23456 );
    if ( fabs( vectora->L2Norm() - vectorb->L2Norm() ) < 0.000001 )
        utils->passes( "AMP interface on native petsc scaling works" );
    else
        utils->failure( "AMP interface on native petsc scaling doesn't work" );
    checkPetscError( utils, VecScale( veca2, 1.234567 ) );
    checkPetscError( utils, VecScale( veca2, 99.99 ) );
    if ( fabs( vectora2->L2Norm() - 99.99 * 1.234567 * vectorb->L2Norm() ) < 0.000001 )
        utils->passes( "Multiple scales working in native petsc" );
    else
        utils->failure( "Multiple scales failing in native petsc" );

    AMP::LinearAlgebra::Vector::shared_ptr vectorc = d_factory->getManagedVector();
    AMP::LinearAlgebra::Vector::shared_ptr vectord = d_factory->getManagedVector();
    vectorc->copyVector( vectora );
    vectord->copyVector( vectora );

    Vec vecc;
    vecc = vectorc->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    double norm3, norm4;
    checkPetscError( utils, VecScale( vecc, 1.23456 ) );
    norm3 = vectorc->L2Norm();
    norm4 = 1.23456 * vectord->L2Norm();
    if ( fabs( norm3 - norm4 ) < 0.000001 )
        utils->passes( "native interface on managed petsc scaling works" );
    else
        utils->failure( "native interface on managed petsc scaling doesn't work" );
    vectorc->scale( 1. / 1.23456 );
    if ( fabs( vectorc->L2Norm() - vectord->L2Norm() ) < 0.000001 )
        utils->passes( "AMP interface on managed petsc scaling works" );
    else
        utils->failure( "AMP interface on managed petsc scaling doesn't work" );
}


void PetscVectorTests::VerifyDotPetscVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr vectora( d_factory->getNativeVector() );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb( d_factory->getNativeVector() );
    Vec veca, vecb;
    veca = vectora->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecb = vectorb->castTo<AMP::LinearAlgebra::PetscVector>().getVec();

    vectora->setRandomValues();
    vectorb->setRandomValues();
    double dot1, dot2, dot12;
    checkPetscError( utils, VecDot( veca, vecb, &dot1 ) );
    checkPetscError( utils, VecDot( veca, vecb, &dot12 ) );
    dot2 = vectora->dot( vectorb );
    if ( dot1 == dot2 ) // These should be identical, since same method called
        utils->passes( "native dot equals interface dot for native vector" );
    else
        utils->failure( "native dot does not equal interface dot for native vector" );
    if ( dot1 == dot12 )
        utils->passes( "multiple native dot passes" );
    else
        utils->failure( "multiple native dot fails" );

    AMP::LinearAlgebra::Vector::shared_ptr vectorc = d_factory->getManagedVector();
    AMP::LinearAlgebra::Vector::shared_ptr vectord = d_factory->getManagedVector();
    vectorc->copyVector( vectora );
    vectord->copyVector( vectorb );
    Vec vecc, vecd;
    vecc = vectorc->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecd = vectord->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    double dot3, dot4;
    checkPetscError( utils, VecDot( vecc, vecd, &dot3 ) );
    dot4 = vectorc->dot( vectord );
    if ( dot3 == dot4 ) // These should be identical, since same method called
        utils->passes( "native dot equals interface dot for managed vector" );
    else
        utils->failure( "native dot does not equal interface dot for managed vector" );
    if ( fabs( dot3 - dot1 ) < 0.00000001 ) // This may test two different implementations
        utils->passes( "native dot equals managed dot" );
    else
        utils->failure( "native dot does not equal managed dot" );

    AMP::LinearAlgebra::Vector::shared_ptr vectore = d_factory->getManagedVector();
    AMP::LinearAlgebra::Vector::shared_ptr vectorf = d_factory->getManagedVector();
    vectore->copyVector( vectora );
    vectorf->copyVector( vectorb );
    Vec vece, vecf;
    vece = vectore->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    vecf = vectorf->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    double dot5, dot6;
    checkPetscError( utils, VecDot( vece, vecf, &dot5 ) );
    dot6 = vectore->dot( vectorf );
    if ( dot5 == dot6 ) // These should be identical, since same method called
        utils->passes( "native dot equals interface dot for managed alloc vector" );
    else
        utils->failure( "native dot does not equal interface dot for managed alloc vector" );
    if ( dot3 == dot5 ) // Again, same function
        utils->passes( "native alloc dot equals managed alloc dot" );
    else
        utils->failure( "native alloc dot does not equal managed alloc dot" );

    /**** Need to test for failures at some point...
    try {
      checkPetscError ( utils , VecDot ( veca , vecc , &dot1 ) );  // This
    should fail
      utils->failure ( "incorrect failure of checkPetscError ( utils , VecDot" )
    );
    } catch ( ... ) {
      utils->passes ( "correct failure of checkPetscError ( utils , VecDot" ) );
    }
    ***/
}



}
}

/// \endcond

#endif
