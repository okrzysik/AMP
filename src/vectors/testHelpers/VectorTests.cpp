#include "AMP/vectors/testHelpers/VectorTests.h"

#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/Vector.h"
#include <algorithm>
#ifdef USE_EXT_SUNDIALS
#include "AMP/vectors/sundials/ManagedSundialsVector.h"
#include "AMP/vectors/sundials/SundialsVector.h"
#endif
#ifdef USE_EXT_PETSC
#include "AMP/vectors/petsc/ManagedPetscVector.h"
#include "AMP/vectors/petsc/PetscVector.h"
#endif


namespace AMP {
namespace LinearAlgebra {


static inline int lround( double x ) { return x >= 0 ? floor( x ) : ceil( x ); }

void VectorTests::InstantiateVector( AMP::UnitTest *utils )
{
    auto vector = d_factory->getVector();
    if ( vector )
        utils->passes( "created " + d_factory->name() );
    else
        utils->failure( "created " + d_factory->name() );
}


void VectorTests::CopyVectorConsistency( AMP::UnitTest *utils )
{
    auto vec1        = d_factory->getVector();
    auto vec2        = vec1->cloneVector();
    auto vec3        = vec1->cloneVector();
    auto commList    = vec1->getCommunicationList();
    double *t1       = nullptr;
    double *t2       = nullptr;
    size_t *ndx      = nullptr;
    size_t numGhosts = commList->getGhostIDList().size();
    vec1->setRandomValues();
    vec2->copyVector( vec1 );
    if ( numGhosts ) {
        t1  = new double[numGhosts];
        t2  = new double[numGhosts];
        ndx = new size_t[numGhosts];
        std::copy( commList->getGhostIDList().begin(), commList->getGhostIDList().end(), ndx );
        vec1->getValuesByGlobalID( numGhosts, ndx, t1 );
        vec2->getValuesByGlobalID( numGhosts, ndx, t2 );
        if ( std::equal( t1, t1 + numGhosts, t2 ) )
            utils->passes( "Ghosts are the same (1) - " + d_factory->name() );
        else
            utils->failure( "Ghosts are different (1) - " + d_factory->name() );
    }

    vec1->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
    vec3->copyVector( vec1 );
    if ( numGhosts ) {
        vec1->getValuesByGlobalID( numGhosts, ndx, t1 );
        vec3->getValuesByGlobalID( numGhosts, ndx, t2 );
        if ( std::equal( t1, t1 + numGhosts, t2 ) )
            utils->passes( "Ghosts are the same (2) - " + d_factory->name() );
        else
            utils->failure( "Ghosts are different (2) - " + d_factory->name() );
        delete[] t1;
        delete[] t2;
        delete[] ndx;
    }
}


void VectorTests::Bug_728( AMP::UnitTest *utils )
{
    auto vector = d_factory->getVector();
    auto var1   = vector->getVariable();
    if ( !var1 )
        return;
    auto var2 = var1->cloneVariable( var1->getName() );
    if ( vector->subsetVectorForVariable( var1 ) )
        utils->passes( "Found vector for same variable pointer " + d_factory->name() );
    else
        utils->failure( "Did not find vector for same variable pointer " + d_factory->name() );

    if ( vector->subsetVectorForVariable( var2 ) )
        utils->passes( "Found vector for cloned variable pointer " + d_factory->name() );
    else
        utils->failure( "Did not find vector for cloned variable pointer " + d_factory->name() );
}


void VectorTests::SetToScalarVector( AMP::UnitTest *utils )
{
    auto vector = d_factory->getVector();
    vector->setToScalar( 0. );
    utils->passes( "setToScalar ran to completion " + d_factory->name() );
    bool fail   = false;
    auto curVec = vector->begin();
    auto endVec = vector->end();
    while ( curVec != endVec ) {
        if ( *curVec != 0. ) {
            fail = true;
            break;
        }
        ++curVec;
    }
    if ( !fail )
        utils->passes( "Set data to 0 " + d_factory->name() );
    else
        utils->failure( "Failed to set scalar to 0 " + d_factory->name() );
    fail = false;
    vector->setToScalar( 5. );
    auto curVal = vector->begin();
    while ( curVal != endVec ) {
        if ( *curVal != 5. ) {
            fail = true;
            break;
        }
        ++curVal;
    }
    if ( !fail )
        utils->passes( "Set data to 5 " + d_factory->name() );
    else
        utils->failure( "Failed to set scalar to 5 " + d_factory->name() );
    auto remoteDofs = vector->getDOFManager()->getRemoteDOFs();
    fail            = false;
    for ( auto &remoteDof : remoteDofs ) {
        if ( vector->getValueByGlobalID( remoteDof ) != 5. )
            fail = true;
    }
    if ( !fail )
        utils->passes( "Set ghost data to 5 " + d_factory->name() );
    else
        utils->failure( "Failed to set ghost scalar values to 5 " + d_factory->name() );
}


void VectorTests::CloneVector( AMP::UnitTest *utils )
{
    auto vector = d_factory->getVector();
    auto clone  = vector->cloneVector( "cloned vector" );
    clone->setToScalar( 0. );
    utils->passes( "Clone created " + d_factory->name() );
    bool pass = true;
    for ( size_t i = 0; i != vector->numberOfDataBlocks(); i++ ) {
        auto *clone_ptr  = clone->getRawDataBlock<double>( i );
        auto *vector_ptr = vector->getRawDataBlock<double>( i );
        if ( clone_ptr == vector_ptr )
            pass = false;
    }
    if ( pass )
        utils->passes( "New data allocated " + d_factory->name() );
    else
        utils->failure( "Data not allocated " + d_factory->name() );
}


void VectorTests::DotProductVector( AMP::UnitTest *utils )
{
    auto vector1( d_factory->getVector() );
    auto vector2( d_factory->getVector() );
    vector1->setToScalar( 1. );
    vector2->setToScalar( 2. );
    auto d11 = vector1->dot( vector1 );
    auto d12 = vector1->dot( vector2 );
    auto d21 = vector2->dot( vector1 );
    auto d22 = vector2->dot( vector2 );
    if ( 2. * d11 == d12 )
        utils->passes( "dot product 1 " + d_factory->name() );
    else
        utils->failure( "dot product 1 " + d_factory->name() );
    if ( 2. * d11 == d21 )
        utils->passes( "dot product 2 " + d_factory->name() );
    else
        utils->failure( "dot product 2 " + d_factory->name() );
    if ( 4. * d11 == d22 )
        utils->passes( "dot product 3 " + d_factory->name() );
    else
        utils->failure( "dot product 3 " + d_factory->name() );
    if ( d11 == vector1->getGlobalSize() )
        utils->passes( "dot product 4 " + d_factory->name() );
    else
        utils->failure( "dot product 4 " + d_factory->name() );
    if ( d21 == d12 )
        utils->passes( "dot product 5 " + d_factory->name() );
    else
        utils->failure( "dot product 5 " + d_factory->name() );
}


void VectorTests::L2NormVector( AMP::UnitTest *utils )
{
    auto vector = d_factory->getVector();
    vector->setToScalar( 1. );
    auto norm  = vector->L2Norm();
    auto norm2 = vector->dot( vector );
    if ( fabs( norm * norm - norm2 ) < 0.000001 )
        utils->passes( "L2 norm 1 " + d_factory->name() );
    else
        utils->failure( "L2 norm 1 " + d_factory->name() );
    vector->setRandomValues();
    norm  = vector->L2Norm();
    norm2 = vector->dot( vector );
    if ( fabs( norm * norm - norm2 ) < 0.000001 )
        utils->passes( "L2 norm 2 " + d_factory->name() );
    else
        utils->failure( "L2 norm 2 " + d_factory->name() );
}


void VectorTests::AbsVector( AMP::UnitTest *utils )
{
    auto vec1 = d_factory->getVector();
    auto vec2 = vec1->cloneVector();
    vec1->setRandomValues();
    vec2->copyVector( vec1 );
    vec2->scale( -1.0 );
    vec2->abs( vec2 );
    if ( vec1->equals( vec2 ) )
        utils->passes( "Abs passes " + d_factory->name() );
    else
        utils->failure( "Abs fails " + d_factory->name() );
}


void VectorTests::L1NormVector( AMP::UnitTest *utils )
{
    auto vector = d_factory->getVector();
    auto vector_1( d_factory->getVector() );
    vector->setRandomValues();
    vector_1->setToScalar( 1. );
    auto norm = vector->L1Norm();
    vector->abs( vector );
    auto norm2 = vector->dot( vector_1 );
    if ( fabs( norm - norm2 ) < 0.000001 )
        utils->passes( "L1 norm " + d_factory->name() );
    else
        utils->failure( "L1 norm " + d_factory->name() );
}


void VectorTests::MaxNormVector( AMP::UnitTest *utils )
{
    auto vector = d_factory->getVector();
    vector->setRandomValues();
    auto infNorm = vector->maxNorm();
    vector->abs( vector );
    auto curData     = vector->begin();
    auto endData     = vector->end();
    double local_ans = *curData;
    while ( curData != endData ) {
        local_ans = std::max( local_ans, *curData );
        ++curData;
    }
    double global_ans = vector->getComm().maxReduce( local_ans );
    if ( global_ans == infNorm )
        utils->passes( "Inf norm " + d_factory->name() );
    else
        utils->failure( "Inf norm " + d_factory->name() );
}


void VectorTests::ScaleVector( AMP::UnitTest *utils )
{
    auto vector1( d_factory->getVector() );
    auto vector2( d_factory->getVector() );
    double beta = 1.2345;
    vector2->setRandomValues();
    vector1->scale( beta, vector2 );
    bool pass     = true;
    auto curData1 = vector1->begin();
    auto endData1 = vector1->end();
    auto curData2 = vector2->begin();
    while ( curData1 != endData1 ) {
        if ( *curData1 != beta * *curData2 )
            pass = false;
        ++curData1;
        ++curData2;
    }
    if ( pass )
        utils->passes( "scale vector 1 " + d_factory->name() );
    else
        utils->failure( "scale vector 1 " + d_factory->name() );
    vector2->scale( beta );
    vector1->subtract( vector2, vector1 );
    if ( vector1->maxNorm() < 0.0000001 )
        utils->passes( "scale vector 2 " + d_factory->name() );
    else
        utils->failure( "scale vector 2 " + d_factory->name() );
}


#ifdef USE_EXT_PETSC
void VectorTests::Bug_491( AMP::UnitTest *utils )
{
    auto vector1( d_factory->getVector() );
    vector1->setRandomValues();
    auto managed_petsc = AMP::LinearAlgebra::PetscVector::view( vector1 );
    auto petsc_vec  = std::dynamic_pointer_cast<AMP::LinearAlgebra::PetscVector>( managed_petsc );
    Vec managed_vec = petsc_vec->getVec();


    // This sets the petsc cache
    double n1, n2, ninf;
    VecNormBegin( managed_vec, NORM_1, &n1 );
    VecNormBegin( managed_vec, NORM_2, &n2 );
    VecNormBegin( managed_vec, NORM_INFINITY, &ninf );
    VecNormEnd( managed_vec, NORM_1, &n1 );
    VecNormEnd( managed_vec, NORM_2, &n2 );
    VecNormEnd( managed_vec, NORM_INFINITY, &ninf );
    VecNorm( managed_vec, NORM_1, &n1 );
    VecNorm( managed_vec, NORM_2, &n2 );
    VecNorm( managed_vec, NORM_INFINITY, &ninf );

    // Now, we perform some math on vector1
    vector1->scale( 100000 );
    double sp_n1  = vector1->L1Norm();
    double sp_n2  = vector1->L2Norm();
    double sp_inf = vector1->maxNorm();

    // Check to see if petsc cache has been invalidated
    VecNormBegin( managed_vec, NORM_1, &n1 );
    VecNormBegin( managed_vec, NORM_2, &n2 );
    VecNormBegin( managed_vec, NORM_INFINITY, &ninf );
    VecNormEnd( managed_vec, NORM_1, &n1 );
    VecNormEnd( managed_vec, NORM_2, &n2 );
    VecNormEnd( managed_vec, NORM_INFINITY, &ninf );

    if ( fabs( n1 - sp_n1 ) < 0.00000001 * n1 )
        utils->passes( "L1 norm -- Petsc interface begin/end " + d_factory->name() );
    else
        utils->failure( "l1 norm -- Petsc interface begin/end " + d_factory->name() );
    if ( fabs( n2 - sp_n2 ) < 0.00000001 * n1 )
        utils->passes( "L2 norm -- Petsc interface begin/end " + d_factory->name() );
    else
        utils->failure( "l2 norm -- Petsc interface begin/end " + d_factory->name() );
    if ( fabs( ninf - sp_inf ) < 0.00000001 * n1 )
        utils->passes( "Linf norm -- Petsc interface begin/end " + d_factory->name() );
    else
        utils->failure( "Linf norm -- Petsc interface begin/end " + d_factory->name() );

    VecNorm( managed_vec, NORM_1, &n1 );
    VecNorm( managed_vec, NORM_2, &n2 );
    VecNorm( managed_vec, NORM_INFINITY, &ninf );

    if ( fabs( n1 - vector1->L1Norm() ) < 0.00000001 * n1 )
        utils->passes( "L1 norm -- Petsc interface begin/end " + d_factory->name() );
    else
        utils->failure( "l1 norm -- Petsc interface begin/end " + d_factory->name() );
    if ( fabs( n2 - vector1->L2Norm() ) < 0.00000001 * n1 )
        utils->passes( "L2 norm -- Petsc interface begin/end " + d_factory->name() );
    else
        utils->failure( "l2 norm -- Petsc interface begin/end " + d_factory->name() );
    if ( fabs( ninf - vector1->maxNorm() ) < 0.00000001 * n1 )
        utils->passes( "inf norm -- Petsc interface begin/end " + d_factory->name() );
    else
        utils->failure( "inf norm -- Petsc interface begin/end " + d_factory->name() );
}
#endif


void VectorTests::AddVector( AMP::UnitTest *utils )
{
    auto vector1( d_factory->getVector() );
    auto vector2( d_factory->getVector() );
    auto vector3( d_factory->getVector() );
    vector1->setRandomValues();
    vector2->setRandomValues();
    vector3->add( vector1, vector2 );
    bool pass     = true;
    auto curData1 = vector1->begin();
    auto endData1 = vector1->end();
    auto curData2 = vector2->begin();
    auto curData3 = vector3->begin();
    while ( curData1 != endData1 ) {
        if ( *curData3 != *curData1 + *curData2 )
            pass = false;
        ++curData1;
        ++curData2;
        ++curData3;
    }

    if ( pass )
        utils->passes( "add vector " + d_factory->name() );
    else
        utils->failure( "add vector " + d_factory->name() );
}


void VectorTests::SubtractVector( AMP::UnitTest *utils )
{
    auto vector1( d_factory->getVector() );
    auto vector2( d_factory->getVector() );
    auto vector3( d_factory->getVector() );
    auto vector4( d_factory->getVector() );
    vector1->setRandomValues();
    vector2->setRandomValues();
    vector3->subtract( vector1, vector2 );
    bool pass     = true;
    auto curData1 = vector1->begin();
    auto endData1 = vector1->end();
    auto curData2 = vector2->begin();
    auto curData3 = vector3->begin();
    while ( curData1 != endData1 ) {
        if ( *curData3 != *curData1 - *curData2 )
            pass = false;
        ++curData1;
        ++curData2;
        ++curData3;
    }
    if ( pass )
        utils->passes( "vector subtract 1 " + d_factory->name() );
    else
        utils->failure( "vector subtract 1 " + d_factory->name() );
    vector2->scale( -1. );
    vector4->add( vector1, vector2 );
    vector4->subtract( vector3, vector4 );
    if ( vector4->maxNorm() < 0.0000001 )
        utils->passes( "vector subtract 2 " + d_factory->name() );
    else
        utils->failure( "vector subtract 2 " + d_factory->name() );
}


void VectorTests::MultiplyVector( AMP::UnitTest *utils )
{
    auto vector1( d_factory->getVector() );
    auto vector2( d_factory->getVector() );
    auto vector3( d_factory->getVector() );
    vector1->setRandomValues();
    vector2->setToScalar( 3. );
    vector3->multiply( vector1, vector2 );
    bool pass     = true;
    auto curData1 = vector1->begin();
    auto endData1 = vector1->end();
    auto curData2 = vector2->begin();
    auto curData3 = vector3->begin();
    while ( curData1 != endData1 ) {
        if ( *curData3 != *curData1 * *curData2 )
            pass = false;
        ++curData1;
        ++curData2;
        ++curData3;
    }
    if ( pass )
        utils->passes( "vector::multiply " + d_factory->name() );
    else
        utils->failure( "vector::multiply " + d_factory->name() );
}


void VectorTests::DivideVector( AMP::UnitTest *utils )
{
    auto vector1( d_factory->getVector() );
    auto vector2( d_factory->getVector() );
    auto vector3( d_factory->getVector() );
    vector1->setRandomValues();
    vector2->setRandomValues();
    vector3->divide( vector1, vector2 );
    bool pass    = true;
    auto curVal1 = vector1->begin();
    auto curVal2 = vector2->begin();
    auto curVal3 = vector3->begin();
    auto endVal3 = vector3->end();
    while ( curVal3 != endVal3 ) {
        if ( *curVal3 != *curVal1 / *curVal2 )
            pass = false;
        ++curVal1;
        ++curVal2;
        ++curVal3;
    }
    if ( pass )
        utils->passes( "vector::divide " + d_factory->name() );
    else
        utils->failure( "vector::divide " + d_factory->name() );
}


void VectorTests::VectorIteratorLengthTest( AMP::UnitTest *utils )
{
    auto vector1( d_factory->getVector() );
    size_t i = 0;
    for ( auto it = vector1->begin(); it != vector1->end(); ++it )
        i++;
    size_t k = vector1->getLocalSize();
    if ( i == k )
        utils->passes( "Iterated over the correct number of entries " + d_factory->name() );
    else
        utils->failure( "Wrong number of entries in iterator " + d_factory->name() );
}


void VectorTests::VectorIteratorTests( AMP::UnitTest *utils )
{
    auto vector1 = d_factory->getVector();
    both_VectorIteratorTests<AMP::LinearAlgebra::Vector::iterator>( vector1, utils );
    both_VectorIteratorTests<AMP::LinearAlgebra::Vector::const_iterator>( vector1, utils );
}


void VectorTests::VerifyVectorMin( AMP::UnitTest *utils )
{
    auto vec = d_factory->getVector();
    vec->setRandomValues();
    vec->scale( -1.0 ); // make negative
    if ( fabs( vec->min() + vec->maxNorm() ) < 1.e-10 )
        utils->passes( "minimum of negative vector == ||.||_infty " + d_factory->name() );
    else
        utils->failure( "minimum of negative vector != ||.||_infty " + d_factory->name() );
}


void VectorTests::VerifyVectorMax( AMP::UnitTest *utils )
{
    auto vec = d_factory->getVector();
    vec->setRandomValues();
    if ( fabs( vec->max() - vec->maxNorm() ) < 1.e-10 )
        utils->passes( "maximum of positive vector == ||.||_infty " + d_factory->name() );
    else
        utils->failure( "maximum of positive vector != ||.||_infty " + d_factory->name() );
}


void VectorTests::VerifyVectorMaxMin( AMP::UnitTest *utils )
{
    auto vec    = d_factory->getVector();
    bool passes = true;
    for ( size_t i = 0; i != 10; i++ ) {
        vec->setRandomValues();
        vec->addScalar( vec, -0.5 );
        vec->scale( 2.0 ); // vec i.i.d [-1,1);
        auto max = vec->max();
        auto min = vec->min();
        auto ans = std::max( fabs( max ), fabs( min ) );
        if ( fabs( ans - vec->maxNorm() ) >= 1.e-20 )
            passes = false;
    }
    if ( passes )
        utils->passes( "Max and min correctly predict maxNorm() " + d_factory->name() );
    else
        utils->failure( "Max and min fail to predict maxNorm() " + d_factory->name() );
}


static void SetRandomValuesVectorVerify( AMP::UnitTest *utils,
                                         AMP::LinearAlgebra::Vector::shared_ptr v )
{
    if ( v->min() >= 0 )
        utils->passes( "SetRandomValuesVector: Min value >= 0" );
    else
        utils->failure( "SetRandomValuesVector: Min value < 0" );
    if ( v->max() < 1 )
        utils->passes( "SetRandomValuesVector: Max value < 1" );
    else
        utils->failure( "SetRandomValuesVector: Max value >= 1" );
    if ( v->L2Norm() > 0 )
        utils->passes( "SetRandomValuesVector: Non-zero vector created" );
    else
        utils->failure( "SetRandomValuesVector: Zero vector created" );
}
void VectorTests::SetRandomValuesVector( AMP::UnitTest *utils )
{
    auto vector  = d_factory->getVector();
    auto l2norm1 = -1;
    for ( size_t i = 0; i < 5; i++ ) {
        vector->setRandomValues();
        auto l2norm2 = vector->L2Norm();
        if ( fabs( l2norm1 - l2norm2 ) > 0.000001 )
            utils->passes( "Distinct vector created " + d_factory->name() );
        else
            utils->failure( "Similar vector created " + d_factory->name() );
        l2norm1 = l2norm2;
        SetRandomValuesVectorVerify( utils, vector );
    }
}


void VectorTests::ReciprocalVector( AMP::UnitTest *utils )
{
    auto vectora = d_factory->getVector();
    auto vectorb = d_factory->getVector();
    auto vectorc = d_factory->getVector();
    auto vectord = d_factory->getVector();
    auto vector1 = d_factory->getVector();
    vectora->setRandomValues();
    vectorb->reciprocal( vectora );
    vector1->setToScalar( 1. );
    vectorc->divide( vector1, vectora );
    vectord->subtract( vectorb, vectorc );
    if ( vectord->maxNorm() < 0.0000001 ) {
        utils->passes( "vector::reciprocal " + d_factory->name() );
    } else {
        utils->failure( "vector::reciprocal " + d_factory->name() );
    }
}


static void LinearSumVectorRun( std::shared_ptr<const VectorFactory> factory,
                                AMP::UnitTest *utils,
                                double alpha,
                                double beta,
                                const char *msg )
{
    auto vectora = factory->getVector();
    auto vectorb = factory->getVector();
    auto vectorc = factory->getVector();
    auto vectord = factory->getVector();
    vectora->setRandomValues();
    vectorb->setRandomValues();
    vectorc->linearSum( alpha, vectora, beta, vectorb );
    vectora->scale( alpha );
    vectorb->scale( beta );
    vectord->add( vectora, vectorb );
    vectord->subtract( vectorc, vectord );
    if ( vectord->maxNorm() < 0.0000001 )
        utils->passes( msg );
    else
        utils->failure( msg );
}
void VectorTests::LinearSumVector( AMP::UnitTest *utils )
{
    LinearSumVectorRun( d_factory, utils, 1.2345, 9.8765, "linear sum 1" );
    LinearSumVectorRun( d_factory, utils, -1.2345, 9.8765, "linear sum 2" );
    LinearSumVectorRun( d_factory, utils, 1.2345, -9.8765, "linear sum 3" );
    LinearSumVectorRun( d_factory, utils, -1.2345, -9.8765, "linear sum 4" );
}


static void AxpyVectorRun( std::shared_ptr<const VectorFactory> factory,
                           AMP::UnitTest *utils,
                           double alpha,
                           const char *msg )
{
    auto vectora = factory->getVector();
    auto vectorb = factory->getVector();
    auto vectorc = factory->getVector();
    auto vectord = factory->getVector();
    vectora->setRandomValues();
    vectorb->setRandomValues();
    vectorc->linearSum( alpha, vectora, 1., vectorb );
    vectord->axpy( alpha, vectora, vectorb );
    vectorc->subtract( vectorc, vectord );
    if ( vectorc->maxNorm() < 0.0000001 )
        utils->passes( msg );
    else
        utils->failure( msg );
}
void VectorTests::AxpyVector( AMP::UnitTest *utils )
{
    AxpyVectorRun( d_factory, utils, 6.38295, "axpy 1" );
    AxpyVectorRun( d_factory, utils, -6.38295, "axpy 2" );
}


static void AxpbyVectorRun( std::shared_ptr<const VectorFactory> factory,
                            AMP::UnitTest *utils,
                            double alpha,
                            double beta,
                            const char *msg )
{
    auto vectora  = factory->getVector();
    auto vectorb  = factory->getVector();
    auto vectorb1 = factory->getVector();
    auto vectorc  = factory->getVector();
    auto vectord  = factory->getVector();
    vectora->setRandomValues();
    vectorb->setRandomValues();
    vectorc->copyVector( vectorb );
    vectorb1->linearSum( alpha, vectora, beta, vectorb );
    vectorb->linearSum( alpha, vectora, beta, vectorb );
    vectorc->axpby( alpha, beta, vectora );
    vectord->subtract( vectorc, vectorb1 );
    auto maxNorm = vectord->maxNorm();
    if ( maxNorm < 0.0000001 )
        utils->passes( msg );
    else
        utils->failure( msg );
    vectord->subtract( vectorc, vectorb );
    maxNorm = vectord->maxNorm();
    if ( maxNorm < 0.0000001 )
        utils->passes( msg );
    else
        utils->failure( msg );
}
void VectorTests::AxpbyVector( AMP::UnitTest *utils )
{
    AxpbyVectorRun( d_factory, utils, 6.38295, 99.273624, "axpby 1" );
    AxpbyVectorRun( d_factory, utils, 6.38295, -99.273624, "axpby 2" );
    AxpbyVectorRun( d_factory, utils, -6.38295, 99.273624, "axpby 3" );
    AxpbyVectorRun( d_factory, utils, -6.38295, -99.273624, "axpby 4" );
}


void VectorTests::CopyVector( AMP::UnitTest *utils )
{
    auto vectora = d_factory->getVector();
    auto vectorb = d_factory->getVector();
    auto vectorc = d_factory->getVector();

    vectora->setRandomValues();
    vectorb->copyVector( vectora );
    vectorc->subtract( vectora, vectorb );
    if ( vectorc->maxNorm() == 0 )
        utils->passes( "copy vector 1" );
    else
        utils->failure( "copy vector 1" );

    vectora->scale( 100. );
    vectorc->subtract( vectora, vectorb );
    auto c_maxNorm = vectorc->maxNorm();
    auto b_maxNorm = vectorb->maxNorm();
    if ( fabs( c_maxNorm - 99 * b_maxNorm ) < 1e-12 * b_maxNorm )
        utils->passes( "copy vector 2" );
    else
        utils->failure( "copy vector 2" );
}


void VectorTests::CopyRawDataBlockVector( AMP::UnitTest *utils )
{
    auto vectora = d_factory->getVector();
    auto vectorb = d_factory->getVector();
    auto vectorc = d_factory->getVector();
    vectora->setRandomValues();
    vectorb->zero();
    auto buf = new double[vectora->getLocalSize()];
    vectora->copyOutRawData( buf );
    vectorb->putRawData( buf );
    delete[] buf;
    vectorc->subtract( vectora, vectorb );
    if ( vectorc->maxNorm() == 0 )
        utils->passes( "copy raw data block" );
    else
        utils->failure( "copy raw data block" );
}


void VectorTests::VerifyVectorGhostCreate( AMP::UnitTest *utils )
{
    AMP_MPI globalComm( AMP_COMM_WORLD );
    auto vector    = d_factory->getVector();
    int num_ghosts = vector->getGhostSize();
    num_ghosts     = globalComm.sumReduce( num_ghosts );
    if ( globalComm.getSize() == 1 )
        utils->expected_failure( "No ghost cells for single processor " + d_factory->name() );
    else if ( num_ghosts > 0 )
        utils->passes( "verify ghosts created " + d_factory->name() );
    else
        utils->failure( "verify ghosts created " + d_factory->name() );
}


void VectorTests::VerifyVectorMakeConsistentAdd( AMP::UnitTest *utils )
{
    AMP_MPI globalComm( AMP_COMM_WORLD );
    auto dofmap = d_factory->getDOFMap();
    auto vector = d_factory->getVector();
    if ( !vector )
        utils->failure( "verify makeConsistent () for add " + d_factory->name() );

    // Zero the vector
    vector->zero();
    if ( vector->getUpdateStatus() != AMP::LinearAlgebra::Vector::UpdateState::UNCHANGED )
        utils->failure( "zero leaves vector in UpdateState::UNCHANGED state " + d_factory->name() );

    // Set and add local values by global id (this should not interfer with the add)
    for ( size_t i = dofmap->beginDOF(); i != dofmap->endDOF(); i++ ) {
        vector->setLocalValueByGlobalID( i, 0.0 );
        vector->addLocalValueByGlobalID( i, 0.0 );
    }
    if ( vector->getUpdateStatus() != AMP::LinearAlgebra::Vector::UpdateState::LOCAL_CHANGED )
        utils->failure( "local set/add leaves vector in UpdateState::LOCAL_CHANGED state " +
                        d_factory->name() );

    // Add values by global id
    for ( size_t i = dofmap->beginDOF(); i != dofmap->endDOF(); i++ )
        vector->addValueByGlobalID( i, (double) i );
    if ( vector->getUpdateStatus() != AMP::LinearAlgebra::Vector::UpdateState::ADDING )
        utils->failure( "addValueByGlobalID leaves vector in UpdateState::ADDING state " +
                        d_factory->name() );

    auto offset = (double) ( 1 << globalComm.getRank() );
    for ( size_t i = 0; i != vector->getGhostSize(); i++ ) {
        size_t ndx = vector->getCommunicationList()->getGhostIDList()[i];
        vector->addValueByGlobalID( ndx, offset );
    }

    // Perform a makeConsistent ADD and check the result
    vector->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_ADD );
    if ( vector->getUpdateStatus() != AMP::LinearAlgebra::Vector::UpdateState::UNCHANGED )
        utils->failure( "makeConsistent leaves vector in UpdateState::UNCHANGED state " +
                        d_factory->name() );
    std::map<int, std::set<size_t>> ghosted_entities;
    for ( size_t i = dofmap->beginDOF(); i != dofmap->endDOF(); i++ ) {
        double diff_double = fabs( vector->getValueByGlobalID( i ) - (double) i );
        if ( diff_double > 0.00001 ) {
            int ioffset  = lround( diff_double );
            int cur_rank = 0;
            while ( ioffset > 0 ) {
                if ( ioffset & 1 )
                    ghosted_entities[cur_rank].insert( i );
                ioffset >>= 1;
                cur_rank++;
            }
        }
    }
    auto cur_replicated = vector->getCommunicationList()->getReplicatedIDList().begin();
    auto end_replicated = vector->getCommunicationList()->getReplicatedIDList().end();
    while ( cur_replicated != end_replicated ) {
        bool found = false;
        for ( int i = 0; i != globalComm.getSize(); i++ ) {
            auto location = ghosted_entities[i].find( *cur_replicated );
            if ( location != ghosted_entities[i].end() ) {
                found = true;
                ghosted_entities[i].erase( location );
                break;
            }
        }
        if ( !found ) {
            utils->failure( "overly ghosted value " + d_factory->name() );
            return;
        }
        ++cur_replicated;
    }
    size_t last_size = 0;
    for ( int i = 0; i != globalComm.getSize(); i++ )
        last_size += ghosted_entities[i].size();
    if ( last_size == 0 )
        utils->passes( "all ghosted values accounted for " + d_factory->name() );
    else
        utils->failure( "some ghosted values not set " + d_factory->name() );
}


void VectorTests::VerifyVectorMakeConsistentSet( AMP::UnitTest *utils )
{
    auto dofmap = d_factory->getDOFMap();
    auto vector = d_factory->getVector();

    // Zero the vector
    vector->zero();
    if ( vector->getUpdateStatus() != AMP::LinearAlgebra::Vector::UpdateState::UNCHANGED )
        utils->failure( "zero leaves vector in UpdateState::UNCHANGED state " + d_factory->name() );

    // Set and add local values by global id (this should not interfere with the add)
    for ( size_t i = dofmap->beginDOF(); i != dofmap->endDOF(); i++ ) {
        vector->setLocalValueByGlobalID( i, 0.0 );
        vector->addLocalValueByGlobalID( i, 0.0 );
    }
    if ( vector->getUpdateStatus() != AMP::LinearAlgebra::Vector::UpdateState::LOCAL_CHANGED )
        utils->failure( "local set/add leaves vector in UpdateState::LOCAL_CHANGED state " +
                        d_factory->name() );

    // Set values by global id
    for ( size_t i = dofmap->beginDOF(); i != dofmap->endDOF(); i++ )
        vector->setValueByGlobalID( i, (double) i );
    if ( vector->getUpdateStatus() != AMP::LinearAlgebra::Vector::UpdateState::LOCAL_CHANGED &&
         vector->getUpdateStatus() != AMP::LinearAlgebra::Vector::UpdateState::SETTING )
        utils->failure( "setValueByGlobalID leaves vector in UpdateState::SETTING or "
                        "UpdateState::LOCAL_CHANGED state " +
                        d_factory->name() );

    // Perform a makeConsistent SET and check the result
    vector->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
    if ( vector->getUpdateStatus() != AMP::LinearAlgebra::Vector::UpdateState::UNCHANGED )
        utils->failure( "makeConsistent leaves vector in UpdateState::UNCHANGED state " +
                        d_factory->name() );
    if ( vector->getGhostSize() > 0 ) {
        auto comm_list = vector->getCommunicationList();
        std::vector<double> ghostList( vector->getGhostSize() );
        auto ghostIDList = comm_list->getGhostIDList();
        vector->getValuesByGlobalID(
            vector->getGhostSize(), (size_t *) &( ghostIDList[0] ), &( ghostList[0] ) );
        bool testPassed = true;
        for ( size_t i = 0; i != vector->getGhostSize(); i++ ) {
            if ( fabs( ghostList[i] - (double) ( ghostIDList[i] ) ) > 0.0000001 )
                testPassed = false;
        }
        if ( testPassed )
            utils->passes( "ghost set correctly in vector " + d_factory->name() );
        else
            utils->failure( "ghost not set correctly in vector " + d_factory->name() );
    }
    if ( vector->getGhostSize() > 0 ) {
        auto comm_list   = vector->getCommunicationList();
        auto ghostIDList = comm_list->getGhostIDList();
        bool testPassed  = true;
        for ( size_t i = 0; i != vector->getGhostSize(); i++ ) {
            size_t ghostNdx = ghostIDList[i];
            double ghostVal = vector->getValueByGlobalID( ghostNdx );
            if ( fabs( ghostVal - (double) ghostNdx ) > 0.0000001 )
                testPassed = false;
        }
        if ( testPassed )
            utils->passes( "ghost set correctly in alias " + d_factory->name() );
        else
            utils->failure( "ghost set correctly in alias " + d_factory->name() );
    }
}


// Test creating a multivector with multiple copies of the data
// This should always return one copy of the superset of the data
void VectorTests::TestMultivectorDuplicate( AMP::UnitTest *utils )
{
    auto vec0 = d_factory->getVector();
    // Create a multivector
    auto var      = std::make_shared<AMP::LinearAlgebra::Variable>( "multivec" );
    auto multiVec = AMP::LinearAlgebra::MultiVector::create( var, vec0->getComm() );
    // Add different views of vec0
    multiVec->addVector( vec0 );
    multiVec->addVector( vec0 );
    multiVec->addVector( multiVec->getVector( 0 ) );
    auto var2 = std::make_shared<AMP::LinearAlgebra::Variable>( "vec2" );
    multiVec->addVector( AMP::LinearAlgebra::MultiVector::create( var2, vec0->getComm() ) );
    // Verify the size of the multivector
    auto dof1 = vec0->getDOFManager();
    auto dof2 = multiVec->getDOFManager();
    bool pass = dof1->numLocalDOF() == dof2->numLocalDOF() &&
                dof1->numGlobalDOF() == dof2->numGlobalDOF() &&
                dof1->beginDOF() == dof2->beginDOF();
    if ( pass )
        utils->passes( "multivector resolves multiple copies of a vector " + d_factory->name() );
    else
        utils->failure( "multivector resolves multiple copies of a vector " + d_factory->name() );
}


} // namespace LinearAlgebra
} // namespace AMP
