#include "AMP/vectors/testHelpers/VectorTests.h"
#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/Vector.h"
#ifdef USE_EXT_SUNDIALS
    #include "AMP/vectors/sundials/ManagedSundialsVector.h"
    #include "AMP/vectors/sundials/SundialsVector.h"
#endif
#ifdef USE_EXT_PETSC
    #include "AMP/vectors/petsc/PetscVector.h"
#endif

#include <algorithm>


#define PASS_FAIL( test, MSG )                                                    \
    do {                                                                          \
        if ( test )                                                               \
            ut->passes( d_factory->name() + " - " + __FUNCTION__ + ": " + MSG );  \
        else                                                                      \
            ut->failure( d_factory->name() + " - " + __FUNCTION__ + ": " + MSG ); \
    } while ( 0 )


namespace AMP {
namespace LinearAlgebra {


static inline int lround( double x ) { return x >= 0 ? floor( x ) : ceil( x ); }

void VectorTests::InstantiateVector( AMP::UnitTest *ut )
{
    auto vector = d_factory->getVector();
    PASS_FAIL( vector, "created " );
}


void VectorTests::CopyVectorConsistency( AMP::UnitTest *ut )
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
        PASS_FAIL( std::equal( t1, t1 + numGhosts, t2 ), "Ghosts are the same (1)" );
    }

    vec1->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
    vec3->copyVector( vec1 );
    if ( numGhosts ) {
        vec1->getValuesByGlobalID( numGhosts, ndx, t1 );
        vec3->getValuesByGlobalID( numGhosts, ndx, t2 );
        PASS_FAIL( std::equal( t1, t1 + numGhosts, t2 ), "Ghosts are the same (2)" );
        delete[] t1;
        delete[] t2;
        delete[] ndx;
    }
}


void VectorTests::Bug_728( AMP::UnitTest *ut )
{
    auto vector = d_factory->getVector();
    auto var1   = vector->getVariable();
    if ( !var1 )
        return;
    auto var2 = var1->cloneVariable( var1->getName() );
    PASS_FAIL( vector->subsetVectorForVariable( var1 ), "Found vector for same variable" );
    PASS_FAIL( vector->subsetVectorForVariable( var2 ), "Found vector for cloned variable" );
}


void VectorTests::SetToScalarVector( AMP::UnitTest *ut )
{
    auto vector = d_factory->getVector();
    vector->setToScalar( 0. );
    vector->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
    ut->passes( "setToScalar ran to completion " + d_factory->name() );
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
    PASS_FAIL( !fail, "Set data to 0" );
    fail = false;
    vector->setToScalar( 5. );
    vector->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
    auto curVal = vector->begin();
    while ( curVal != endVec ) {
        if ( *curVal != 5. ) {
            fail = true;
            break;
        }
        ++curVal;
    }
    PASS_FAIL( !fail, "Set data to 5" );
    auto remoteDofs = vector->getDOFManager()->getRemoteDOFs();
    fail            = false;
    for ( auto &remoteDof : remoteDofs ) {
        if ( vector->getValueByGlobalID( remoteDof ) != 5. )
            fail = true;
    }
    PASS_FAIL( !fail, "Set ghost data to 5" );
}


void VectorTests::CloneVector( AMP::UnitTest *ut )
{
    auto vector = d_factory->getVector();
    auto clone  = vector->cloneVector( "cloned vector" );
    vector->setToScalar( 3.0 );
    clone->setToScalar( 0.0 );
    ut->passes( "Clone created " + d_factory->name() );
    bool pass = true;
    for ( size_t i = 0; i != vector->numberOfDataBlocks(); i++ ) {
        auto *clone_ptr  = clone->getRawDataBlock<double>( i );
        auto *vector_ptr = vector->getRawDataBlock<double>( i );
        if ( clone_ptr == vector_ptr )
            pass = false;
    }
    PASS_FAIL( pass, "CloneVector: allocated" );
    clone->setToScalar( 1. );
    double t = static_cast<double>( clone->L1Norm() );
    PASS_FAIL( clone->getGlobalSize() == vector->getGlobalSize(),
               "CloneVector: global size equality" );
    PASS_FAIL( clone->getLocalSize() == vector->getLocalSize(),
               "CloneVector: local size equality" );
    PASS_FAIL( fabs( t - (double) clone->getGlobalSize() ) < 0.0000001,
               "CloneVector: trivial set data" );
}


void VectorTests::DotProductVector( AMP::UnitTest *ut )
{
    auto vector1( d_factory->getVector() );
    auto vector2( d_factory->getVector() );
    vector1->setToScalar( 1. );
    vector2->setToScalar( 2. );
    auto d11 = vector1->dot( *vector1 );
    auto d12 = vector1->dot( *vector2 );
    auto d21 = vector2->dot( *vector1 );
    auto d22 = vector2->dot( *vector2 );
    PASS_FAIL( 2 * d11 == d12, "dot product 1" );
    PASS_FAIL( 2 * d11 == d21, "dot product 2" );
    PASS_FAIL( 4 * d11 == d22, "dot product 3" );
    PASS_FAIL( d11 == vector1->getGlobalSize(), "dot product 4" );
    PASS_FAIL( d21 == d12, "dot product 5" );
}


void VectorTests::L2NormVector( AMP::UnitTest *ut )
{
    auto vector = d_factory->getVector();
    vector->setToScalar( 1. );
    auto norm  = static_cast<double>( vector->L2Norm() );
    auto norm2 = static_cast<double>( vector->dot( *vector ) );
    PASS_FAIL( fabs( norm * norm - norm2 ) < 0.000001, "L2 norm 1" );
    vector->setRandomValues();
    norm  = static_cast<double>( vector->L2Norm() );
    norm2 = static_cast<double>( vector->dot( *vector ) );
    PASS_FAIL( fabs( norm * norm - norm2 ) < 0.000001, "L2 norm 2" );
}


void VectorTests::AbsVector( AMP::UnitTest *ut )
{
    auto vec1 = d_factory->getVector();
    auto vec2 = vec1->cloneVector();
    vec1->setRandomValues();
    vec2->copyVector( vec1 );
    vec2->scale( -1.0 );
    vec2->abs( *vec2 );
    PASS_FAIL( vec1->equals( *vec2 ), "Abs passes" );
}


void VectorTests::L1NormVector( AMP::UnitTest *ut )
{
    auto vector = d_factory->getVector();
    auto vector_1( d_factory->getVector() );
    vector->setRandomValues();
    vector_1->setToScalar( 1. );
    auto norm = static_cast<double>( vector->L1Norm() );
    vector->abs( *vector );
    auto norm2 = static_cast<double>( vector->dot( *vector_1 ) );
    PASS_FAIL( fabs( norm - norm2 ) < 0.000001, "L1 norm" );
}


void VectorTests::MaxNormVector( AMP::UnitTest *ut )
{
    auto vector = d_factory->getVector();
    vector->setRandomValues();
    auto infNorm = vector->maxNorm();
    vector->abs( *vector );
    auto curData     = vector->begin();
    auto endData     = vector->end();
    double local_ans = *curData;
    while ( curData != endData ) {
        local_ans = std::max( local_ans, *curData );
        ++curData;
    }
    double global_ans = vector->getComm().maxReduce( local_ans );
    PASS_FAIL( global_ans == infNorm, "Inf norm" );
}


void VectorTests::ScaleVector( AMP::UnitTest *ut )
{
    auto vector1( d_factory->getVector() );
    auto vector2( d_factory->getVector() );
    double beta = 1.2345;
    vector2->setRandomValues();
    vector1->scale( beta, *vector2 );
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
    PASS_FAIL( pass, "scale vector 1" );
    vector2->scale( beta );
    vector1->subtract( *vector2, *vector1 );
    PASS_FAIL( vector1->maxNorm() < 0.0000001, "scale vector 2" );
}


#ifdef USE_EXT_PETSC
void VectorTests::Bug_491( AMP::UnitTest *ut )
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
    double sp_n1  = static_cast<double>( vector1->L1Norm().get<double>() );
    double sp_n2  = static_cast<double>( vector1->L2Norm().get<double>() );
    double sp_inf = static_cast<double>( vector1->maxNorm().get<double>() );

    // Check to see if petsc cache has been invalidated
    VecNormBegin( managed_vec, NORM_1, &n1 );
    VecNormBegin( managed_vec, NORM_2, &n2 );
    VecNormBegin( managed_vec, NORM_INFINITY, &ninf );
    VecNormEnd( managed_vec, NORM_1, &n1 );
    VecNormEnd( managed_vec, NORM_2, &n2 );
    VecNormEnd( managed_vec, NORM_INFINITY, &ninf );

    PASS_FAIL( fabs( n1 - sp_n1 ) < 0.00000001 * n1, "L1 norm -- Petsc interface begin/end" );
    PASS_FAIL( fabs( n2 - sp_n2 ) < 0.00000001 * n1, "L2 norm -- Petsc interface begin/end" );
    PASS_FAIL( fabs( ninf - sp_inf ) < 0.00000001 * n1, "Linf norm -- Petsc interface begin/end" );

    VecNorm( managed_vec, NORM_1, &n1 );
    VecNorm( managed_vec, NORM_2, &n2 );
    VecNorm( managed_vec, NORM_INFINITY, &ninf );

    double L1Norm( vector1->L1Norm() );
    double L2Norm( vector1->L2Norm() );
    double maxNorm( vector1->maxNorm() );
    PASS_FAIL( fabs( n1 - L1Norm ) < 0.00000001 * n1, "L1 norm -- Petsc interface begin/end " );
    PASS_FAIL( fabs( n2 - L2Norm ) < 0.00000001 * n1, "L2 norm -- Petsc interface begin/end " );
    PASS_FAIL( fabs( ninf - maxNorm ) < 0.00000001 * n1, "inf norm -- Petsc interface begin/end " );
}
#endif


void VectorTests::AddVector( AMP::UnitTest *ut )
{
    auto vector1( d_factory->getVector() );
    auto vector2( d_factory->getVector() );
    auto vector3( d_factory->getVector() );
    vector1->setRandomValues();
    vector2->setRandomValues();
    vector3->add( *vector1, *vector2 );
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
    PASS_FAIL( pass, "add vector" );
}


void VectorTests::SubtractVector( AMP::UnitTest *ut )
{
    auto vector1( d_factory->getVector() );
    auto vector2( d_factory->getVector() );
    auto vector3( d_factory->getVector() );
    auto vector4( d_factory->getVector() );
    vector1->setRandomValues();
    vector2->setRandomValues();
    vector3->subtract( *vector1, *vector2 );
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
    PASS_FAIL( pass, "vector subtract 1" );
    vector2->scale( -1. );
    vector4->add( *vector1, *vector2 );
    vector4->subtract( *vector3, *vector4 );
    PASS_FAIL( vector4->maxNorm() < 0.0000001, "vector subtract 2" );
}


void VectorTests::MultiplyVector( AMP::UnitTest *ut )
{
    auto vector1( d_factory->getVector() );
    auto vector2( d_factory->getVector() );
    auto vector3( d_factory->getVector() );
    vector1->setRandomValues();
    vector2->setToScalar( 3. );
    vector3->multiply( *vector1, *vector2 );
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
    PASS_FAIL( pass, "vector::multiply" );
}


void VectorTests::DivideVector( AMP::UnitTest *ut )
{
    auto vector1( d_factory->getVector() );
    auto vector2( d_factory->getVector() );
    auto vector3( d_factory->getVector() );
    vector1->setRandomValues();
    vector2->setRandomValues();
    vector3->divide( *vector1, *vector2 );
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
    PASS_FAIL( pass, "vector::divide" );
}


void VectorTests::VectorIteratorLengthTest( AMP::UnitTest *ut )
{
    auto vector1( d_factory->getVector() );
    size_t i = 0;
    for ( auto it = vector1->begin(); it != vector1->end(); ++it )
        i++;
    size_t k = vector1->getLocalSize();
    PASS_FAIL( i == k, "Iterated over the correct number of entries" );
}


void VectorTests::VectorIteratorTests( AMP::UnitTest *ut )
{
    auto vector1 = d_factory->getVector();
    both_VectorIteratorTests<AMP::LinearAlgebra::VectorDataIterator<double>>( vector1, ut );
    both_VectorIteratorTests<AMP::LinearAlgebra::VectorDataIterator<const double>>( vector1, ut );
}


void VectorTests::VerifyVectorMin( AMP::UnitTest *ut )
{
    auto vec = d_factory->getVector();
    vec->setRandomValues();
    vec->scale( -1.0 ); // make negative
    PASS_FAIL( ( vec->min() + vec->maxNorm() ).abs() < 1.e-10,
               "minimum of negative vector == ||.||_infty" );
}


void VectorTests::VerifyVectorMax( AMP::UnitTest *ut )
{
    auto vec = d_factory->getVector();
    vec->setRandomValues();
    PASS_FAIL( ( vec->max() - vec->maxNorm().abs() ) < 1.e-10,
               "maximum of positive vector == ||.||_infty" );
}


void VectorTests::VerifyVectorMaxMin( AMP::UnitTest *ut )
{
    auto vec    = d_factory->getVector();
    bool passes = true;
    for ( size_t i = 0; i != 10; i++ ) {
        vec->setRandomValues();
        vec->addScalar( *vec, -0.5 );
        vec->scale( 2.0 ); // vec i.i.d [-1,1);
        double max( vec->max() );
        double min( vec->min() );
        auto ans = std::max( fabs( max ), fabs( min ) );
        double norm( vec->maxNorm() );
        if ( fabs( ans - norm ) >= 1.e-20 )
            passes = false;
    }
    PASS_FAIL( passes, "Max and min correctly predict maxNorm()" );
}


void VectorTests::SetRandomValuesVector( AMP::UnitTest *ut )
{
    auto vector  = d_factory->getVector();
    auto l2norm1 = -1;
    for ( size_t i = 0; i < 5; i++ ) {
        vector->setRandomValues();
        auto l2norm2 = static_cast<double>( vector->L2Norm() );
        PASS_FAIL( fabs( l2norm1 - l2norm2 ) > 0.000001, "Distinct vector created" );
        l2norm1 = l2norm2;
        PASS_FAIL( vector->min() >= 0, "SetRandomValuesVector: Min value >= 0" );
        PASS_FAIL( vector->max() < 1, "SetRandomValuesVector: Max value < 1" );
        PASS_FAIL( vector->L2Norm() > 0, "SetRandomValuesVector: Non-zero vector created" );
    }
}


void VectorTests::ReciprocalVector( AMP::UnitTest *ut )
{
    auto vectora = d_factory->getVector();
    auto vectorb = d_factory->getVector();
    auto vectorc = d_factory->getVector();
    auto vectord = d_factory->getVector();
    auto vector1 = d_factory->getVector();
    vectora->setRandomValues();
    vectorb->reciprocal( *vectora );
    vector1->setToScalar( 1. );
    vectorc->divide( *vector1, *vectora );
    vectord->subtract( *vectorb, *vectorc );
    PASS_FAIL( vectord->maxNorm() < 0.0000001, "vector::reciprocal" );
}


static void LinearSumVectorRun( std::shared_ptr<const VectorFactory> d_factory,
                                AMP::UnitTest *ut,
                                double alpha,
                                double beta,
                                const char *msg )
{
    auto vectora = d_factory->getVector();
    auto vectorb = d_factory->getVector();
    auto vectorc = d_factory->getVector();
    auto vectord = d_factory->getVector();
    vectora->setRandomValues();
    vectorb->setRandomValues();
    vectorc->linearSum( alpha, *vectora, beta, *vectorb );
    vectora->scale( alpha );
    vectorb->scale( beta );
    vectord->add( *vectora, *vectorb );
    vectord->subtract( *vectorc, *vectord );
    PASS_FAIL( vectord->maxNorm() < 0.0000001, msg );
}
void VectorTests::LinearSumVector( AMP::UnitTest *ut )
{
    LinearSumVectorRun( d_factory, ut, 1.2345, 9.8765, "linear sum 1" );
    LinearSumVectorRun( d_factory, ut, -1.2345, 9.8765, "linear sum 2" );
    LinearSumVectorRun( d_factory, ut, 1.2345, -9.8765, "linear sum 3" );
    LinearSumVectorRun( d_factory, ut, -1.2345, -9.8765, "linear sum 4" );
}


static void AxpyVectorRun( std::shared_ptr<const VectorFactory> d_factory,
                           AMP::UnitTest *ut,
                           double alpha,
                           const char *msg )
{
    auto vectora = d_factory->getVector();
    auto vectorb = d_factory->getVector();
    auto vectorc = d_factory->getVector();
    auto vectord = d_factory->getVector();
    vectora->setRandomValues();
    vectorb->setRandomValues();
    vectorc->linearSum( alpha, *vectora, 1., *vectorb );
    vectord->axpy( alpha, *vectora, *vectorb );
    vectorc->subtract( *vectorc, *vectord );
    PASS_FAIL( vectorc->maxNorm() < 0.0000001, msg );
}
void VectorTests::AxpyVector( AMP::UnitTest *ut )
{
    AxpyVectorRun( d_factory, ut, 6.38295, "axpy 1" );
    AxpyVectorRun( d_factory, ut, -6.38295, "axpy 2" );
}


static void AxpbyVectorRun( std::shared_ptr<const VectorFactory> d_factory,
                            AMP::UnitTest *ut,
                            double alpha,
                            double beta,
                            const char *msg )
{
    auto vectora  = d_factory->getVector();
    auto vectorb  = d_factory->getVector();
    auto vectorb1 = d_factory->getVector();
    auto vectorc  = d_factory->getVector();
    auto vectord  = d_factory->getVector();
    vectora->setRandomValues();
    vectorb->setRandomValues();
    vectorc->copyVector( vectorb );
    vectorb1->linearSum( alpha, *vectora, beta, *vectorb );
    vectorb->linearSum( alpha, *vectora, beta, *vectorb );
    vectorc->axpby( alpha, beta, *vectora );
    vectord->subtract( *vectorc, *vectorb1 );
    auto maxNorm = vectord->maxNorm();
    PASS_FAIL( maxNorm < 0.0000001, msg );
    vectord->subtract( *vectorc, *vectorb );
    maxNorm = vectord->maxNorm();
    PASS_FAIL( maxNorm < 0.0000001, msg );
}
void VectorTests::AxpbyVector( AMP::UnitTest *ut )
{
    AxpbyVectorRun( d_factory, ut, 6.38295, 99.273624, "axpby 1" );
    AxpbyVectorRun( d_factory, ut, 6.38295, -99.273624, "axpby 2" );
    AxpbyVectorRun( d_factory, ut, -6.38295, 99.273624, "axpby 3" );
    AxpbyVectorRun( d_factory, ut, -6.38295, -99.273624, "axpby 4" );
}


void VectorTests::CopyVector( AMP::UnitTest *ut )
{
    auto vectora = d_factory->getVector();
    auto vectorb = d_factory->getVector();
    auto vectorc = d_factory->getVector();

    vectora->setRandomValues();
    vectorb->copyVector( vectora );
    vectorc->subtract( *vectora, *vectorb );
    PASS_FAIL( vectorc->maxNorm() == 0, "copy vector 1" );

    vectora->scale( 100. );
    vectorc->subtract( *vectora, *vectorb );
    auto c_maxNorm = vectorc->maxNorm().get<double>();
    auto b_maxNorm = vectorb->maxNorm().get<double>();
    PASS_FAIL( fabs( c_maxNorm - 99 * b_maxNorm ) < 1e-12 * b_maxNorm, "copy vector 2" );
}


void VectorTests::CopyRawDataBlockVector( AMP::UnitTest *ut )
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
    vectorc->subtract( *vectora, *vectorb );
    PASS_FAIL( vectorc->maxNorm() == 0, "copy raw data block" );
}


void VectorTests::VerifyVectorGhostCreate( AMP::UnitTest *ut )
{
    AMP_MPI globalComm( AMP_COMM_WORLD );
    auto vector    = d_factory->getVector();
    int num_ghosts = vector->getGhostSize();
    num_ghosts     = globalComm.sumReduce( num_ghosts );
    if ( globalComm.getSize() == 1 ) {
        ut->expected_failure( "No ghost cells for single processor " + d_factory->name() );
    } else {
        PASS_FAIL( num_ghosts > 0, "verify ghosts created " );
    }
}


void VectorTests::VerifyVectorMakeConsistentAdd( AMP::UnitTest *ut )
{
    AMP_MPI globalComm( AMP_COMM_WORLD );
    auto vector = d_factory->getVector();
    auto dofmap = vector->getDOFManager();
    if ( !vector )
        ut->failure( "verify makeConsistent () for add " + d_factory->name() );

    // Zero the vector
    vector->zero();
    if ( vector->getUpdateStatus() != AMP::LinearAlgebra::VectorData::UpdateState::UNCHANGED )
        ut->failure( "zero leaves vector in UpdateState::UNCHANGED state " + d_factory->name() );

    // Set and add local values by global id (this should not interfer with the add)
    const double val = 0.0;
    for ( size_t i = dofmap->beginDOF(); i != dofmap->endDOF(); i++ ) {
        vector->setLocalValuesByGlobalID( 1, &i, &val );
        vector->addLocalValuesByGlobalID( 1, &i, &val );
    }
    if ( vector->getUpdateStatus() != AMP::LinearAlgebra::VectorData::UpdateState::LOCAL_CHANGED )
        ut->failure( "local set/add leaves vector in UpdateState::LOCAL_CHANGED state " +
                     d_factory->name() );

    // Add values by global id
    for ( size_t i = dofmap->beginDOF(); i != dofmap->endDOF(); i++ ) {
        const double val = double( i );
        vector->addValuesByGlobalID( 1, &i, &val );
    }
    if ( vector->getUpdateStatus() != AMP::LinearAlgebra::VectorData::UpdateState::ADDING )
        ut->failure( "addValueByGlobalID leaves vector in UpdateState::ADDING state " +
                     d_factory->name() );

    auto offset = (double) ( 1 << globalComm.getRank() );
    for ( size_t i = 0; i != vector->getGhostSize(); i++ ) {
        size_t ndx = vector->getCommunicationList()->getGhostIDList()[i];
        vector->addValuesByGlobalID( 1, &ndx, &offset );
    }

    // Perform a makeConsistent ADD and check the result
    vector->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_ADD );
    if ( vector->getUpdateStatus() != AMP::LinearAlgebra::VectorData::UpdateState::UNCHANGED )
        ut->failure( "makeConsistent leaves vector in UpdateState::UNCHANGED state " +
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
            ut->failure( "overly ghosted value " + d_factory->name() );
            return;
        }
        ++cur_replicated;
    }
    size_t last_size = 0;
    for ( int i = 0; i != globalComm.getSize(); i++ )
        last_size += ghosted_entities[i].size();
    PASS_FAIL( last_size == 0, "all ghosted values accounted for" );
}


void VectorTests::VerifyVectorMakeConsistentSet( AMP::UnitTest *ut )
{
    auto vector = d_factory->getVector();
    auto dofmap = vector->getDOFManager();

    // Zero the vector
    vector->zero();
    if ( vector->getUpdateStatus() != AMP::LinearAlgebra::VectorData::UpdateState::UNCHANGED )
        ut->failure( "zero leaves vector in UpdateState::UNCHANGED state " + d_factory->name() );

    // Set and add local values by global id (this should not interfere with the add)
    const double val = 0.0;
    for ( size_t i = dofmap->beginDOF(); i != dofmap->endDOF(); i++ ) {
        vector->setLocalValuesByGlobalID( 1, &i, &val );
        vector->addLocalValuesByGlobalID( 1, &i, &val );
    }
    if ( vector->getUpdateStatus() != AMP::LinearAlgebra::VectorData::UpdateState::LOCAL_CHANGED )
        ut->failure( "local set/add leaves vector in UpdateState::LOCAL_CHANGED state " +
                     d_factory->name() + " - " + vector->type() );

    // Set values by global id
    for ( size_t i = dofmap->beginDOF(); i != dofmap->endDOF(); i++ ) {
        const double val = double( i );
        vector->setValuesByGlobalID( 1, &i, &val );
    }
    if ( vector->getUpdateStatus() != AMP::LinearAlgebra::VectorData::UpdateState::LOCAL_CHANGED &&
         vector->getUpdateStatus() != AMP::LinearAlgebra::VectorData::UpdateState::SETTING )
        ut->failure( "setValueByGlobalID leaves vector in UpdateState::SETTING or "
                     "UpdateState::LOCAL_CHANGED state " +
                     d_factory->name() );

    // Perform a makeConsistent SET and check the result
    vector->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
    if ( vector->getUpdateStatus() != AMP::LinearAlgebra::VectorData::UpdateState::UNCHANGED )
        ut->failure( "makeConsistent leaves vector in UpdateState::UNCHANGED state " +
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
        PASS_FAIL( testPassed, "ghost set correctly in vector" );
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
        PASS_FAIL( testPassed, "ghost set correctly in alias " );
    }
}


// Test creating a multivector with multiple copies of the data
// This should always return one copy of the superset of the data
void VectorTests::TestMultivectorDuplicate( AMP::UnitTest *ut )
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
    PASS_FAIL( pass, "multivector resolves multiple copies of a vector" );
}


} // namespace LinearAlgebra
} // namespace AMP
