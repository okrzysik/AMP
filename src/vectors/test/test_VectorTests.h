#ifndef included_test_VectorTests
#define included_test_VectorTests

#include <algorithm>
#include "ampmesh/Mesh.h"
#include "utils/AMP_MPI.h"
#include "utils/UnitTest.h"
#include "discretization/DOF_Manager.h"
#include "vectors/Vector.h"
#include "vectors/MultiVector.h"
#ifdef USE_EXT_SUNDIALS
    #include "vectors/sundials/ManagedSundialsVector.h"
    #include "vectors/sundials/SundialsVector.h"
#endif
#ifdef USE_EXT_PETSC
    #include "vectors/petsc/ManagedPetscVector.h"
    #include "vectors/petsc/PetscVector.h"
#endif


/// \cond UNDOCUMENTED
namespace AMP {
namespace unit_test {


static inline int lround(double x) { return x>=0 ? floor(x):ceil(x); }


template <typename VECTOR_FACTORY>
void InstantiateVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr  vector ( VECTOR_FACTORY::getVector() );
    if ( vector )
        utils->passes( "created" );
    else
        utils->failure( "created" );
}


template <typename VECTOR_FACTORY>
void CopyVectorConsistency( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr  vec1 ( VECTOR_FACTORY::getVector() );
    AMP::LinearAlgebra::Vector::shared_ptr  vec2 = vec1->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr  vec3 = vec1->cloneVector();
    AMP::LinearAlgebra::CommunicationList::shared_ptr  commList = vec1->getCommunicationList();
    double *t1=NULL;
    double *t2=NULL;
    size_t *ndx=NULL;
    size_t numGhosts = commList->getGhostIDList().size();

    vec1->setRandomValues ();
    vec2->copyVector ( vec1 );
    if ( numGhosts ) {
        t1 = new double[ numGhosts ];
        t2 = new double[ numGhosts ];
        ndx = new size_t[ numGhosts ];
        std::copy ( commList->getGhostIDList().begin() ,
                    commList->getGhostIDList().end() ,
                    ndx );
        vec1->getValuesByGlobalID ( numGhosts , ndx , t1 );
        vec2->getValuesByGlobalID ( numGhosts , ndx , t2 );
        if ( std::equal ( t1 , t1 + numGhosts , t2 ) )
            utils->passes ( "Ghosts are the same" );
        else
            utils->failure ( "Ghosts are different" );
    }

    vec1->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    vec3->copyVector ( vec1 );
    if ( numGhosts ) {
        vec1->getValuesByGlobalID ( numGhosts , ndx , t1 );
        vec3->getValuesByGlobalID ( numGhosts , ndx , t2 );
        if ( std::equal ( t1 , t1 + numGhosts , t2 ) )
            utils->passes ( "Ghosts are the same" );
        else
            utils->failure ( "Ghosts are different" );
        delete [] t1;
        delete [] t2;
        delete [] ndx;
    }

}


template <typename VECTOR_FACTORY, typename VIEWER>
void DeepCloneOfView( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr  vector1 ( VECTOR_FACTORY::getVector() );
    if ( !vector1->isA<AMP::LinearAlgebra::MultiVector>() ) return;
    vector1 = VIEWER::view ( vector1 );
    AMP::LinearAlgebra::Vector::shared_ptr  vector2 = vector1->cloneVector ();
    bool pass = true;
    for ( size_t i = 0 ; i != vector1->numberOfDataBlocks() ; i++ ) {
        pass &= (vector1->getRawDataBlock<double>(i) != vector2->getRawDataBlock<double>(i));
    }
    if ( pass )
        utils->passes ( "Deep clone succeeded" );
    else
        utils->failure ( "Deep clone failed" );
}


template <typename VECTOR_FACTORY>
void Bug_728( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr  vector ( VECTOR_FACTORY::getVector() );
    AMP::LinearAlgebra::Variable::shared_ptr var1 = vector->getVariable ();

    // Exit if there is no associated variable
    if ( !var1 )
        return;

    AMP::LinearAlgebra::Variable::shared_ptr var2 = var1->cloneVariable ( var1->getName() );

    if ( vector->subsetVectorForVariable ( var1 ) )
        utils->passes ( "Found vector for same variable pointer" );
    else
        utils->failure ( "Did not find vector for same variable pointer" );

    if ( vector->subsetVectorForVariable ( var2 ) )
        utils->passes ( "Found vector for cloned variable pointer" );
    else
        utils->failure ( "Did not find vector for cloned variable pointer" );
}


template <typename VECTOR_FACTORY>
void SetToScalarVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr  vector ( VECTOR_FACTORY::getVector() );
    vector->setToScalar ( 0. );
    utils->passes( "setToScalar ran to completion" );
    bool fail = false;
    AMP::LinearAlgebra::Vector::iterator  curVec = vector->begin();
    AMP::LinearAlgebra::Vector::iterator  endVec = vector->end();
    while ( curVec != endVec ) {
        if ( *curVec != 0. ) {
            fail = true;
            break;
        }
        ++curVec;
    }
    if ( !fail )
        utils->passes ( "Set data to 0" );
    else
        utils->failure ( "Failed to set scalar to 0" );
    fail = false;
    vector->setToScalar ( 5. );
    AMP::LinearAlgebra::Vector::iterator curVal = vector->begin();
    while ( curVal != endVec ) {
        if ( *curVal != 5. ) {
            fail = true;
            break;
        }
        ++curVal;
    }
    if ( !fail )
      utils->passes ( "Set data to 5" );
    else 
      utils->failure ( "Failed to set scalar to 5" );
    std::vector<size_t> remoteDofs = vector->getDOFManager()->getRemoteDOFs();
    fail = false;
    for (size_t i=0; i<remoteDofs.size(); i++) {
        if ( vector->getValueByGlobalID(remoteDofs[i])!=5. )
            fail = true;
    }
    if ( !fail )
      utils->passes ( "Set ghost data to 5" );
    else 
      utils->failure ( "Failed to set ghost scalar values to 5" );
}


template <typename VECTOR_FACTORY>
void CloneVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr  vector ( VECTOR_FACTORY::getVector() );
    AMP::LinearAlgebra::Vector::shared_ptr  clone = vector->cloneVector ( "cloned vector" );
    clone->setToScalar ( 0. );
    utils->passes( "Clone created" );
    bool  pass = true;
    for ( size_t i = 0 ; i != vector->numberOfDataBlocks() ; i++ ) {
        double *clone_ptr = clone->getRawDataBlock<double>(i);
        double *vector_ptr = vector->getRawDataBlock<double>(i);
        if ( clone_ptr == vector_ptr )
            pass = false;
    }
    if ( pass )
        utils->passes ( "New data allocated" );
    else
        utils->failure ( "Data not allocated" );
}


template <typename VECTOR_FACTORY>
void DotProductVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr  vector1 ( VECTOR_FACTORY::getVector() );
    AMP::LinearAlgebra::Vector::shared_ptr  vector2 ( VECTOR_FACTORY::getVector() );
    vector1->setToScalar ( 1. );
    vector2->setToScalar ( 2. );
    double  d11 , d12 , d21 , d22;
    d11 = vector1->dot ( vector1 );
    d12 = vector1->dot ( vector2 );
    d21 = vector2->dot ( vector1 );
    d22 = vector2->dot ( vector2 );
    if ( 2.*d11 == d12 )
        utils->passes ( "dot product 1" );
    else
        utils->failure ( "dot product 1" );
    if ( 2.*d11 == d21 )
        utils->passes ( "dot product 2" );
    else
        utils->failure ( "dot product 2" );
    if ( 4.*d11 == d22 )
        utils->passes ( "dot product 3" );
    else
        utils->failure ( "dot product 3" );
    if ( d11 == vector1->getGlobalSize() )
        utils->passes ( "dot product 4" );
    else
        utils->failure ( "dot product 4" );
    if ( d21 == d12 )
        utils->passes ( "dot product 5" );
    else
        utils->failure ( "dot product 5" );
}


template <typename VECTOR_FACTORY>
void L2NormVector( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr  vector ( VECTOR_FACTORY::getVector() );
    vector->setToScalar ( 1. );
    double  norm , norm2;
    norm = vector->L2Norm();
    norm2 = vector->dot ( vector );
    if ( fabs ( norm * norm - norm2 ) < 0.000001 )
        utils->passes ( "L2 norm 1" );
    else
        utils->failure ( "L2 norm 1" );
    vector->setRandomValues ();
    norm = vector->L2Norm();
    norm2 = vector->dot ( vector );
    if ( fabs ( norm * norm - norm2 ) < 0.000001 )
        utils->passes ( "L2 norm 2" );
    else
        utils->failure ( "L2 norm 2" );
}


template <typename VECTOR_FACTORY>
void AbsVector( AMP::UnitTest *utils )
{
        AMP::LinearAlgebra::Vector::shared_ptr  vec1 ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vec2 = vec1->cloneVector ();
        vec1->setRandomValues ();
        vec2->copyVector ( vec1 );
        vec2->scale ( -1.0 );
        vec2->abs ( vec2 );
        if ( vec1->equals ( vec2 ) )
          utils->passes ( "Abs passes" );
        else
          utils->failure ( "Abs fails" );
}


template <typename VECTOR_FACTORY>
void L1NormVector( AMP::UnitTest *utils )
{
        AMP::LinearAlgebra::Vector::shared_ptr  vector ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vector_1 ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vector_abs;
        vector->setRandomValues ();
        vector_1->setToScalar ( 1. );
        double  norm , norm2;
        norm = vector->L1Norm();
        vector->abs(vector);
        norm2 = vector->dot ( vector_1 );
        if ( fabs ( norm - norm2 ) < 0.000001 )
          utils->passes ( "L1 norm" );
        else
          utils->failure ( "L1 norm" );
}


template <typename VECTOR_FACTORY>
void MaxNormVector( AMP::UnitTest *utils )
{
        AMP::LinearAlgebra::Vector::shared_ptr  vector ( VECTOR_FACTORY::getVector() );
        vector->setRandomValues ();
        double infNorm = vector->maxNorm();
        vector->abs(vector);
        AMP::LinearAlgebra::Vector::iterator curData = vector->begin();
        AMP::LinearAlgebra::Vector::iterator endData = vector->end();
        double local_ans = *curData;
        while ( curData != endData )
        {
          local_ans = std::max ( local_ans , *curData );
          ++curData;
        }
        double global_ans = vector->getComm().maxReduce(local_ans);
        if ( global_ans == infNorm )
          utils->passes ( "Inf norm" );
        else
          utils->failure ( "Inf norm" );
}


template <typename VECTOR_FACTORY>
void ScaleVector( AMP::UnitTest *utils )
{
        AMP::LinearAlgebra::Vector::shared_ptr  vector1 ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vector2 ( VECTOR_FACTORY::getVector() );
        double beta = 1.2345;
        vector2->setRandomValues ();
        vector1->scale ( beta , vector2 );
        bool pass = true;
        AMP::LinearAlgebra::Vector::iterator  curData1 = vector1->begin();
        AMP::LinearAlgebra::Vector::iterator  endData1 = vector1->end();
        AMP::LinearAlgebra::Vector::iterator  curData2 = vector2->begin();
        while ( curData1 != endData1 )
        {
          if ( *curData1 != beta * *curData2 )
            pass = false;
          ++curData1;
          ++curData2;
        }
        if ( pass )
          utils->passes ( "scale vector 1" );
        else
          utils->failure ( "scale vector 1" );
        vector2->scale ( beta );
        vector1->subtract ( vector2 , vector1 );
        if ( vector1->maxNorm() < 0.0000001 )
          utils->passes ( "scale vector 2" );
        else
          utils->failure ( "scale vector 2" );
}


#ifdef USE_EXT_PETSC
template <typename VECTOR_FACTORY>
void Bug_491( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr  vector1 ( VECTOR_FACTORY::getVector() );
    vector1->setRandomValues ();
    AMP::LinearAlgebra::Vector::shared_ptr  managed_petsc = AMP::LinearAlgebra::PetscVector::view ( vector1 );
    Vec managed_vec = managed_petsc->castTo<AMP::LinearAlgebra::PetscVector>().getVec();

    double n1, n2, ninf;
    double sp_n1 , sp_n2 , sp_inf;

    // This sets the petsc cache
    VecNormBegin ( managed_vec , NORM_1 , &n1 );
    VecNormBegin ( managed_vec , NORM_2 , &n2 );
    VecNormBegin  ( managed_vec , NORM_INFINITY , &ninf );
    VecNormEnd ( managed_vec , NORM_1 , &n1 );
    VecNormEnd ( managed_vec , NORM_2 , &n2 );
    VecNormEnd ( managed_vec , NORM_INFINITY , &ninf );
    VecNorm ( managed_vec , NORM_1 , &n1 );
    VecNorm ( managed_vec , NORM_2 , &n2 );
    VecNorm ( managed_vec , NORM_INFINITY , &ninf );

    // Now, we perform some math on vector1
    vector1->scale ( 100000 );
    sp_n1 = vector1->L1Norm();
    sp_n2 = vector1->L2Norm();
    sp_inf = vector1->maxNorm();

    // Check to see if petsc cache has been invalidated
    VecNormBegin ( managed_vec , NORM_1 , &n1 );
    VecNormBegin ( managed_vec , NORM_2 , &n2 );
    VecNormBegin  ( managed_vec , NORM_INFINITY , &ninf );
    VecNormEnd ( managed_vec , NORM_1 , &n1 );
    VecNormEnd ( managed_vec , NORM_2 , &n2 );
    VecNormEnd ( managed_vec , NORM_INFINITY ,&ninf );

    if ( fabs ( n1 - sp_n1 ) < 0.00000001*n1 )
        utils->passes ( "L1 norm -- Petsc interface begin/end" );
    else
        utils->failure ( "l1 norm -- Petsc interface begin/end" );
    if ( fabs ( n2 - sp_n2 ) < 0.00000001*n1 )
        utils->passes ( "L2 norm -- Petsc interface begin/end" );
    else
        utils->failure ( "l2 norm -- Petsc interface begin/end" );
    if ( fabs ( ninf - sp_inf ) < 0.00000001*n1 )
        utils->passes ( "Linf norm -- Petsc interface begin/end" );
    else
        utils->failure ( "Linf norm -- Petsc interface begin/end" );

    VecNorm ( managed_vec , NORM_1 , &n1 );
    VecNorm ( managed_vec , NORM_2 , &n2 );
    VecNorm ( managed_vec , NORM_INFINITY , &ninf );

    if ( fabs ( n1 - vector1->L1Norm() ) < 0.00000001*n1 )
        utils->passes ( "L1 norm -- Petsc interface begin/end" );
    else
        utils->failure ( "l1 norm -- Petsc interface begin/end" );
    if ( fabs ( n2 - vector1->L2Norm() ) < 0.00000001*n1 )
        utils->passes ( "L2 norm -- Petsc interface begin/end" );
    else
        utils->failure ( "l2 norm -- Petsc interface begin/end" );
    if ( fabs ( ninf - vector1->maxNorm() ) < 0.00000001*n1 )
        utils->passes ( "inf norm -- Petsc interface begin/end" );
    else
        utils->failure ( "inf norm -- Petsc interface begin/end" );
}
#endif


template <typename VECTOR_FACTORY>
void AddVector( AMP::UnitTest *utils )
{
        AMP::LinearAlgebra::Vector::shared_ptr  vector1 ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vector2 ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vector3 ( VECTOR_FACTORY::getVector() );
        vector1->setRandomValues ();
        vector2->setRandomValues ();
        vector3->add ( vector1 , vector2 );
        bool pass = true;
        AMP::LinearAlgebra::Vector::iterator curData1 = vector1->begin();
        AMP::LinearAlgebra::Vector::iterator endData1 = vector1->end();
        AMP::LinearAlgebra::Vector::iterator curData2 = vector2->begin();
        AMP::LinearAlgebra::Vector::iterator curData3 = vector3->begin();
        while ( curData1 != endData1 )
        {
          if ( *curData3 != *curData1 + *curData2 )
            pass = false;
          ++curData1;
          ++curData2;
          ++curData3;
        }

        if ( pass )
          utils->passes ( "add vector" );
        else
          utils->failure ( "add vector" );
}


template <typename VECTOR_FACTORY>
void SubtractVector( AMP::UnitTest *utils )
{
        AMP::LinearAlgebra::Vector::shared_ptr  vector1 ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vector2 ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vector3 ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vector4 ( VECTOR_FACTORY::getVector() );
        vector1->setRandomValues ();
        vector2->setRandomValues ();
        vector3->subtract ( vector1 , vector2 );
        bool pass = true;
        AMP::LinearAlgebra::Vector::iterator curData1 = vector1->begin();
        AMP::LinearAlgebra::Vector::iterator endData1 = vector1->end();
        AMP::LinearAlgebra::Vector::iterator curData2 = vector2->begin();
        AMP::LinearAlgebra::Vector::iterator curData3 = vector3->begin();
        while ( curData1 != endData1 )
        {
          if ( *curData3 != *curData1 - *curData2 )
            pass = false;
          ++curData1;
          ++curData2;
          ++curData3;
        }
        if ( pass )
          utils->passes ( "vector subtract 1" );
        else
          utils->failure ( "vector subtract 1" );
        vector2->scale ( -1. );
        vector4->add ( vector1 , vector2 );
        vector4->subtract ( vector3 , vector4 );
        if ( vector4->maxNorm() < 0.0000001 )
          utils->passes ( "vector subtract 2" );
        else
          utils->passes ( "vector subtract 2" );
}


template <typename VECTOR_FACTORY>
void MultiplyVector( AMP::UnitTest *utils )
{
        AMP::LinearAlgebra::Vector::shared_ptr  vector1 ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vector2 ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vector3 ( VECTOR_FACTORY::getVector() );
        vector1->setRandomValues ();
        vector2->setToScalar ( 3. );
        vector3->multiply ( vector1 , vector2 );
        bool pass = true;
        AMP::LinearAlgebra::Vector::iterator curData1 = vector1->begin();
        AMP::LinearAlgebra::Vector::iterator endData1 = vector1->end();
        AMP::LinearAlgebra::Vector::iterator curData2 = vector2->begin();
        AMP::LinearAlgebra::Vector::iterator curData3 = vector3->begin();
        while ( curData1 != endData1 )
        {
          if ( *curData3 != *curData1 * *curData2 )
            pass = false;
          ++curData1;
          ++curData2;
          ++curData3;
        }
        if ( pass )
          utils->passes  ( "vector::multiply" );
        else
          utils->failure ( "vector::multiply" );
}


template <typename VECTOR_FACTORY>
void DivideVector( AMP::UnitTest *utils )
{
        AMP::LinearAlgebra::Vector::shared_ptr  vector1 ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vector2 ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vector3 ( VECTOR_FACTORY::getVector() );
        vector1->setRandomValues ();
        vector2->setRandomValues ();
        vector3->divide ( vector1 , vector2 );
        bool pass = true;
        AMP::LinearAlgebra::Vector::iterator curVal1 = vector1->begin();
        AMP::LinearAlgebra::Vector::iterator curVal2 = vector2->begin();
        AMP::LinearAlgebra::Vector::iterator curVal3 = vector3->begin();
        AMP::LinearAlgebra::Vector::iterator endVal3 = vector3->end();
        while ( curVal3 != endVal3 )
        {
          if ( *curVal3 != *curVal1 / *curVal2 )
            pass = false;
          ++curVal1;
          ++curVal2;
          ++curVal3;
        }
        if ( pass )
          utils->passes ( "vector::divide" );
        else
          utils->failure ( "vector::divide" );

        //if ( utils->rank() == 2 )
        //{
        //  std::cout << vector2 << std::endl;
        //}
}


template <typename VECTOR_FACTORY>
void VectorIteratorLengthTest( AMP::UnitTest *utils )
{
        AMP::LinearAlgebra::Vector::shared_ptr  vector1 ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::iterator  curEntry = vector1->begin();
        AMP::LinearAlgebra::Vector::iterator  endEntry = vector1->end();
        size_t i = 0;
        while ( curEntry != endEntry )
        {
          i++;
          ++curEntry;
        }
        size_t k = vector1->getLocalSize();
        if ( i == k )
          utils->passes ( "Iterated over the correct number of entries" );
        else
          utils->failure ( "Wrong number of entries in iterator" );
}

template <typename ITERATOR>
void both_VectorIteratorTests ( AMP::LinearAlgebra::Vector::shared_ptr p , AMP::UnitTest *utils )
{
        typename ITERATOR::vector_type  &ref = p->castTo<typename ITERATOR::vector_type> ();

        int kk = p->getLocalSize();
        if ( ( p->end() - p->begin() ) == (int)p->getLocalSize() )
          utils->passes ( "Subtracting begin from end" );
        else
          utils->failure ( "Subtracting begin from end" );

        if ( (int)( p->begin() - p->end() ) == -(int)p->getLocalSize() )
          utils->passes ( "Subtracting end from beginning" );
        else
          utils->failure ( "Subtracting end from beginning" );

        ITERATOR cur1 , cur2;
        cur1 = cur2 = ref.begin();
        ITERATOR end = ref.end();
        ++cur1;
        ++cur2;
        int i = 0;
        while ( cur2 != end )
        {
          if ( i == 10 ) break;
          ++cur2;
          i++;
        }
        int tt = (cur2 - cur1);
        if ( i == tt )
          utils->passes ( "Subtracting arbitrary iterators" );
        else
          utils->failure ( "Subtracting arbitrary iterators" );

        p->setToScalar ( 5.0 );
        i = 0;
        for ( cur1 = ref.begin() ; cur1 != end ; cur1++ )
        {
          if ( (*cur1) != 5.0 )
            break;
          i++;
        }
        if ( i == (int)p->getLocalSize() )
          utils->passes ( "Iterating data access" );
        else
          utils->failure ( "Iterating data access" );

        cur1 = end;
        i = 0;
        do
        {
          --cur1;
          if ( (*cur1) != 5.0 )
            break;
          i++;
        }
        while ( cur1 != ref.begin() );

        if ( i == kk )
          utils->passes ( "Iterating backward data access" );
        else
          utils->failure ( "Iterating backward data access" );

        if ( p->getLocalSize() > 7 )
        {
          cur1 = ref.begin();
          cur2 = cur1 + 5;
          if ( (cur2 - cur1) == 5 )
            utils->passes ( "Adding and subtracting" );
          else
            utils->failure ( "Adding and subtracting" );
          i = 0;
          while ( cur2 != end )
          {
            i++;
            ++cur2;
          }
          if ( i == ( (int)p->getLocalSize() - 5 ) )
            utils->passes ( "Adding and iterating" );
          else
            utils->failure ( "Adding and iterating" );

          cur1 += 5;
          i = 0;
          while ( cur1 != end )
          {
            i++;
            ++cur1;
          }
          if ( i == ( (int)p->getLocalSize() - 5 ) )
            utils->passes ( "Add-equal and iterating" );
          else
            utils->failure ( "Add-equal and iterating" );
        }
}


template <typename VECTOR_FACTORY>
void VectorIteratorTests( AMP::UnitTest *utils )
{
      AMP::LinearAlgebra::Vector::shared_ptr  vector1 ( VECTOR_FACTORY::getVector() );
      both_VectorIteratorTests<AMP::LinearAlgebra::Vector::iterator> ( vector1, utils );
      both_VectorIteratorTests<AMP::LinearAlgebra::Vector::const_iterator> ( vector1, utils );
}


template <typename VECTOR_FACTORY>
void VerifyVectorMin( AMP::UnitTest *utils )
{
        AMP::LinearAlgebra::Vector::shared_ptr  vec ( VECTOR_FACTORY::getVector() );
        vec->setRandomValues ();
        vec->scale ( -1.0 );           // make negative
        if ( fabs ( vec->min() + vec->maxNorm() ) < 1.e-10 )
          utils->passes ( "minimum of negative vector == ||.||_infty" );
        else
          utils->failure ( "minimum of negative vector != ||.||_infty" );
}


template <typename VECTOR_FACTORY>
void VerifyVectorMax( AMP::UnitTest *utils )
{
        AMP::LinearAlgebra::Vector::shared_ptr  vec ( VECTOR_FACTORY::getVector() );
        vec->setRandomValues ();
        if ( fabs ( vec->max() - vec->maxNorm() ) < 1.e-10 )
          utils->passes ( "maximum of positive vector == ||.||_infty" );
        else
          utils->failure ( "maximum of positive vector != ||.||_infty" );
}


template <typename VECTOR_FACTORY>
void VerifyVectorMaxMin( AMP::UnitTest *utils )
{
    AMP::LinearAlgebra::Vector::shared_ptr  vec ( VECTOR_FACTORY::getVector() );
    bool passes = true;
    for ( size_t i=0; i!=10; i++) {
        vec->setRandomValues ();
        vec->addScalar ( vec , -0.5 );
        vec->scale ( 2.0 );  // vec i.i.d [-1,1);
        double  max = vec->max ();
        double  min = vec->min ();
        double  ans = std::max ( fabs ( max ) , fabs ( min ) );
        if ( fabs( ans - vec->maxNorm() ) >= 1.e-20 )
            passes = false;
    }
    if ( passes )
        utils->passes ( "Max and min correctly predict maxNorm()" );
    else
        utils->failure ( "Max and min fail to predict maxNorm()" );
}


template <typename VECTOR_FACTORY>
class SetRandomValuesVector
{
    public:
      static const char * get_test_name () { return "vector::setRandomValues"; }

      static  void verify_vector ( AMP::UnitTest *utils , AMP::LinearAlgebra::Vector::shared_ptr v )
      {
        if ( v->min() >= 0 )
          utils->passes ( "Min value >= 0 " );
        else
          utils->failure ( "Min value < 0 " );
        if ( v->max() < 1 )
          utils->passes ( "Max value < 1" );
        else
          utils->failure ( "Max value >= 1" );
        if ( v->L2Norm() > 0 )
          utils->passes ( "Non-zero vector created" );
        else
          utils->failure ( "Zero vector created" );
      }

      static void run_test( AMP::UnitTest *utils )
      {
        AMP::LinearAlgebra::Vector::shared_ptr  vector ( VECTOR_FACTORY::getVector() );
        double l2norm = -1;
        for ( size_t i = 0 ; i != 5 ; i++ )
        {
          vector->setRandomValues ();
          if ( fabs ( l2norm - vector->L2Norm() ) > 0.00001 )
            utils->passes ( "Distinct vector created" );
          else
            utils->failure ( "Similar vector created" );
          l2norm = vector->L2Norm();
          verify_vector ( utils , vector );
        }

      }
};


template <typename VECTOR_FACTORY>
void ReciprocalVector( AMP::UnitTest *utils )
{
        AMP::LinearAlgebra::Vector::shared_ptr  vectora ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vectorb ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vectorc ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vectord ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vector1 ( VECTOR_FACTORY::getVector() );
        vectora->setRandomValues ();
        vectorb->reciprocal ( vectora );
        vector1->setToScalar ( 1. );
        vectorc->divide ( vector1 , vectora );
        vectord->subtract ( vectorb , vectorc );
        if ( vectord->maxNorm() < 0.0000001 )
        {
          utils->passes ("vector::reciprocal");
        }
        else
        {
          utils->failure ("vector::reciprocal");
        }
}


template <typename VECTOR_FACTORY>
class LinearSumVector
{
    public:
      static const char * get_test_name () { return "vector::linearSum"; }

      static  void do_instance ( AMP::UnitTest *utils , double alpha , double beta , const char *msg )
      {
        AMP::LinearAlgebra::Vector::shared_ptr  vectora ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vectorb ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vectorc ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vectord ( VECTOR_FACTORY::getVector() );
        vectora->setRandomValues ();
        vectorb->setRandomValues ();
        vectorc->linearSum ( alpha , vectora , beta , vectorb );
        vectora->scale ( alpha );
        vectorb->scale ( beta );
        vectord->add ( vectora , vectorb );
        vectord->subtract ( vectorc , vectord );
        if ( vectord->maxNorm() < 0.0000001 )
          utils->passes ( msg );
        else
          utils->failure ( msg );
      }

      static void run_test( AMP::UnitTest *utils )
      {
        do_instance ( utils , 1.2345 , 9.8765 , "linear sum 1" );
        do_instance ( utils , -1.2345 , 9.8765 , "linear sum 2" );
        do_instance ( utils , 1.2345 , -9.8765 , "linear sum 3" );
        do_instance ( utils , -1.2345 , -9.8765 , "linear sum 4" );
      }
};


template <typename VECTOR_FACTORY>
class AxpyVector
{
    public:
      static const char * get_test_name () { return "vector::axpy"; }

      static void do_instance ( AMP::UnitTest *utils , double alpha , const char *msg )
      {
        AMP::LinearAlgebra::Vector::shared_ptr  vectora ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vectorb ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vectorc ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vectord ( VECTOR_FACTORY::getVector() );
        vectora->setRandomValues ();
        vectorb->setRandomValues ();
        vectorc->linearSum ( alpha , vectora , 1. , vectorb );
        vectord->axpy ( alpha , vectora , vectorb );
        vectorc->subtract ( vectorc , vectord );
        if ( vectorc->maxNorm() < 0.0000001 )
          utils->passes ( msg );
        else
          utils->failure ( msg );
      }

      static void  run_test ( AMP::UnitTest *utils )
      {
        do_instance ( utils , 6.38295 , "axpy 1" );
        do_instance ( utils , -6.38295 , "axpy 2" );
    }
};


template <typename VECTOR_FACTORY>
class AxpbyVector
{
    public:
      static const char * get_test_name () { return "vector::axpby"; }

      static void do_instance ( AMP::UnitTest *utils , double alpha , double beta , const char *msg )
      {
        AMP::LinearAlgebra::Vector::shared_ptr  vectora ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vectorb ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vectorb1 ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vectorc ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vectord ( VECTOR_FACTORY::getVector() );
        vectora->setRandomValues ();
        vectorb->setRandomValues ();
        vectorc->copyVector ( vectorb );
        vectorb1->linearSum ( alpha , vectora , beta , vectorb );
        vectorb->linearSum ( alpha , vectora , beta , vectorb );
        vectorc->axpby ( alpha , beta , vectora );
        vectord->subtract ( vectorc , vectorb1 );
        double maxNorm = vectord->maxNorm();
        if ( maxNorm < 0.0000001 )
          utils->passes ( msg );
        else
          utils->failure ( msg );
        vectord->subtract ( vectorc , vectorb );
        maxNorm = vectord->maxNorm();
        if ( maxNorm < 0.0000001 )
          utils->passes ( msg );
        else
          utils->failure ( msg );
      }

      static void  run_test ( AMP::UnitTest *utils )
      {
        do_instance ( utils , 6.38295 , 99.273624 , "axpby 1" );
        do_instance ( utils , 6.38295 , -99.273624 , "axpby 2" );
        do_instance ( utils , -6.38295 , 99.273624 , "axpby 3" );
        do_instance ( utils , -6.38295 , -99.273624 , "axpby 4" );
    }
};


template <typename VECTOR_FACTORY>
class CopyVector
{
    public:
      static const char * get_test_name () { return "vector::copyVector"; }

      static void  do_instance ( AMP::UnitTest *utils , const char *msg1 , const char *msg2 )
      {
        AMP::LinearAlgebra::Vector::shared_ptr  vectora ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vectorb ( VECTOR_FACTORY::getVector() );
        AMP::LinearAlgebra::Vector::shared_ptr  vectorc ( VECTOR_FACTORY::getVector() );

        vectora->setRandomValues ();
        vectorb->copyVector ( vectora );
        vectorc->subtract ( vectora , vectorb );
        if ( vectorc->maxNorm() < 0.0000001 )
          utils->passes ( msg1 );
        else
          utils->failure ( msg1 );
        vectora->scale ( 100. );
        vectorc->subtract ( vectora , vectorb );

        double c_maxNorm = vectorc->maxNorm();
        double b_maxNorm = vectorb->maxNorm();
        if ( c_maxNorm >= 98.999 * b_maxNorm )
          utils->passes ( msg2 );
        else
          utils->failure ( msg2 );
      }

      static void  run_test ( AMP::UnitTest *utils )
      {
        do_instance ( utils , "copy vector 1" , "copy vector 2" );
        do_instance ( utils , "copy vector 3" , "copy vector 4" );
        do_instance ( utils , "copy vector 5" , "copy vector 6" );
        do_instance ( utils , "copy vector 7" , "copy vector 8" );
    }
};


template <typename VECTOR_FACTORY>
void VerifyVectorGhostCreate( AMP::UnitTest *utils )
      {
        AMP::LinearAlgebra::Vector::shared_ptr  vector = VECTOR_FACTORY::getVector();
        int num_ghosts = vector->getGhostSize();
        AMP_MPI globalComm = AMP_MPI(AMP_COMM_WORLD);
        num_ghosts = globalComm.sumReduce( num_ghosts );

        if ( utils->size()==1 )
          utils->expected_failure ("No ghost cells for single processor");
        else if ( num_ghosts > 0 )
          utils->passes ("verify ghosts created");
        else
          utils->failure ("verify ghosts created");

}


template <typename VECTOR_FACTORY>
void VerifyVectorMakeConsistentAdd( AMP::UnitTest *utils )
{
    AMP::Discretization::DOFManager::shared_ptr  dofmap = VECTOR_FACTORY::getDOFMap();
    AMP::LinearAlgebra::Vector::shared_ptr  vector = VECTOR_FACTORY::getVector();
    if ( !vector )
        utils->failure ( "verify makeConsistent () for add" );

    // Zero the vector
    vector->zero();
    if ( vector->getUpdateStatus() != AMP::LinearAlgebra::Vector::UNCHANGED )
        utils->failure ( "zero leaves vector in UNCHANGED state" );
    
    // Set and add local values by global id (this should not interfer with the add)
    for (size_t i = dofmap->beginDOF() ; i != dofmap->endDOF() ; i++ ) {
        vector->setLocalValueByGlobalID( i, 0.0 );
        vector->addLocalValueByGlobalID( i, 0.0 );
    }
    if ( vector->getUpdateStatus() != AMP::LinearAlgebra::Vector::LOCAL_CHANGED )
        utils->failure ( "local set/add leaves vector in LOCAL_CHANGED state" );

    // Add values by global id
    for (size_t i = dofmap->beginDOF() ; i != dofmap->endDOF() ; i++ )
        vector->addValueByGlobalID ( i , (double) i );
    if ( vector->getUpdateStatus() != AMP::LinearAlgebra::Vector::ADDING )
        utils->failure ( "addValueByGlobalID leaves vector in ADDING state" );

    double offset = (double) (1 << utils->rank() );
    for ( size_t i = 0 ; i != vector->getGhostSize() ; i++ ) {
        size_t ndx = vector->getCommunicationList()->getGhostIDList()[i];
        vector->addValueByGlobalID ( ndx , offset );
    }

    // Perform a makeConsistent ADD and check the result
    vector->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_ADD );
    if ( vector->getUpdateStatus() != AMP::LinearAlgebra::Vector::UNCHANGED )
        utils->failure ( "makeConsistent leaves vector in UNCHANGED state" );
    std::map<int,std::set<size_t> >  ghosted_entities;
    for ( size_t i = dofmap->beginDOF() ; i != dofmap->endDOF() ; i++ ) {
        double diff_double = fabs ( vector->getValueByGlobalID ( i ) - (double)i );
        if ( diff_double > 0.00001 ) {
            int ioffset = lround(diff_double);
            int cur_rank = 0;
            while ( ioffset > 0 ) {
                if ( ioffset & 1 )
                    ghosted_entities[cur_rank].insert ( i );
                ioffset >>= 1;
                cur_rank++;
            }
        }
    }
    std::vector<size_t>::const_iterator  cur_replicated = vector->getCommunicationList()->getReplicatedIDList().begin();
    std::vector<size_t>::const_iterator  end_replicated = vector->getCommunicationList()->getReplicatedIDList().end();
    while ( cur_replicated != end_replicated ) {
        bool found = false;
        for ( int i = 0 ; i != utils->size() ; i++ ) {
            std::set<size_t>::iterator  location = ghosted_entities[i].find ( *cur_replicated );
            if ( location != ghosted_entities[i].end() ) {
               found = true;
                ghosted_entities[i].erase ( location );
                break;
            }
        }
        if ( !found ) {
            utils->failure ( "overly ghosted value" );
            return;
        }
        ++cur_replicated;
    }
    size_t last_size = 0;
    for ( int i = 0 ; i != utils->size() ; i++ )
        last_size += ghosted_entities[i].size();
    if ( last_size == 0 )
        utils->passes ( "all ghosted values accounted for" );
    else
        utils->failure ( "some ghosted values not set" );
}


template <typename VECTOR_FACTORY>
void VerifyVectorMakeConsistentSet( AMP::UnitTest *utils )
{
    AMP::Discretization::DOFManager::shared_ptr  dofmap = VECTOR_FACTORY::getDOFMap();
    AMP::LinearAlgebra::Vector::shared_ptr  vector = VECTOR_FACTORY::getVector();

    // Zero the vector
    vector->zero();
    if ( vector->getUpdateStatus() != AMP::LinearAlgebra::Vector::UNCHANGED )
        utils->failure ( "zero leaves vector in UNCHANGED state" );
    
    // Set and add local values by global id (this should not interfer with the add)
    for (size_t i = dofmap->beginDOF() ; i != dofmap->endDOF() ; i++ ) {
        vector->setLocalValueByGlobalID( i, 0.0 );
        vector->addLocalValueByGlobalID( i, 0.0 );
    }
    if ( vector->getUpdateStatus() != AMP::LinearAlgebra::Vector::LOCAL_CHANGED )
        utils->failure ( "local set/add leaves vector in LOCAL_CHANGED state" );

    // Set values by global id
    for (size_t i = dofmap->beginDOF() ; i != dofmap->endDOF() ; i++ )
        vector->setValueByGlobalID ( i , (double) i );
    if ( vector->getUpdateStatus()!=AMP::LinearAlgebra::Vector::LOCAL_CHANGED &&
         vector->getUpdateStatus()!=AMP::LinearAlgebra::Vector::SETTING )
        utils->failure ( "setValueByGlobalID leaves vector in SETTING or LOCAL_CHANGED state" );

    // Perform a makeConsistent SET and check the result
    vector->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    if ( vector->getUpdateStatus() != AMP::LinearAlgebra::Vector::UNCHANGED )
        utils->failure ( "makeConsistent leaves vector in UNCHANGED state" );
    if ( vector->getGhostSize() > 0 ) {
        AMP::LinearAlgebra::CommunicationList::shared_ptr comm_list = vector->getCommunicationList();
        std::vector<double> ghostList ( vector->getGhostSize() );
        std::vector<size_t> ghostIDList = comm_list->getGhostIDList();
        vector->getValuesByGlobalID ( vector->getGhostSize() , (size_t*) &(ghostIDList[0]) , &(ghostList[0]) );
        bool testPassed = true;
        for ( size_t i = 0 ; i != vector->getGhostSize() ; i++ ) {
            if ( fabs ( ghostList[i] - (double)(ghostIDList[i]) ) > 0.0000001 )
                testPassed = false;
        }
        if ( testPassed )
            utils->passes ( "ghost set correctly in vector" );
        else
            utils->failure ( "ghost not set correctly in vector" );
    }
    if ( vector->getGhostSize() > 0 ) {
        AMP::LinearAlgebra::CommunicationList::shared_ptr comm_list = vector->getCommunicationList();
        std::vector<size_t> ghostIDList = comm_list->getGhostIDList();
        bool testPassed = true;
        for ( size_t i = 0 ; i != vector->getGhostSize() ; i++ ) {
            size_t  ghostNdx = ghostIDList[i];
            double  ghostVal = vector->getValueByGlobalID ( ghostNdx );
            if ( fabs ( ghostVal - (double)ghostNdx ) > 0.0000001 )
                testPassed = false;
        }
        if ( testPassed )
            utils->passes ( "ghost set correctly in alias" );
        else
            utils->failure ( "ghost set correctly in alias" );
    }
}


}
}


/// \endcond
#endif
