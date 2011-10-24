
#include "MultiVector.h"
#include "ManagedVector.h"
#include <stdexcept>
#include <algorithm>
#include <math.h>

#include "utils/Utilities.h"


namespace AMP {
namespace LinearAlgebra {


Vector::shared_ptr  MultiVector::encapsulate ( Vector::shared_ptr &vec , AMP_MPI comm )
{
    if ( comm.isNull() )
    {
      comm = vec->getComm();
    }
    Vector::shared_ptr retval = create ( vec->getVariable()->getName() , comm );
    retval->castTo<MultiVector>().addVector ( vec );
    if ( vec->isA<DataChangeFirer>() )
    {
      vec->castTo<DataChangeFirer>().registerListener ( &(retval->castTo<DataChangeListener>()) );
    }
    return retval;
}


Vector::shared_ptr  MultiVector::view ( Vector::shared_ptr &vec , AMP_MPI comm )
{
    Vector::shared_ptr  retval = Vector::shared_ptr ();

    // Check to see if this is a multivector
    if ( vec->isA<MultiVector>() )
    {
      if ( !comm.isNull() )
      {
        AMP_ASSERT ( comm.compare(vec->getComm()) != 0 );
      }
      retval = vec;
    }

    // Check to see if the engine is a multivector
    if ( vec->isA<ManagedVector>() )
    {
      if ( vec->castTo<ManagedVector>().getVectorEngine()->isA<MultiVector>() )
      {
        if ( !comm.isNull() )
        {
          AMP_ASSERT ( comm.compare(vec->getComm()) != 0 );
        }
        retval = boost::dynamic_pointer_cast<Vector> ( vec->castTo<ManagedVector>().getVectorEngine() );
      }
    }

    // If still don't have a multivector, make one
    if ( !retval )
    {
      if ( comm.isNull() )
      {
        comm = vec->getComm();
      }
      retval = create ( vec->getVariable()->getName() , comm );
      MultiVector &tretval = retval->castTo<MultiVector> ();
      tretval.d_vVectors.push_back ( vec );
      tretval.d_vGlobalOffsets.push_back ( vec->getGlobalMaxID() + tretval.d_vGlobalOffsets.back() );
      tretval.d_vLocalOffsets.push_back ( vec->getLocalMaxID() + tretval.d_vLocalOffsets.back() );
      if ( vec->isA<DataChangeFirer>() )
      {
        vec->castTo<DataChangeFirer>().registerListener ( &tretval );
      }
//      vec = retval;  This line of code will be important one day...
    }

    return retval;
}


bool MultiVector::containsPointer ( const Vector::shared_ptr p ) const
{
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
    {
      if ( d_vVectors[i].get() == p.get() )
      {
        return true;
      }
    }
    return false;
}


void MultiVector::selectInto ( const VectorSelector &s , Vector::shared_ptr retVal )
{
    vector_iterator  cur = beginVector();
    while ( cur != endVector() )
    {
      (*cur)->selectInto ( s , retVal );
      cur++;
    }
}

void MultiVector::subtract ( const VectorOperations &x , const VectorOperations &y )
{
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
      d_vVectors[i]->subtract ( getVector ( x , i ) , getVector ( y , i ) );
}


void MultiVector::multiply ( const VectorOperations &x , const VectorOperations &y )
{
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
      d_vVectors[i]->multiply ( getVector ( x , i ) , getVector ( y , i ) );
}


void MultiVector::divide ( const VectorOperations &x , const VectorOperations &y )
{
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
      d_vVectors[i]->divide ( getVector ( x , i ) , getVector ( y , i ) );
}


void MultiVector::reciprocal ( const VectorOperations &x )
{
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
      d_vVectors[i]->reciprocal ( getVector ( x , i ) );
}


void MultiVector::linearSum(double alpha, const VectorOperations &x, double beta, const VectorOperations &y)
{
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
      d_vVectors[i]->linearSum ( alpha , getVector ( x , i ) , beta , getVector ( y , i ) );
}


void MultiVector::axpy(double alpha, const VectorOperations &x, const VectorOperations &y)
{
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
      d_vVectors[i]->axpy ( alpha , getVector ( x , i ) , getVector ( y , i ) );
}


void MultiVector::axpby(double alpha, double beta, const VectorOperations &x)
{
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
      d_vVectors[i]->axpby ( alpha , beta , getVector ( x , i ) );
}


void MultiVector::abs ( const VectorOperations &x )
{
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
      d_vVectors[i]->abs ( getVector ( x , i ) );
}


double MultiVector::min(void) const
{
    double retVal = 0.0;
    if ( d_vVectors.size() > 0 )
    {
      retVal = d_vVectors[0]->min();
      for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
        retVal = std::min ( retVal , d_vVectors[i]->min () );
    }
    return retVal;
}


double MultiVector::max(void) const
{
    double retVal = 0.0;
    if ( d_vVectors.size() > 0 )
    {
      retVal = d_vVectors[0]->max();
      for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
        retVal = std::max ( retVal , d_vVectors[i]->max () );
    }
    return retVal;
}


void MultiVector::setRandomValues ()
{
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
      d_vVectors[i]->setRandomValues ();
}


void MultiVector::makeConsistent ( ScatterType t )
{
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
      d_vVectors[i]->makeConsistent ( t );
}


void MultiVector::assemble ()
{
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
      d_vVectors[i]->assemble ();
}


double MultiVector::L1Norm () const
{
    double ans = 0.0;
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
      ans += d_vVectors[i]->L1Norm ();
    return ans;
}


double MultiVector::L2Norm () const
{
    double ans = 0.0;
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
    {
      double t = d_vVectors[i]->L2Norm();
      ans += t*t;
    }
    return sqrt ( ans );
}


double MultiVector::maxNorm () const
{
    double ans = 0.0;
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
      ans = std::max ( ans , d_vVectors[i]->maxNorm () );
    return ans;
}


double MultiVector::dot ( const VectorOperations &rhs ) const
{
    double ans = 0.0;
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
      ans += d_vVectors[i]->dot ( getVector ( rhs , i ) );
    return ans;
}


size_t MultiVector::getLocalSize () const
{
    size_t ans = 0;
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
      ans += d_vVectors[i]->getLocalSize();
    return ans;
}


size_t MultiVector::getGlobalSize () const
{
    size_t ans = 0;
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
      ans += d_vVectors[i]->getGlobalSize();
    return ans;
}


size_t MultiVector::getGhostSize () const
{
    size_t ans = 0;
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
      ans += d_vVectors[i]->getGhostSize();
    return ans;
}


void  MultiVector::putRawData ( double *in )
{
    int cur_off = 0;
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
    {
      d_vVectors[i]->putRawData ( in + cur_off );
      cur_off += d_vVectors[i]->getLocalSize();
    }
}


/*
void MultiVector::copyVector ( const Vector &src )
{
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
      d_vVectors[i]->copyVector ( getVector ( src , i ) );
}
*/


void MultiVector::swapVectors(Vector &other)
{
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
      d_vVectors[i]->swapVectors ( getVector ( other , i ) );
}


void MultiVector::aliasVector (Vector &other)
{
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
      d_vVectors[i]->aliasVector ( getVector ( other , i ) );
}


void MultiVector::setToScalar ( double alpha )
{
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
      d_vVectors[i]->setToScalar ( alpha );
}


void MultiVector::scale ( double alpha , const VectorOperations &x )
{
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
      d_vVectors[i]->scale ( alpha , getVector ( x , i ) );
}


void MultiVector::scale ( double alpha )
{
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
      d_vVectors[i]->scale ( alpha );
}


void MultiVector::add ( const VectorOperations &x , const VectorOperations &y )
{
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
      d_vVectors[i]->add ( getVector ( x , i ) , getVector ( y , i ) );
}


VectorEngine::BufferPtr  MultiVector::getNewBuffer()
{
    return VectorEngine::BufferPtr ();
}


const void              *MultiVector::getDataBlock ( size_t i ) const
{
    return getRawDataBlockAsVoid ( i );
}


void                    *MultiVector::getDataBlock ( size_t i )
{
    return getRawDataBlockAsVoid ( i );
}


bool                     MultiVector::sameEngine ( VectorEngine &rhs ) const
{
    return rhs.isA<MultiVector>();
}


void                     MultiVector::swapEngines ( VectorEngine::shared_ptr p )
{
    return swapVectors ( p->castTo<MultiVector>() );
}


VectorEngine::shared_ptr MultiVector::cloneEngine ( VectorEngine::BufferPtr  ) const
{
    return boost::dynamic_pointer_cast<VectorEngine> ( Vector::cloneVector ( "engine_clone" ) );
}


size_t MultiVector::numberOfDataBlocks () const
{
    size_t ans = 0;
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
      ans += d_vVectors[i]->numberOfDataBlocks();
    return ans;
}


Vector::shared_ptr  MultiVector::subsetVectorForVariable ( const Variable::shared_ptr  &name )
{
    Vector::shared_ptr retVal_1;

    AMP_ASSERT ( name );
    retVal_1 = Vector::subsetVectorForVariable ( name );
    if ( !retVal_1 )
    {
      for ( unsigned int i = 0 ; i != d_vVectors.size() ; i++ )
      {
        retVal_1 = d_vVectors[i]->subsetVectorForVariable ( name );
        if ( retVal_1 )
        {
          break;
        }
      }
    }
    if ( !retVal_1 )
    {
      if ( name->isA<MultiVariable> () )
      {
        retVal_1 = create ( name , getComm() );
        MultiVector  &out_vec = retVal_1->castTo<MultiVector> ();
        MultiVariable &in_var = name->castTo<MultiVariable> ();
        for ( size_t i = 0 ; i != in_var.numVariables () ; i++ )
        {
          Vector::shared_ptr  t = subsetVectorForVariable ( in_var.getVariable ( i ) );
          if ( !t )
          {
            retVal_1.reset ();
            break;
          }
          out_vec.addVector ( t );
        }
      }
    }

    return retVal_1;
}


size_t MultiVector::sizeOfDataBlock ( size_t i ) const
{
    size_t retVal = 0;
    size_t rightOffset , leftOffset;
    rightOffset = leftOffset = 0;
    for ( size_t j = 0 ; j != d_vVectors.size() ; j++ )
    {
      rightOffset += d_vVectors[j]->numberOfDataBlocks();
      if ( i < rightOffset )
      {
        retVal = d_vVectors[j]->sizeOfDataBlock ( i - leftOffset );
        break;
      }
     leftOffset = rightOffset;
    }
    return retVal;
}


void  MultiVector::copyOutRawData ( double ** out )
{
    *out = new double [ getLocalSize() ];
    size_t curOffset = 0;
    for ( size_t j = 0 ; j != d_vVectors.size() ; j++ )
    {
      curOffset += d_vVectors[j]->getLocalSize();
      d_vVectors[j]->copyOutRawData ( out + curOffset );
    }
}


void MultiVector::addVector ( Vector::shared_ptr  v )
{
    Vector::shared_ptr  in_vec = view ( v );

    for ( size_t i = 0 ; i != in_vec->castTo<MultiVector>().getNumberOfSubvectors() ; i++ )
    {
      Vector::shared_ptr  curvec = in_vec->castTo<MultiVector>().getVector ( i );
      d_vVectors.push_back ( curvec );
      d_vGlobalOffsets.push_back ( curvec->getGlobalMaxID() + d_vGlobalOffsets.back() );
      d_vLocalOffsets.push_back ( curvec->getLocalMaxID() + d_vLocalOffsets.back() );
      if ( curvec->isA<DataChangeFirer>() )
      {
        curvec->castTo<DataChangeFirer>().registerListener ( this );
      }
    }
}


void * MultiVector::getRawDataBlockAsVoid ( size_t i )
{
    size_t curOffset = 0;
    for ( size_t j = 0 ; j != d_vVectors.size() ; j++ )
    {
      curOffset += d_vVectors[j]->numberOfDataBlocks();
      if ( i < curOffset )
      {
        return d_vVectors[j]->getRawDataBlock<double> ( i - curOffset + d_vVectors[j]->numberOfDataBlocks () );
      }
    }
    return 0;
}


 const void * MultiVector::getRawDataBlockAsVoid ( size_t i ) const
 {
    size_t curOffset = 0;
    for ( size_t j = 0 ; j != d_vVectors.size() ; j++ )
    {
      curOffset += d_vVectors[j]->numberOfDataBlocks();
      if ( i < curOffset )
      {
        return d_vVectors[j]->getRawDataBlock<double> ( i - curOffset + d_vVectors[j]->numberOfDataBlocks () );
      }
    }
    return 0;
 }


void MultiVector::dumpOwnedData ( std::ostream &out , size_t GIDoffset , size_t LIDoffset ) const
{
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
    {
      if ( d_vVectors[i]->getVariable() )
      {
        out << "[ " << d_vVectors[i]->getVariable()->getName() << " ]\n";
      }
      d_vVectors[i]->dumpOwnedData ( out , GIDoffset + d_vGlobalOffsets[i] , LIDoffset + d_vLocalOffsets[i] );
    }
}


void MultiVector::dumpGhostedData ( std::ostream &out , size_t offset ) const
{
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
    {
      if ( d_vVectors[i]->getVariable() )
      {
        out << "[ " << d_vVectors[i]->getVariable()->getName() << " ]\n";
      }
      d_vVectors[i]->dumpGhostedData ( out , offset + d_vGlobalOffsets[i] );
    }
}


Vector::shared_ptr MultiVector::cloneVector(const Variable::shared_ptr name) const
{
    Vector::shared_ptr  retVal = Vector::shared_ptr ( new MultiVector ( name ) );
    MultiVector  &ret = retVal->castTo<MultiVector> ();
    ret.d_Comm = d_Comm;
    ret.d_vVectors.resize ( d_vVectors.size() );
    ret.d_vGlobalOffsets.resize ( d_vGlobalOffsets.size() );
    ret.d_vLocalOffsets.resize ( d_vLocalOffsets.size() );
    for ( size_t i = 0 ; i != d_vGlobalOffsets.size() ; i++ )
    {
      ret.d_vGlobalOffsets[i] = d_vGlobalOffsets[i];
      ret.d_vLocalOffsets[i] = d_vLocalOffsets[i];
    }
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
    {
      ret.d_vVectors[i] = d_vVectors[i]->cloneVector ();
    }
    return retVal;
}


void  MultiVector::partitionValues ( int num , int *indices , const double *vals ,
                        std::vector<std::vector<int> >    &out_indices ,
                        std::vector<std::vector<double> >  &out_vals ,
                        const std::vector<size_t>    &offset_vec ,
                        std::vector<std::vector<size_t> >   *remap ) const
{
    out_indices.resize ( offset_vec.size()-1 );
    out_vals.resize ( offset_vec.size()-1 );
    std::vector<size_t>  cur_spot ( offset_vec.size()-1 );
    if ( remap )
    {
      remap->resize ( offset_vec.size()-1 );
    }

    for ( int i = 0 ; i != num ; i++ )
    {
      std::vector<size_t>::const_iterator  lower_bound = std::lower_bound ( offset_vec.begin() , offset_vec.end() , indices[i] );
      if ( indices[i] == (int)*lower_bound )  // Silly fence posts
      {
        lower_bound++;
      }
      size_t  vec_num = lower_bound - offset_vec.begin() - 1;
      cur_spot[vec_num]++;
    }

    for ( size_t i = 0 ; i != offset_vec.size()-1 ; i++ )
    {
      out_indices[i].resize ( cur_spot[i] );
      out_vals[i].resize ( cur_spot[i] );
      if ( remap )
      {
        (*remap)[i].resize ( cur_spot[i] );
      }
      cur_spot[i] = 0;
    }

    for ( int i = 0 ; i != num ; i++ )
    {
      std::vector<size_t>::const_iterator  lower_bound = std::lower_bound ( offset_vec.begin() , offset_vec.end() , indices[i] );
      if ( indices[i] == (int)*lower_bound )  // Silly fence posts
      {
        lower_bound++;
      }
      size_t  vec_num = lower_bound - offset_vec.begin() - 1;  // Who put this fence post here?
      out_indices[vec_num][cur_spot[vec_num]] = indices[i] - offset_vec[vec_num];
      out_vals[vec_num][cur_spot[vec_num]] = vals[i];
      if ( remap )
      {
        (*remap)[vec_num][cur_spot[vec_num]] = i;
      }
      cur_spot[vec_num]++;
    }

}


void MultiVector::setValuesByLocalID ( int num , int *indices , const double *in_vals )
{
    INCREMENT_COUNT("Virtual");
    std::vector<std::vector<int> >     ndxs;
    std::vector<std::vector<double> >  vals;
    partitionValues ( num , indices , in_vals , ndxs , vals , d_vLocalOffsets );
    for ( size_t i = 0 ; i != ndxs.size() ; i++ )
    {
      if ( ndxs[i].size() )
      {
        d_vVectors[i]->setValuesByLocalID ( ndxs[i].size() , &(ndxs[i][0]) , &(vals[i][0]) );
      }
    }
}


void MultiVector::setLocalValuesByGlobalID ( int num , int *indices , const double *in_vals )
{
    INCREMENT_COUNT("Virtual");
    std::vector<std::vector<int> >     ndxs;
    std::vector<std::vector<double> >  vals;
    partitionValues ( num , indices , in_vals , ndxs , vals , d_vGlobalOffsets );
    for ( size_t i = 0 ; i != ndxs.size() ; i++ )
    {
      if ( ndxs[i].size() )
      {
        d_vVectors[i]->setLocalValuesByGlobalID ( ndxs[i].size() , &(ndxs[i][0]) , &(vals[i][0]) );
      }
    }
}


void MultiVector::setValuesByGlobalID ( int num , int *indices , const double *in_vals )
{
    INCREMENT_COUNT("Virtual");
    std::vector<std::vector<int> >     ndxs;
    std::vector<std::vector<double> >  vals;
    partitionValues ( num , indices , in_vals , ndxs , vals , d_vGlobalOffsets );
    for ( size_t i = 0 ; i != ndxs.size() ; i++ )
    {
      if ( ndxs[i].size() )
      {
        d_vVectors[i]->setValuesByGlobalID ( ndxs[i].size() , &(ndxs[i][0]) , &(vals[i][0]) );
      }
    }
}


void MultiVector::addValuesByLocalID ( int num , int *indices , const double *in_vals )
{
    INCREMENT_COUNT("Virtual");
    std::vector<std::vector<int> >     ndxs;
    std::vector<std::vector<double> >  vals;
    partitionValues ( num , indices , in_vals , ndxs , vals , d_vLocalOffsets );
    for ( size_t i = 0 ; i != ndxs.size() ; i++ )
    {
      if ( ndxs[i].size() )
      {
        d_vVectors[i]->addValuesByLocalID ( ndxs[i].size() , &(ndxs[i][0]) , &(vals[i][0]) );
      }
    }
}


void MultiVector::addLocalValuesByGlobalID ( int num , int *indices , const double *in_vals )
{
    INCREMENT_COUNT("Virtual");
    std::vector<std::vector<int> >     ndxs;
    std::vector<std::vector<double> >  vals;
    partitionValues ( num , indices , in_vals , ndxs , vals , d_vGlobalOffsets );
    for ( size_t i = 0 ; i != ndxs.size() ; i++ )
    {
      if ( ndxs[i].size() )
      {
        d_vVectors[i]->addLocalValuesByGlobalID ( ndxs[i].size() , &(ndxs[i][0]) , &(vals[i][0]) );
      }
    }
}


void MultiVector::addValuesByGlobalID ( int num , int *indices , const double *in_vals )
{
    INCREMENT_COUNT("Virtual");
    std::vector<std::vector<int> >     ndxs;
    std::vector<std::vector<double> >  vals;
    partitionValues ( num , indices , in_vals , ndxs , vals , d_vGlobalOffsets );
    for ( size_t i = 0 ; i != ndxs.size() ; i++ )
    {
      if ( ndxs[i].size() )
      {
        d_vVectors[i]->addValuesByGlobalID ( ndxs[i].size() , &(ndxs[i][0]) , &(vals[i][0]) );
      }
    }
}


void MultiVector::getValuesByGlobalID ( int num , int *indices , double *out_vals ) const
{
    INCREMENT_COUNT("Virtual");
    std::vector<std::vector<int> >     ndxs;
    std::vector<std::vector<double> >  vals;
    std::vector<std::vector<size_t> >  remap;
    partitionValues ( num , indices , out_vals , ndxs , vals , d_vGlobalOffsets , &remap );

    for ( size_t i = 0 ; i != ndxs.size() ; i++ )
    {
      if ( ndxs[i].size() )
      {
        d_vVectors[i]->getValuesByGlobalID ( ndxs[i].size() , &(ndxs[i][0]) , &(vals[i][0]) );
      }
    }

    for ( size_t i = 0 ; i != remap.size() ; i++ )
    {
      for ( size_t j = 0 ; j != remap[i].size() ; j++ )
      {
        out_vals[remap[i][j]] = vals[i][j];
      }
    }
}


void MultiVector::getLocalValuesByGlobalID ( int num , int *indices , double *out_vals ) const
{
    INCREMENT_COUNT("Virtual");
    std::vector<std::vector<int> >     ndxs;
    std::vector<std::vector<double> >  vals;
    std::vector<std::vector<size_t> >  remap;
    partitionValues ( num , indices , out_vals , ndxs , vals , d_vGlobalOffsets , &remap );

    for ( size_t i = 0 ; i != ndxs.size() ; i++ )
    {
      if ( ndxs[i].size() )
      {
        d_vVectors[i]->getLocalValuesByGlobalID ( ndxs[i].size() , &(ndxs[i][0]) , &(vals[i][0]) );
      }
    }

    for ( size_t i = 0 ; i != remap.size() ; i++ )
    {
      for ( size_t j = 0 ; j != remap[i].size() ; j++ )
      {
        out_vals[remap[i][j]] = vals[i][j];
      }
    }
}


void MultiVector::getValuesByLocalID ( int num , int *indices , double *out_vals ) const
{
    INCREMENT_COUNT("Virtual");
    std::vector<std::vector<int> >     ndxs;
    std::vector<std::vector<double> >  vals;
    std::vector<std::vector<size_t> >  remap;
    partitionValues ( num , indices , out_vals , ndxs , vals , d_vLocalOffsets , &remap );

    for ( size_t i = 0 ; i != ndxs.size() ; i++ )
    {
      if ( ndxs[i].size() )
      {
        d_vVectors[i]->getValuesByLocalID ( ndxs[i].size() , &(ndxs[i][0]) , &(vals[i][0]) );
      }
    }

    for ( size_t i = 0 ; i != remap.size() ; i++ )
    {
      for ( size_t j = 0 ; j != remap[i].size() ; j++ )
      {
        out_vals[remap[i][j]] = vals[i][j];
      }
    }
}


}
}

