#include "MultiVector.h"
#include "VectorSelector.h"

namespace AMP {
namespace LinearAlgebra {

  inline
  SimpleVector::SimpleVector () : Vector ()
  {
  }

  inline
  void SimpleVector::copyVector(const Vector &src_vec)
  {
    if ( getLocalSize() != src_vec.getLocalSize() )
    {
      AMP_ERROR( "Mismatched vectors" );
    }
    std::copy ( src_vec.begin() , src_vec.end() , begin() );
  }

  inline
  Vector::shared_ptr  SimpleVector::create ( size_t localSize , Variable::shared_ptr var )
  {
    SimpleVector *retVal = new SimpleVector;
    retVal->d_Data.resize ( localSize );
    retVal->setVariable ( var );
    return Vector::shared_ptr ( retVal );
  }

  inline
  double SimpleVector::min(void) const
  {
    return *std::min_element ( d_Data.begin() , d_Data.end() );
  }

  inline
  double SimpleVector::max(void) const
  {
    return *std::max_element ( d_Data.begin() , d_Data.end() );
  }


  inline
  boost::shared_ptr<ParameterBase> SimpleVector::getParameters ()
  {
    AMP_ERROR( "Not implemented" );
    return boost::shared_ptr<ParameterBase> ();
  }

  inline
  Vector::shared_ptr SimpleVector::cloneVector(const Variable::shared_ptr name) const
  {
    return create ( d_Data.size() , name );
  }

  inline
  size_t  SimpleVector::numberOfDataBlocks () const
  {
    return 1;
  }

  inline
  size_t  SimpleVector::sizeOfDataBlock ( size_t i ) const
  {
    if ( i > 0 )
    {
      return 0;
    }
    return d_Data.size();
  }

  inline
  void SimpleVector::swapVectors(Vector &rhs)
  {
    d_Data.swap ( rhs.castTo<SimpleVector>().d_Data );
  }

  inline
  void SimpleVector::aliasVector(Vector &)
  {
    AMP_ERROR( "Not implemented" );
  }

  inline
  void SimpleVector::setValuesByLocalID ( int num , size_t *indices , const double *vals )
  {
    INCREMENT_COUNT("Virtual");
    for ( int i = 0 ; i != num ; i++ )
    {
      d_Data[indices[i]] = vals[i];
    }
  }

  inline
  void SimpleVector::setLocalValuesByGlobalID ( int , size_t * , const double * )
  {
    INCREMENT_COUNT("Virtual");
    AMP_ERROR( "Not implemented" );
  }

  inline
  void SimpleVector::addValuesByLocalID ( int num , size_t *indices , const double *vals )
  {
    INCREMENT_COUNT("Virtual");
    for ( int i = 0 ; i != num ; i++ )
    {
      d_Data[indices[i]] += vals[i];
    }
  }

  inline
  void SimpleVector::addLocalValuesByGlobalID ( int , size_t * , const double * )
  {
    INCREMENT_COUNT("Virtual");
    AMP_ERROR( "Not implemented" );
  }

  inline
  void SimpleVector::getLocalValuesByGlobalID ( int , size_t * , double * ) const
  {
    INCREMENT_COUNT("Virtual");
    AMP_ERROR( "Not implemented" );
  }

  inline
  void SimpleVector::assemble()
  {
    AMP_ERROR( "Not implemented" );
  }

  inline
  void   SimpleVector::putRawData ( double *in )
  {
    std::copy ( in , in + d_Data.size() , d_Data.begin() );
  }

  inline
  size_t SimpleVector::getLocalSize() const
  {
    return d_Data.size();
  }

  inline
  size_t SimpleVector::getGlobalSize() const
  {
    return getLocalSize();
  }

  inline
  void *SimpleVector::getRawDataBlockAsVoid ( size_t i )
  {
    if ( i != 0 )
    {
      return 0;
    }
    return &(d_Data[0]);
  }

  inline
  const void *SimpleVector::getRawDataBlockAsVoid ( size_t i ) const
  {
    if ( i != 0 )
    {
      return 0;
    }
    return &(d_Data[0]);
  }

  inline
  double &SimpleVector::operator[] ( size_t i )
  {
    AMP_ASSERT ( i < d_Data.size() );
    return d_Data[i];
  }

  inline
  double  SimpleVector::operator[] ( size_t i ) const 
  {
    AMP_ASSERT ( i < d_Data.size() );
    return d_Data[i];
  }

  inline
  void SimpleVector::resize ( size_t i )
  {
    d_Data.resize ( i );
  }

}
}
