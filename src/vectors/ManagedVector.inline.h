#include "VectorSelector.h"
#include "utils/Utilities.h"

namespace AMP {
namespace LinearAlgebra {

  inline
  VectorEngine::shared_ptr  ManagedVector::getVectorEngine()
  {
    return d_Engine;
  }

  inline
  std::string ManagedVector::type () const
  {
    if ( d_vBuffer ) return " ( managed data )";
    std::string retVal = " ( managed view of ";
    retVal += d_Engine->castTo<Vector>().type();
    retVal += " )";
    return retVal;
  }

  inline
  Vector::shared_ptr  ManagedVector::getRootVector()
  {
    if ( d_vBuffer )
      return shared_from_this();
    Vector *t = &(d_Engine->castTo<Vector>());
    if ( t->isA<ManagedVector>() )
      return t->castTo<ManagedVector>().getRootVector();
    return t->shared_from_this();
  }

  inline
  Vector::const_iterator ManagedVector::begin() const
  {
    if ( d_vBuffer )
      return Vector::begin();
    else
      return d_Engine->castTo<const Vector>().begin();
  }

  inline
  Vector::const_iterator ManagedVector::end() const
  {
    if ( d_vBuffer )
      return Vector::end();
    else
      return d_Engine->castTo<const Vector>().end();
  }

  inline
  Vector::iterator ManagedVector::begin()
  {
    if ( d_vBuffer )
      return Vector::begin();
    else
      return d_Engine->castTo<Vector>().begin();
  }

  inline
  Vector::iterator ManagedVector::end()
  {
    if ( d_vBuffer )
      return Vector::end();
    else
      return d_Engine->castTo<Vector>().end();
  }

  inline
  void ManagedVector::dataChanged ()
  {
     if ( *d_UpdateState == UNCHANGED ) 
        *d_UpdateState = LOCAL_CHANGED;
  }


  inline
  void ManagedVector::selectInto ( const VectorSelector &s , Vector::shared_ptr retVal )
  {
    if ( d_vBuffer )
    {
      Vector::selectInto ( s , retVal );
    }
    else
    {
      d_Engine->castTo<Vector>().selectInto ( s , retVal );
    }
  }


  inline
  ManagedVectorParameters::ManagedVectorParameters ()
  {
    d_CloneEngine = true;
    d_Buffer = VectorEngine::BufferPtr();
  }

  inline
  void *ManagedVector::getRawDataBlockAsVoid ( size_t i )
  {
    return d_Engine->getDataBlock ( i );
  }

  inline
  const void *ManagedVector::getRawDataBlockAsVoid ( size_t i ) const
  {
    return d_Engine->getDataBlock ( i );
  }

  inline
  void ManagedVector::addCommunicationListToParameters ( CommunicationList::shared_ptr comm )
  {
    d_pParameters->d_CommList = comm;
  }

  inline
  size_t  ManagedVector::numberOfDataBlocks () const
  {
    return d_Engine->numberOfDataBlocks();
  }

  inline
  size_t  ManagedVector::sizeOfDataBlock ( size_t i ) const
  {
    return d_Engine->sizeOfDataBlock ( i );
  }

  inline
  bool ManagedVector::isAnAliasOf ( Vector::shared_ptr rhs )
  {
    return isAnAliasOf ( *rhs );
  }

  inline
  boost::shared_ptr<ParameterBase>  ManagedVector::getParameters ()
  {
    return boost::dynamic_pointer_cast<ParameterBase> ( d_pParameters );
  }

  inline
   boost::shared_ptr<ManagedVectorParameters>  ManagedVector::getManagedVectorParameters ()
  {
    return d_pParameters;
  }

  inline
  size_t ManagedVector::getLocalSize() const
  {
    return d_Engine->getLocalSize();
  }

  inline
  size_t ManagedVector::getGlobalSize() const
  {
    return d_Engine->getGlobalSize();
  }

  inline
  void ManagedVector::putRawData ( double *in )
  {
    d_Engine->putRawData ( in );
  }

  inline
  void ManagedVector::copyOutRawData ( double **in )
  {
    d_Engine->copyOutRawData ( in );
  }

  inline
  void ManagedVector::setToScalar(double alpha)
  {
    d_Engine->setToScalar ( alpha );
    dataChanged();
    this->makeConsistent(CONSISTENT_SET);
  }

  inline
  double ManagedVector::L1Norm(void) const
  {
    return d_Engine->L1Norm ();
  }

  inline
  double ManagedVector::L2Norm(void) const
  {
    return d_Engine->L2Norm ();
  }

  inline
  double ManagedVector::maxNorm(void) const
  {
    return d_Engine->maxNorm ();
  }

  inline
  double ManagedVector::dot(const VectorOperations &x) const
  {
    if ( x.isA<ManagedVector> () )
      return d_Engine->dot ( *x.castTo<ManagedVector>().d_Engine );
    return Vector::dot (*this);
  }


  inline
  void ManagedVector::scale(double alpha, const VectorOperations &x)
  {
    if ( x.isA<ManagedVector> () )
    {
      d_Engine->scale ( alpha , *x.castTo<ManagedVector>().d_Engine );
    }
    else
    {
      Vector::scale ( alpha , x );
    }
    dataChanged();
  }

  inline
  void ManagedVector::scale(double alpha)
  {
    d_Engine->scale ( alpha );
    dataChanged();
  }

  inline
  void ManagedVector::add(const VectorOperations &x, const VectorOperations &y)
  {
    if ( x.isA<ManagedVector> () && y.isA<ManagedVector>() )
    {
      d_Engine->add ( *x.castTo<ManagedVector>().d_Engine , *y.castTo<ManagedVector>().d_Engine );
    }
    else
    {
      Vector::add ( x , y );
    }
    dataChanged();
  }

  inline
  void ManagedVector::subtract(const VectorOperations &x, const VectorOperations &y)
  {
    if ( x.isA<ManagedVector> () && y.isA<ManagedVector>() )
    {
      d_Engine->subtract ( *x.castTo<ManagedVector>().d_Engine , *y.castTo<ManagedVector>().d_Engine );
    }
    else
    {
      Vector::subtract ( x , y );
    }
    dataChanged();
  }

  inline
  void ManagedVector::multiply(const VectorOperations &x, const VectorOperations &y)
  {
    if ( x.isA<ManagedVector> () && y.isA<ManagedVector>() )
    {
      d_Engine->multiply ( *x.castTo<ManagedVector>().d_Engine , *y.castTo<ManagedVector>().d_Engine );
    }
    else
    {
      Vector::multiply ( x , y );
    }
    dataChanged();
  }

  inline
  void ManagedVector::divide(const VectorOperations &x, const VectorOperations &y)
  {
    if ( x.isA<ManagedVector> () && y.isA<ManagedVector>() )
    {
      d_Engine->divide ( *x.castTo<ManagedVector>().d_Engine , *y.castTo<ManagedVector>().d_Engine );
    }
    else
    {
      Vector::divide ( x , y );
    }
    dataChanged();
  }

  inline
  void ManagedVector::reciprocal(const VectorOperations &x)
  {
    d_Engine->reciprocal ( *x.castTo<ManagedVector>().d_Engine );
    dataChanged();
  }

  inline
  void ManagedVector::linearSum(double alpha, const VectorOperations &x,
          double beta, const VectorOperations &y)
  {
    if ( x.isA<ManagedVector> () && y.isA<ManagedVector>() )
    {
      d_Engine->linearSum ( alpha , *x.castTo<ManagedVector>().d_Engine , beta , *y.castTo<ManagedVector>().d_Engine );
    }
    else
    {
      Vector::linearSum ( alpha , x , beta , y );
    }
    dataChanged();
  }

  inline
  void ManagedVector::axpy(double alpha, const VectorOperations &x, const VectorOperations &y)
  {
    if ( x.isA<ManagedVector> () && y.isA<ManagedVector>() )
    {
      d_Engine->axpy ( alpha , *x.castTo<ManagedVector>().d_Engine , *y.castTo<ManagedVector>().d_Engine );
    }
    else
    {
      Vector::axpy ( alpha , x , y );
    }
    dataChanged();
  }

  inline
  void ManagedVector::axpby(double alpha, double beta, const VectorOperations &x)
  {
    if ( x.isA<ManagedVector> () )
    {
      d_Engine->axpby ( alpha , beta , *x.castTo<ManagedVector>().d_Engine );
    }
    else
    {
      Vector::axpby ( alpha , beta , x );
    }
    dataChanged();
  }

  inline
  void ManagedVector::abs(const VectorOperations &x)
  {
    if ( x.isA<ManagedVector> () )
    {
      d_Engine->abs ( *x.castTo<ManagedVector>().d_Engine );
    }
    else
    {
      Vector::abs ( x );
    }
    dataChanged();
  }

  inline
  double ManagedVector::min(void) const
  {
    return d_Engine->min ();
  }

  inline
  double ManagedVector::max(void) const
  {
    return d_Engine->max ();
  }

  inline
  void ManagedVector::setRandomValues(void)
  {
    d_Engine->setRandomValues ();
    dataChanged();
    this->makeConsistent(CONSISTENT_SET);
  }

  inline
  boost::shared_ptr<Vector>  ManagedVector::cloneVector ( const Variable::shared_ptr name ) const
  {
    boost::shared_ptr<Vector> retVal ( getNewRawPtr() );
    if ( !d_vBuffer )
    {
      retVal->castTo<ManagedVector>().d_Engine = d_Engine->cloneEngine ( VectorEngine::BufferPtr() );
    }
    retVal->setVariable ( name );
    return retVal;
  }


  inline
  ManagedVector::~ManagedVector ()
  {
  }

}
}

