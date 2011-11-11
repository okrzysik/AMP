#include "utils/Counter.h"

namespace AMP {
namespace LinearAlgebra {

  inline
  size_t Vector::getGlobalMaxID() const
  {
    return getGlobalSize();
  }

  inline
  size_t Vector::getLocalMaxID() const
  {
    return getLocalSize();
  }

  inline
  size_t Vector::getLocalStartID() const
  {
    return getCommunicationList()->getStartGID();
  }

  inline
  void  Vector::dataChanged()
  {
  }

  inline
  boost::shared_ptr<ParameterBase> Vector::getParameters ()
  {
    return  boost::shared_ptr<ParameterBase> ( );
  }

  inline
  void Vector::setDefaultRNG ( RNG::shared_ptr p )
  {
    d_DefaultRNG = p;
  }

  inline
  RNG::shared_ptr Vector::getDefaultRNG ()
  {
    if ( !d_DefaultRNG )
    {
      AMP_MPI globalComm(AMP_COMM_WORLD);
      int rank = globalComm.getRank();
      RNGParameters::shared_ptr params ( new RNGParameters ( RNGParameters::USE_GLOBAL_SEED ,
                                                             static_cast<size_t> ( rank ) ) );
      d_DefaultRNG = RNG::shared_ptr ( new RNG ( params ) );
    }
    return d_DefaultRNG;
  }

  inline
  void Vector::copyGhostValues ( const Vector &rhs )
  {
    std::copy ( rhs.d_Ghosts->begin() , rhs.d_Ghosts->end() , d_Ghosts->begin() );
  }


  inline
  bool Vector::containsGlobalElement ( size_t i )
  {
    if ( ( i >= d_CommList->getStartGID() ) && ( i < d_CommList->getStartGID() + d_CommList->numLocalRows() ) )
      return true;
    return std::find ( d_CommList->getGhostIDList().begin() ,
                       d_CommList->getGhostIDList().end() ,
                       (unsigned int) i ) != d_CommList->getGhostIDList().end();
  }


  inline
  void Vector::requireSameSize ( Vector &rhs )
  {
    if ( rhs.getLocalSize() != getLocalSize() )
    {
      AMP_ERROR( "Vectors are not of compatible size" );
    }
    if ( rhs.getGlobalSize() != getGlobalSize() )
    {
      AMP_ERROR( "Vectors are not of compatible size" );
    }
  }

  inline
  ConstVectorDataIterator  Vector::end () const
  {
    return ConstVectorDataIterator ( this , numberOfDataBlocks() , 0 );
  }

  inline
  ConstVectorDataIterator  Vector::begin () const
  {
    return ConstVectorDataIterator ( this , 0 , 0 );
  }

  inline
  VectorDataIterator  Vector::begin ()
  {
    return VectorDataIterator ( this , 0 , 0 );
  }

  inline
  VectorDataIterator  Vector::end ()
  {
    return VectorDataIterator ( this , numberOfDataBlocks() , 0 );
  }

  inline
  CommunicationList::shared_ptr  Vector::getCommunicationList () const
  {
    return d_CommList;
  }

  inline
  AMP::Discretization::DOFManager::shared_ptr  Vector::getDOFManager () const
  {
    return d_DOFManager;
  }

  inline
  size_t Vector::getGhostSize() const
  {
    return d_Ghosts->size();
  }

  inline
  void Vector::copyOutRawData ( double **in )
  {
    *in = new double [ getLocalSize() ];
    std::copy ( *in , *in + getLocalSize() , begin() );
  }

  inline
  const Variable::shared_ptr Vector::getVariable() const
  {
    return d_pVariable;
  }

  inline
  Variable::shared_ptr Vector::getVariable()
  {
    return d_pVariable;  // Fix this!
  }

  inline
  Vector::shared_ptr Vector::cloneVector () const
  {
    return cloneVector ( getVariable() );
  }


  inline
  Vector::~Vector()
  {
  }

  inline
  Vector::Vector(
      const Vector&rhs ) : VectorOperations ()
                         , boost::enable_shared_from_this<Vector> (rhs)
  {
  }


  inline
  void Vector::operator=(
      const Vector& )
  {
  }

  inline
  void Vector::setVariable(const Variable::shared_ptr name)
  {
    d_pVariable = name;
  }


  inline
  bool  Vector::equals ( Vector::shared_ptr &rhs , double  tol )
  {
    return equals (*rhs,tol);
  }

  inline
  void  Vector::copyVector ( const shared_ptr &src_vec )
  {
    copyVector ( *src_vec );
  }

  inline
  void  Vector::swapVectors ( shared_ptr &other )
  {
    swapVectors ( *other );
  }

  inline
  void  Vector::aliasVector ( shared_ptr &other )
  {
    aliasVector ( *other );
  }

  inline
  void  Vector::setValueByLocalID(int i, const double val)
  {
    setValuesByLocalID ( 1 , &i , &val );
  }

  inline
  void  Vector::setLocalValueByGlobalID(int i, const double val)
  {
    setLocalValuesByGlobalID ( 1 , &i , &val );
  }

  inline
  void  Vector::setValueByGlobalID(int i, const double val)
  {
    setValuesByGlobalID ( 1 , &i , &val );
  }

  inline
  void  Vector::addValueByLocalID(int i, const double val)
  {
    addValuesByLocalID ( 1 , &i , &val );
  }

  inline
  void  Vector::addLocalValueByGlobalID(int i, const double val)
  {
    addLocalValuesByGlobalID ( 1 , &i , &val );
  }

  inline
  void  Vector::addValueByGlobalID(int i, const double val)
  {
    addValuesByGlobalID ( 1 , &i , &val );
  }

  inline
  double  Vector::getValueByGlobalID ( int i ) const
  {
    double ans;
    getValuesByGlobalID ( 1 , &i , &ans );
    return ans;
  }

  inline
  double Vector::getLocalValueByGlobalID ( int i ) const
  {
    double ans;
    getLocalValuesByGlobalID ( 1 , &i , &ans );
    return ans;
  }

  inline
  double Vector::getValueByLocalID ( int ndx ) const
  {
    double ans;
    getValuesByLocalID ( 1 , &ndx , &ans );
    return ans;
  }

  inline
  void Vector::zero()
  {
    setToScalar (0.0);
    for ( size_t i = 0 ; i != d_Ghosts->size() ; i++ )
    {
      (*d_Ghosts)[i] = 0.0;
    }
    (*getUpdateStatus()) = NOT_UPDATING;
  }

  inline
  void Vector::scale ( double alpha , const shared_ptr &x )
  {
    scale ( alpha , *x );
  }

  inline
  void Vector::addScalar ( const shared_ptr &x , double alpha )
  {
    Vector::shared_ptr  one_vec = cloneVector ();
    one_vec->setToScalar ( 1. );
    axpy ( alpha , one_vec , x );
  }

  inline
  void Vector::add ( const shared_ptr &x , const shared_ptr &y )
  {
    add ( *x , *y );
  }

  inline
  void Vector::subtract ( const shared_ptr &x , const shared_ptr &y )
  {
    subtract ( *x , *y );
  }

  inline
  void Vector::multiply ( const shared_ptr &x , const shared_ptr &y )
  {
    multiply ( *x , *y );
  }

  inline
  void Vector::divide ( const shared_ptr &x , const shared_ptr &y )
  {
    divide ( *x , *y );
  }

  inline
  void Vector::reciprocal ( const shared_ptr &x )
  {
    reciprocal ( *x );
  }

  inline
  double Vector::minQuotient(const shared_ptr &x, const shared_ptr &y)
  {
    return(minQuotient(*x, *y));
  }

  inline
  double Vector::wrmsNorm(const shared_ptr &x, const shared_ptr &y)
  {
    return(wrmsNorm(*x, *y));
  }

  inline
  void Vector::linearSum ( double alpha , const shared_ptr &x , double beta , const shared_ptr &y )
  {
    linearSum ( alpha , *x , beta , *y );
  }

  inline
  void Vector::axpy ( double alpha , const shared_ptr &x , const shared_ptr &y )
  {
    axpy ( alpha , *x , *y );
  }

  inline
  void Vector::axpby(double alpha , double beta , const shared_ptr &x )
  {
    axpby ( alpha , beta , *x );
  }

  inline
  void Vector::abs ( const shared_ptr &x )
  {
    abs ( *x );
  }

  inline
  double Vector::dot ( const shared_ptr &x )
  {
    return dot ( *x );
  }

  inline
  AMP_MPI Vector::getComm() const
  {
    return d_CommList->getComm();
  }

  inline
  boost::shared_ptr<Vector::UpdateState>  Vector::getUpdateStatus () const
  {
    return d_UpdateState;
  }

  inline
  void Vector::setUpdateStatus ( boost::shared_ptr<UpdateState> rhs )
  {
    d_UpdateState = rhs;
  }

  inline
  void Vector::addCommunicationListToParameters ( CommunicationList::shared_ptr )
  {
  }

  inline
  void  Vector::aliasGhostBuffer ( shared_ptr in )
  {
    d_Ghosts = in->d_Ghosts;
  }

  inline
  std::ostream &operator << ( std::ostream &out , const Vector::shared_ptr p )
  {
    return operator << ( out , *p );
  }

}
}

