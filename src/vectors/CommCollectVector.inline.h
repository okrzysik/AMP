
namespace AMP {
namespace LinearAlgebra {

  inline
  void * CommCollectVector::getRawDataBlockAsVoid ( size_t i )
  {
    return getVectorEngine()->getDataBlock ( i );
  }

  inline
  const void * CommCollectVector::getRawDataBlockAsVoid ( size_t i ) const
  {
    const VectorEngine::shared_ptr t = getVectorEngine();
    return boost::const_pointer_cast<VectorEngine>(t)->getDataBlock ( i );
  }

  inline
  const VectorOperations & CommCollectVector::getOps ( const VectorOperations &ops ) const
  {
    if ( ops.isA<const CommCollectVector>() ) return *(ops.castTo<const CommCollectVector>().d_SmallCommVector);
    return ops;
  }

  inline
  VectorOperations & CommCollectVector::getOps ( VectorOperations &ops ) const
  {
    if ( ops.isA<CommCollectVector>() ) return *(ops.castTo<CommCollectVector>().d_SmallCommVector);
    return ops;
  }

  inline
  Vector::shared_ptr   CommCollectVector::internalView ( Vector::shared_ptr p , AMP_MPI lc , AMP_MPI rc )
  {
    CommCollectVector *r = new CommCollectVector;
    r->d_SmallCommVector = p->isA<CommCollectVector>()? p->castTo<CommCollectVector>().d_SmallCommVector : p;
    r->d_UpdateState = r->d_SmallCommVector->getUpdateStatus();
    if ( p->isA<DataChangeFirer>() )
    {
      p->castTo<DataChangeFirer>().registerListener ( r );
    }
    r->d_LargeComm = lc;
    r->d_RankComm = rc;
    r->setVariable ( p->getVariable() );
    return Vector::shared_ptr ( r );
  }

  inline
  Vector::shared_ptr   CommCollectVector::getSmallCommVector () 
  { 
    return d_SmallCommVector; 
  }

  inline
  const Vector::shared_ptr   CommCollectVector::getSmallCommVector () const 
  { 
    return d_SmallCommVector; 
  }

  inline
  Vector::shared_ptr  CommCollectVector::view ( Vector::shared_ptr p , AMP_MPI c )
  {
    CommCollectVector *r = new CommCollectVector;
    r->d_SmallCommVector = p->isA<CommCollectVector>()? p->castTo<CommCollectVector>().d_SmallCommVector : p;
    r->d_UpdateState = r->d_SmallCommVector->getUpdateStatus();
    if ( p->isA<DataChangeFirer>() )
    {
      p->castTo<DataChangeFirer>().registerListener ( r );
    }
    r->d_LargeComm = c;
    r->setVariable ( p->getVariable() );
    r->makeRankComm ();
    return Vector::shared_ptr ( r );
  }

  inline
  const Vector::shared_ptr  CommCollectVector::constView ( Vector::shared_ptr p , AMP_MPI c )
  {
    CommCollectVector *r = new CommCollectVector;
    r->d_SmallCommVector = p->isA<CommCollectVector>()? p->castTo<CommCollectVector>().d_SmallCommVector : p;
    r->d_UpdateState = r->d_SmallCommVector->getUpdateStatus();
    r->d_LargeComm = c;
    r->setVariable ( p->getVariable() );
    r->makeRankComm ();
    AMP_ERROR( "Cannot make constView of CommCollectVector.... yet" );
    return Vector::shared_ptr ( r );
  }


  inline
  Vector::iterator    CommCollectVector::begin() 
  { 
    return d_SmallCommVector->begin(); 
  }

  inline
  Vector::iterator    CommCollectVector::end() 
  { 
    return d_SmallCommVector->end(); 
  }

  inline
  Vector::const_iterator    CommCollectVector::begin() const 
  { 
    return d_SmallCommVector->begin(); 
  }

  inline
  Vector::const_iterator    CommCollectVector::end() const 
  { 
    return d_SmallCommVector->end(); 
  }

  inline
  void        CommCollectVector::selectInto ( const VectorSelector &a , Vector::shared_ptr b )
  {
    d_SmallCommVector->selectInto ( a , b );
  }

  inline
  void        CommCollectVector::dumpOwnedData ( std::ostream &out , size_t GIDoffset , size_t LIDoffset ) const
  {
    d_SmallCommVector->dumpOwnedData ( out , GIDoffset , LIDoffset );
  }

  inline
  void        CommCollectVector::dumpGhostedData ( std::ostream &out , size_t offset ) const
  {
    d_SmallCommVector->dumpGhostedData ( out , offset );
  }

  inline
  std::string CommCollectVector::type() const
  {
    return d_SmallCommVector->type() + " in a CommCollectVector";
  }

  inline
  AMP_MPI  CommCollectVector::getComm () const
  {
    return d_LargeComm;
  }

  inline
  size_t CommCollectVector::numberOfDataBlocks () const
  {
    return d_SmallCommVector->numberOfDataBlocks();
  }

  inline
  size_t CommCollectVector::sizeOfDataBlock ( size_t i ) const
  {
    return d_SmallCommVector->sizeOfDataBlock ( i );
  }

  inline
  Vector::shared_ptr  CommCollectVector::subsetVectorForVariable ( const Variable::shared_ptr &name )
  {
    Vector::shared_ptr  sub = d_SmallCommVector->subsetVectorForVariable ( name );
    /*
       inline
    if ( sub )
    {
      sub = internalView ( d_SmallCommVector->subsetVectorForVariable ( name ) , d_LargeComm , d_RankComm);
    }
    inline
    else
    {
      sub = Vector::shared_ptr ();
    }
    */
    return sub;
  }

  inline
  Vector::shared_ptr CommCollectVector::cloneVector(const Variable::shared_ptr name) const
  {
    return internalView ( d_SmallCommVector->cloneVector ( name ) , d_LargeComm , d_RankComm );
  }

  inline
  void CommCollectVector::swapVectors(Vector &other)
  {
    d_SmallCommVector->swapVectors ( other.castTo<CommCollectVector>().d_SmallCommVector );
  }

  inline
  void CommCollectVector::aliasVector(Vector &other)
  {
    d_SmallCommVector->aliasVector ( other.castTo<CommCollectVector>().d_SmallCommVector );
  }

  inline
  void CommCollectVector::setToScalar(double alpha)
  {
    d_SmallCommVector->setToScalar ( alpha );
  }

  inline
  void CommCollectVector::scale(double alpha, const VectorOperations &x)
  {
    d_SmallCommVector->scale ( alpha , x.castTo<CommCollectVector>().d_SmallCommVector );
  }

  inline
  void CommCollectVector::scale(double alpha)
  {
    d_SmallCommVector->scale ( alpha );
  }

  inline
  void CommCollectVector::add(const VectorOperations &x, const VectorOperations &y)
  {
    d_SmallCommVector->add ( getOps ( x ) , getOps ( y ) );
  }

  inline
  void CommCollectVector::subtract(const VectorOperations &x, const VectorOperations &y)
  {
    d_SmallCommVector->subtract ( getOps ( x ) , getOps ( y ) );
  }

  inline
  void CommCollectVector::multiply( const VectorOperations &x, const VectorOperations &y)
  {
    d_SmallCommVector->multiply ( getOps ( x ) , getOps ( y ) );
  }

  inline
  void CommCollectVector::divide( const VectorOperations &x, const VectorOperations &y)
  {
    d_SmallCommVector->divide ( getOps ( x ) , getOps ( y ) );
  }

  inline
  void CommCollectVector::reciprocal(const VectorOperations &x)
  {
    d_SmallCommVector->reciprocal ( getOps ( x ) );
  }

  inline
  void CommCollectVector::linearSum(double alpha, const VectorOperations &x,
          double beta, const VectorOperations &y)
  {
    d_SmallCommVector->linearSum ( alpha , getOps ( x ) , beta , getOps ( y ));
  }

  inline
  void CommCollectVector::axpy(double alpha, const VectorOperations &x, const VectorOperations &y)
  {
    d_SmallCommVector->axpy ( alpha , getOps ( x ) , getOps ( y ) );
  }

  inline
  void CommCollectVector::axpby(double alpha, double beta, const VectorOperations &x)
  {
    d_SmallCommVector->axpby ( alpha , beta , getOps ( x ) );
  }

  inline
  void CommCollectVector::abs(const VectorOperations &x)
  {
    d_SmallCommVector->abs ( getOps ( x ) );
  }

  inline
  void CommCollectVector::setRandomValues(void)
  {
    d_SmallCommVector->setRandomValues ();
  }

  inline
  void CommCollectVector::setValuesByLocalID ( int num , int *indices , const double *vals )
  {
    d_SmallCommVector->setValuesByLocalID ( num , indices , vals );
  }

  inline
  void CommCollectVector::setLocalValuesByGlobalID ( int num , int *indices , const double *vals )
  {
    d_SmallCommVector->setLocalValuesByGlobalID ( num , indices , vals );
  }

  inline
  void CommCollectVector::setValuesByGlobalID ( int num , int *indices , const double *vals )
  {
    d_SmallCommVector->setValuesByGlobalID ( num , indices , vals );
  }

  inline
  void CommCollectVector::addValuesByLocalID ( int num , int *indices , const double *vals )
  {
    d_SmallCommVector->addValuesByLocalID ( num , indices , vals );
  }

  inline
  void CommCollectVector::addLocalValuesByGlobalID ( int num , int *indices , const double *vals )
  {
    d_SmallCommVector->addLocalValuesByGlobalID ( num , indices , vals );
  }

  inline
  void CommCollectVector::addValuesByGlobalID ( int num , int *indices , const double *vals )
  {
    d_SmallCommVector->addValuesByGlobalID ( num , indices , vals );
  }

  inline
  void CommCollectVector::getValuesByGlobalID ( int numVals , int *ndx , double *vals ) const
  {
    d_SmallCommVector->getValuesByGlobalID ( numVals , ndx , vals );
  }

  inline
  void CommCollectVector::getLocalValuesByGlobalID ( int numVals , int *ndx , double *vals ) const
  {
    d_SmallCommVector->getLocalValuesByGlobalID ( numVals , ndx , vals );
  }

  inline
  void CommCollectVector::getValuesByLocalID ( int numVals , int *ndx , double *vals ) const
  {
    d_SmallCommVector->getValuesByLocalID ( numVals , ndx , vals );
  }

  inline
  void CommCollectVector::makeConsistent ( ScatterType  t )
  {
    d_SmallCommVector->makeConsistent ( t );
  }

  inline
  void CommCollectVector::assemble()
  {
    d_SmallCommVector->assemble();
  }

  inline
  size_t CommCollectVector::getLocalSize() const
  {
    return d_SmallCommVector->getLocalSize();
  }


  inline
  size_t CommCollectVector::getGhostSize() const
  {
    return d_SmallCommVector->getGhostSize();
  }

  inline
  void   CommCollectVector::putRawData ( double *d )
  {
    d_SmallCommVector->putRawData ( d );
  }

  inline
  VectorEngine::BufferPtr  CommCollectVector::getNewBuffer()
  {
    return VectorEngine::BufferPtr ();
  }
 
  inline
  bool                     CommCollectVector::sameEngine ( VectorEngine &rhs ) const
  {
    return getVectorEngine()->sameEngine ( rhs );
  }

  inline
  VectorEngine::shared_ptr CommCollectVector::cloneEngine ( VectorEngine::BufferPtr p ) const
  {
    return boost::dynamic_pointer_cast<VectorEngine> ( cloneVector ( getVariable () ) );
  }

  inline
  void                     CommCollectVector::swapEngines ( VectorEngine::shared_ptr p )
  {
    getVectorEngine()->swapEngines ( p );
  }

  inline
  const void              *CommCollectVector::getDataBlock ( size_t i ) const
  {
    return getVectorEngine()->getDataBlock ( i );
  }

  inline
  void                    *CommCollectVector::getDataBlock ( size_t i )
  {
    return getVectorEngine()->getDataBlock ( i );
  }


  inline
  void   CommCollectVector::copyOutRawData ( double ** out )
  {
    d_SmallCommVector->copyOutRawData ( out );
  }

}
}

