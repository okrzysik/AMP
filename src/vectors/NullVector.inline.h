
namespace AMP {
namespace LinearAlgebra {

  inline
  Vector::shared_ptr   NullVector::create ( const std::string &name )
  {
    return Vector::shared_ptr  ( new NullVector ( Variable::shared_ptr ( new Variable ( name ) ) ) );
  }

  inline
  Vector::shared_ptr  NullVector::create ( const Variable::shared_ptr var )
  {
//    if ( !var->isA<NullVariable>() )
//    {
//      AMP_ERROR( "Invalid variable to create a NullVector" );
//    }
    return Vector::shared_ptr ( new NullVector ( var ) );
  }

  /// \cond UNDOCUMENTED
  inline
  NullVector::NullVector( Variable::shared_ptr var )
  {
    setVariable ( var );
  }


  inline
  NullVector::~NullVector() 
  {
  }

  inline
  boost::shared_ptr<ParameterBase> NullVector::getParameters () 
  { 
    return boost::shared_ptr<ParameterBase> (); 
  }

  inline
  Vector::shared_ptr NullVector::cloneVector(const Variable::shared_ptr name) const 
  { 
    return create ( name ); 
  }

  inline
  void NullVector::copyVector(const Vector &) 
  {
  }

  inline
  void NullVector::swapVectors(Vector &) 
  {
  }

  inline
  void NullVector::aliasVector(Vector & ) 
  {
  }


  inline
  void NullVector::setToScalar(double ) 
  {
  }

  inline
  void NullVector::scale(double , const VectorOperations &) 
  {
  }

  inline
  void NullVector::scale(double ) 
  {
  }

  inline
  void NullVector::addScalar(const VectorOperations &, double ) 
  {
  }

  inline
  void NullVector::add(const VectorOperations &, const VectorOperations &) 
  {
  }

  inline
  void NullVector::subtract(const VectorOperations &, const VectorOperations &) 
  {
  }

  inline
  void NullVector::multiply( const VectorOperations &, const VectorOperations &) 
  {
  }

  inline
  void NullVector::divide( const VectorOperations &, const VectorOperations &) 
  {
  }

  inline
  void NullVector::reciprocal(const VectorOperations &) 
  {
  }

  inline
  void NullVector::linearSum(double , const VectorOperations &,
          double , const VectorOperations &) 
  {
  }

  inline
  void NullVector::axpy(double , const VectorOperations &, const VectorOperations &) 
  {
  }

  inline
  void NullVector::axpby(double , double, const VectorOperations &) 
  {
  }

  inline
  void NullVector::abs(const VectorOperations &) 
  {
  }

  inline
  double NullVector::min(void) const 
  { 
    return 0.0; 
  }

  inline
  double NullVector::max(void) const 
  { 
    return 0.0; 
  }

  inline
  void NullVector::setRandomValues(void) 
  {
  }

  inline
  void NullVector::setValuesByLocalID ( int , size_t * , const double * ) 
  { 
    AMP_ERROR( "Can't set values for NullVector" ); 
  }

  inline
  void NullVector::setLocalValuesByGlobalID ( int , size_t * , const double * ) 
  { 
    AMP_ERROR( "Can't set values for NullVector" ); 
  }

  inline
  void NullVector::addValuesByLocalID ( int , size_t * , const double * ) 
  { 
    AMP_ERROR( "Can't set values for NullVector" ); 
  }

  inline
  void NullVector::addLocalValuesByGlobalID ( int , size_t * , const double * ) 
  { 
    AMP_ERROR( "Can't set values for NullVector" ); 
  }

  inline
  void NullVector::getLocalValuesByGlobalID ( int , size_t * , double * ) const 
  { 
    AMP_ERROR( "Can't set values for NullVector" ); 
  }


  inline
  void NullVector::makeConsistent ( ScatterType  ) 
  {
  }

  inline
  void NullVector::assemble() 
  {
  }

  inline
  double NullVector::L1Norm(void) const 
  { 
    return 0.0; 
  }

  inline
  double NullVector::L2Norm(void) const 
  { 
    return 0.0; 
  }

  inline
  double NullVector::maxNorm(void) const 
  { 
    return 0.0; 
  }

  inline
  double NullVector::dot(const VectorOperations &) const 
  { 
    return 0.0; 
  }

  inline
  void NullVector::putRawData ( double * ) 
  {
  }

  inline
  size_t NullVector::getLocalSize() const 
  { 
    return 0; 
  }

  inline
  size_t NullVector::getGlobalSize() const 
  { 
    return 0; 
  }

  inline
  size_t NullVector::getGhostSize() const 
  { 
    return 0; 
  }

  inline
  void NullVector::setCommunicationList ( CommunicationList::shared_ptr  )
  { 
  }

  inline
  size_t NullVector::numberOfDataBlocks () const 
  { 
    return 0; 
  }

  inline
  size_t NullVector::sizeOfDataBlock ( size_t ) const 
  { 
    return 0; 
  }

  inline
  void *NullVector::getRawDataBlockAsVoid ( size_t ) 
  { 
    return 0; 
  }

  inline
  const void *NullVector::getRawDataBlockAsVoid ( size_t ) const 
  { 
    return 0; 
  }

  inline
  void NullVector::addCommunicationListToParameters ( CommunicationList::shared_ptr ) 
  {
  }

      /// \endcond

}
}

