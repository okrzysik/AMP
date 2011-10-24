namespace AMP {
namespace LinearAlgebra {

  inline
  EpetraVectorEngineParameters::EpetraVectorEngineParameters ( ManagedDataMap &m , AMP_MPI c )
       : VectorEngineParameters ( )
       , ManagedDataMap ( m )
       , d_comm ( c.getCommunicator() ) 
  {
  }

  inline
  EpetraVectorEngineParameters::EpetraVectorEngineParameters ( int local_size , int global_size , AMP_MPI c )
       : VectorEngineParameters ()
       , ManagedDataMap ( local_size , global_size ) ,
         d_comm ( c.getCommunicator() ) 
  {
  }

  inline
  EpetraVectorEngineParameters::EpetraVectorEngineParameters ( int local_size , int global_size , boost::shared_ptr<Epetra_Map> emap , AMP_MPI ecomm )
       : VectorEngineParameters ()
       , ManagedDataMap ( local_size , global_size )
       , d_emap ( emap )
       , d_comm ( ecomm )
  {
  }

  inline
  AMP_MPI  EpetraVectorEngineParameters::getEpetraComm () 
  { 
    return d_comm; 
  }

  inline
  EpetraVectorEngine::~EpetraVectorEngine () 
  { 
    ADD_COUNT ( "FLOPS" , d_epetraVector.Flops() ); 
  }

  inline
  Epetra_Vector  &EpetraVectorEngine::getEpetra_Vector()       
  { 
    return d_epetraVector; 
  }

  inline
  const Epetra_Vector  &EpetraVectorEngine::getEpetra_Vector() const 
  { 
    return d_epetraVector; 
  }

  inline
  AMP_MPI  EpetraVectorEngine::getComm () const 
  { 
    return getEngineParameters()->castTo<EpetraVectorEngineParameters>().getEpetraComm(); 
  }

  inline
  void  *EpetraVectorEngine::getDataBlock ( size_t i )
  {
    if ( i > 1 ) return 0;
    double *p;
    getEpetra_Vector().ExtractView ( &p );
    return p;
  }

  inline
  const void  *EpetraVectorEngine::getDataBlock ( size_t i ) const
  {
    if ( i > 1 ) return 0;
    double *p;
    getEpetra_Vector().ExtractView ( &p );
    return p;
  }

  inline
  size_t   EpetraVectorEngine::getLocalSize() const 
  { 
    return d_Params->castTo<EpetraVectorEngineParameters>().getLocalSize(); 
  }

  inline
  size_t   EpetraVectorEngine::getGlobalSize() const 
  { 
    return d_Params->castTo<EpetraVectorEngineParameters>().getGlobalSize(); 
  }

  inline
  bool EpetraVectorEngine::sameEngine ( VectorEngine &e ) const 
  { 
    return e.isA<EpetraVectorEngine>(); 
  }

  inline
  size_t EpetraVectorEngine::numberOfDataBlocks () const 
  { 
    return 1; 
  }

  inline
  size_t EpetraVectorEngine::sizeOfDataBlock ( size_t i ) const
  {
    if ( i != 0 )
      return 0;
    return getLocalSize();
  }

  inline
  void EpetraVectorEngine::putRawData ( double *in )
  {
    double *p;
    getEpetra_Vector().ExtractView ( &p );
    std::copy ( in , in+getLocalSize() , p );
  }

  inline
  void EpetraVectorEngine::copyOutRawData ( double **out )
  {
    *out = new double [ getLocalSize() ];
    getEpetra_Vector().ExtractCopy ( *out );
  }

}
}

