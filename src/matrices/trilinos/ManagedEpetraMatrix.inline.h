namespace AMP {
namespace LinearAlgebra {

  inline
  void ManagedEpetraMatrix::fillComplete ()
  {
    EpetraMatrix::fillComplete();
  }

  inline
  const int * ManagedEpetraMatrixParameters::entryList () const       
  { 
    return &*d_vEntriesPerRow.begin(); 
  }

  inline
  int * ManagedEpetraMatrixParameters::entryList()              
  { 
    return &(d_vEntriesPerRow[0]); 
  }

  inline
  void   ManagedEpetraMatrixParameters::setEntriesInRow ( int row , int entries ) 
  { 
    d_vEntriesPerRow[row] = entries; 
  }

  inline
  int  & ManagedEpetraMatrixParameters::entriesInRow ( int i )       
  { 
    return d_vEntriesPerRow[i]; 
  }

  inline
  int    ManagedEpetraMatrixParameters::entriesInRow ( int i ) const 
  { 
    return d_vEntriesPerRow[i]; 
  }

  inline
  int    ManagedEpetraMatrixParameters::maxEntitiesInRow () const 
  { 
    return *std::max_element ( d_vEntriesPerRow.begin() , d_vEntriesPerRow.end() ); 
  }

  inline
  bool  ManagedEpetraMatrixParameters::isSquare () 
  { 
    return d_ColBase < 0; 
  }

  inline
  boost::shared_ptr < Epetra_Map >   ManagedEpetraMatrixParameters::getEpetraRowMapPtr () 
  { 
    return d_eRowMap; 
  }

  inline
  boost::shared_ptr < Epetra_Map >   ManagedEpetraMatrixParameters::getEpetraColMapPtr () 
  { 
    return d_eColMap; 
  }

  inline
  AMP_MPI   ManagedEpetraMatrixParameters::getEpetraComm () 
  { 
    return d_comm; 
  }

  inline
  ManagedEpetraMatrix::ManagedEpetraMatrix ( Epetra_CrsMatrix *m , bool dele ) 
     : EpetraMatrix ( m , dele ) 
     , ManagedMatrix ( MatrixParameters::shared_ptr() ) 
  {
  }

  inline
  double  ManagedEpetraMatrix::L1Norm() const
  { 
    return d_epetraMatrix->NormOne(); 
  }

  inline
  Matrix::shared_ptr ManagedEpetraMatrix::cloneMatrix () const
  {
    ManagedEpetraMatrix *r = new ManagedEpetraMatrix ( d_epetraMatrix );
    r->d_DeleteMatrix = true;
    return Matrix::shared_ptr ( r );
  }

}
}

