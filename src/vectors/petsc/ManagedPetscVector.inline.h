
namespace AMP {
namespace LinearAlgebra {

  inline
  ManagedVector *ManagedPetscVector::getNewRawPtr () const
  { 
    return new ManagedPetscVector ( AMP::dynamic_pointer_cast<VectorParameters> ( d_pParameters ) ); 
  }

  inline
  bool  ManagedPetscVector::constructedWithPetscDuplicate ()
  { 
    return d_bMadeWithPetscDuplicate; 
  }

  inline
  ManagedPetscVector  *ManagedPetscVector::rawClone () const
  {
    AMP::shared_ptr<ManagedVectorParameters> p ( new ManagedPetscVectorParameters );
    if ( !d_vBuffer )
    {
      p->d_Engine = d_Engine->cloneEngine ( VectorEngine::BufferPtr () );
    }
    else
    {
      p->d_Buffer = VectorEngine::BufferPtr ( new VectorEngine::Buffer ( d_vBuffer->size() ) );
      p->d_Engine = d_Engine->cloneEngine ( p->d_Buffer );
    }
    p->d_CommList = getCommunicationList();
    p->d_DOFManager = getDOFManager();
    ManagedPetscVector *retVal = new ManagedPetscVector ( AMP::dynamic_pointer_cast<VectorParameters> ( p ) );
    return retVal;
  }

  inline
  Vector::shared_ptr  ManagedPetscVector::cloneVector ( const Variable::shared_ptr p ) const
  {
    Vector::shared_ptr  retVal ( rawClone() );
    retVal->setVariable ( p );
    return retVal;
  }

  inline
  std::string ManagedPetscVector::type() const
  {
    std::string retVal = "Managed PETSc Vector";
    retVal += ManagedVector::type();
    return retVal;
  }

  inline
  void  ManagedPetscVector::assemble() 
  {
  }

}
}
