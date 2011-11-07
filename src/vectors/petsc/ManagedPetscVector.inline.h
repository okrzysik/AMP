
namespace AMP {
namespace LinearAlgebra {

  inline
  ManagedVector *ManagedPetscVector::getNewRawPtr () const
  { 
    return new ManagedPetscVector ( boost::dynamic_pointer_cast<VectorParameters> ( d_pParameters ) ); 
  }

  inline
  bool  ManagedPetscVector::constructedWithPetscDuplicate ()
  { 
    return d_bMadeWithPetscDuplicate; 
  }

  inline
  ManagedPetscVector  *ManagedPetscVector::rawClone () const
  {
    parameters_ptr p ( new parameters );
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
    ManagedPetscVector *retVal = new ManagedPetscVector ( boost::dynamic_pointer_cast<VectorParameters> ( p ) );
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
