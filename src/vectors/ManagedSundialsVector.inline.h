
#ifdef USE_SUNDIALS

namespace AMP {
namespace LinearAlgebra {

  inline
  ManagedVector *ManagedSundialsVector::getNewRawPtr () const
  { 
    return new ManagedSundialsVector ( boost::dynamic_pointer_cast<VectorParameters> ( d_pParameters ) ); 
  }

  inline
  ManagedSundialsVector  *ManagedSundialsVector::rawClone () const
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
    ManagedSundialsVector *retVal = new ManagedSundialsVector ( boost::dynamic_pointer_cast<VectorParameters> ( p ) );
    return retVal;
  }

  inline
  std::string ManagedSundialsVector::type() const
  {
    std::string retVal = "Managed SUNDIALS Vector";
    retVal += ManagedVector::type();
    return retVal;
  }

  inline
  Vector::shared_ptr  ManagedSundialsVector::cloneVector ( const Variable::shared_ptr var ) const
  {
      Vector::shared_ptr  retVal ( rawClone() );
      retVal->setVariable ( var );
      return retVal;
  }

  inline
  void ManagedSundialsVector::assemble () 
  {
  }

}
}

#endif
