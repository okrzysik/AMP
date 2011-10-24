namespace AMP {
namespace LinearAlgebra {

  inline
  PetscVector::PetscVector ()
  {
    d_PetscRandom = 0;
  }

  inline
  PetscRandom &PetscVector::getPetscRandom ( AMP_MPI comm )
  {
    if ( d_PetscRandom == 0 )
    {
      d_PetscRandom = new PetscRandom;
      PetscRandomCreate ( comm.getCommunicator() , d_PetscRandom );
      PetscRandomSetType ( *d_PetscRandom , PETSCRAND48 );  // This is a horrible RNG for
                                                            // stochastic simulation.  Do not 
                                                            // use.
    }
    return *d_PetscRandom;
  }

  inline
  Vec &PetscVector::getVec() 
  { 
    return d_petscVec; 
  }

  inline
  Vec  PetscVector::getVec() const 
  { 
    return d_petscVec; 
  }

  inline
  Vector::shared_ptr  PetscVector::createView ( Vector::shared_ptr AmpVector ) 
  { 
    DEPRECATED("createView","view");
    return view ( AmpVector ); 
  }

  inline
  const Vector::shared_ptr  PetscVector::createConstView ( Vector::shared_ptr AmpVector ) 
  { 
    DEPRECATED("createConstView","constView");
    return constView ( AmpVector ); 
  }

  inline
  PetscVector::~PetscVector()
  {
    if ( d_PetscRandom )
    {
      PetscRandomDestroy ( *d_PetscRandom );
      delete d_PetscRandom;
    }
  }

}
}

