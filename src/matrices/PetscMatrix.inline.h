
namespace AMP {
namespace LinearAlgebra {

  inline PetscMatrix::PetscMatrix() {
  }

  inline
  PetscMatrix::PetscMatrix( ParametersPtr params ) 
    : Matrix ( params ) 
  {
  }

  inline
  PetscMatrix::~PetscMatrix() 
  {
  }

  inline
  Mat &PetscMatrix::getMat()       
  { 
    return d_Mat; 
  }

  inline
  Mat  PetscMatrix::getMat() const 
  { 
    return d_Mat; 
  }

}
}

