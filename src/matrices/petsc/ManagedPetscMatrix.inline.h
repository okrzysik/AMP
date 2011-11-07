

namespace AMP {
namespace LinearAlgebra {

  inline
  ManagedPetscMatrix::~ManagedPetscMatrix() 
  { 
    MatDestroy ( d_Mat ); 
  }

  inline
  Matrix::shared_ptr  ManagedPetscMatrix::cloneMatrix () const 
  { 
    return shared_ptr ( new ManagedPetscMatrix ( *this ) ); 
  }

}
}

