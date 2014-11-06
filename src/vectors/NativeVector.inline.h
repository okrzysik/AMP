

namespace AMP {
namespace LinearAlgebra {

  inline
  NativeVector::NativeVector ()
  {
  }

  inline
  NativeVector::NativeVector ( parameters_ptr params ) 
        : Vector ( AMP::dynamic_pointer_cast<VectorParameters> (params) )
  {
  }


}
}

