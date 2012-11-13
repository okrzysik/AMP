

namespace AMP {
namespace LinearAlgebra {

  inline
  NativeVector::NativeVector ()
  {
  }

  inline
  NativeVector::NativeVector ( parameters_ptr params ) 
        : Vector ( boost::dynamic_pointer_cast<VectorParameters> (params) )
  {
  }


}
}

