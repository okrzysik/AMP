

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

  inline
  size_t  NativeVector::numberOfDataBlocks () const 
  { 
    return 1; 
  }

  inline
  size_t  NativeVector::sizeOfDataBlock ( size_t i ) const 
  { 
    return i?0:getLocalSize(); 
  }

}
}

