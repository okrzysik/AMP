

namespace AMP {
namespace LinearAlgebra {

     /// \cond

  template <typename RETURN_TYPE>
  RETURN_TYPE *  NullVector::getRawDataBlock ()
  {
    return static_cast <RETURN_TYPE *> (0);
  }

  template <typename RETURN_TYPE>
  const RETURN_TYPE *  NullVector::getRawDataBlock() const
  {
    return static_cast <const RETURN_TYPE *> (0);
  }


      /// \endcond

}
}

