

namespace AMP {
namespace LinearAlgebra {

    // this will generate an error. This is already taken care of
    // by calling the getRawDataBlockAsVoid function. The latter 
    // can be virtual
#if 0
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
#endif

}
}

