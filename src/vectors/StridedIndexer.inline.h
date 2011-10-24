

namespace AMP {
namespace LinearAlgebra {

  inline
  StridedIndexer::StridedIndexer ( size_t  offset ,   size_t  stride )
                     : d_Offset ( offset )
                     , d_Stride ( stride ) 
  {
  }

  inline
  bool     StridedIndexer::isInSub ( size_t input ) const
  {
    return (input >= d_Offset) && ((( input - d_Offset ) % d_Stride) == 0);
  }

  inline
  size_t   StridedIndexer::getSubID ( size_t input ) const
  {
    return ( input - d_Offset ) / d_Stride;
  }

  inline
  size_t   StridedIndexer::getSuperID ( size_t input ) const
  {
    return d_Stride * input + d_Offset;
  }

  inline
  size_t   StridedIndexer::getIncrement ( size_t ) const
  {
    return d_Stride;
  }

  inline
  size_t   StridedIndexer::getNumLocalElements ( Vector::shared_ptr p ) const
  {
    return p->getLocalSize() / d_Stride;
  }

}
}

