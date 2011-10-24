#ifndef  included_AMP_ContactTypes_h
#define  included_AMP_ContactTypes_h

#include "ComputeNormals.h"

namespace AMP {
namespace Mesh {

  namespace Contact
  {
    template <typename MESH>
    class SerialPenaltyMethod
    {
      public:
        typedef typename MESH::BoundaryNodeIterator  Iterator1;
        typedef typename MESH::BoundarySideIterator  Iterator2;

        typedef typename Iterator1::value_type       Type1;
        typedef typename Iterator2::value_type       Type2;

        typedef Iterator2                            FaceIterator;

        enum { NormalType = ComputeNormals<MESH,FaceIterator>::Node };
    }
    ;

    template <typename MESH>
    class SerialMortarMethod
    {
      public:
        typedef typename MESH::BoundarySideIterator  Iterator1;
        typedef typename MESH::BoundarySideIterator  Iterator2;

        typedef typename Iterator1::value_type       Type1;
        typedef typename Iterator2::value_type       Type2;

        typedef Iterator2                            FaceIterator;

        enum { NormalType = ComputeNormals<MESH,FaceIterator>::Facet };
    }
    ;

  }

}
}

#endif
