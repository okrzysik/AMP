#ifndef  included_AMP_ContactManagerTmpl_h
#define  included_AMP_ContactManagerTmpl_h

#include <iostream>
#include <set>

namespace AMP {
namespace Mesh {

  template <typename MESH , 
            template <typename M> class CONTACT_ALGO , 
            template <typename M1 , typename I1> class NORMAL_FCN ,
            template <typename M2 , typename I2 , typename I4> class COARSE_SEARCH , 
            template <typename M3 , typename MAP , typename I3 , typename I5> class FINE_SEARCH>
  class ContactManagerTmpl
  {
    public:
      typedef MESH                                              Adapter;
      typedef CONTACT_ALGO<MESH>                                ContactAlgorithm;

      typedef typename ContactAlgorithm::Type1                  Side1Type;
      typedef typename ContactAlgorithm::Type2                  Side2Type;

      typedef typename ContactAlgorithm::Iterator1              Side1Iterator;
      typedef typename ContactAlgorithm::Iterator2              Side2Iterator;

      typedef typename ContactAlgorithm::FaceIterator           Side1FaceIterator;

      typedef std::set<typename ContactAlgorithm::Type1>        Side1Set;
      typedef std::set<typename ContactAlgorithm::Type2>        Side2Set;

      typedef NORMAL_FCN<MESH,Side1FaceIterator>                NormalFunction;
      typedef typename NormalFunction::NormalMap                NormalMap;
      typedef typename NormalFunction::Normal                   Normal;

      typedef COARSE_SEARCH<MESH,Side1Iterator,Side2Iterator>             CoarseSearch;
      typedef FINE_SEARCH<MESH,NormalMap,Side1Iterator,Side2Iterator>     FineSearch;

      typedef typename CoarseSearch::FirstToSecondMap           Side1ToSide2MultiMap;
      typedef typename CoarseSearch::SecondToFirstMap           Side2ToSide1MultiMap;


      static void  computeContactSets ( typename MESH::shared_ptr  mesh1 ,
                                        typename MESH::shared_ptr  mesh2 ,
                                        Side1Iterator      side1Begin ,
                                        Side1Iterator      side1End ,
                                        Side1FaceIterator  side1FaceBegin ,
                                        Side1FaceIterator  side1FaceEnd ,
                                        Side2Iterator      side2Begin ,
                                        Side2Iterator      side2End ,
                                        Side1Set          &side1ActiveSet ,
                                        Side2Set          &side2ActiveSet ,
                                        Side1ToSide2MultiMap  &side1ToSide2Map ,
                                        Side2ToSide1MultiMap  &side2ToSide1Map );
  };

}
}

#include "ContactManagerTmpl.tmpl.h"

#include "ContactTypes.h"
#include "CoarseSearch.h"
#include "FineSearch.h"
#include "MeshType.h"
namespace AMP {
namespace Mesh {

  typedef ContactManagerTmpl<MESH_TYPE,Contact::SerialPenaltyMethod,ComputeNormals,CoarseSearch,FineSearch>
                 SerialPenaltyContactSearch;

  typedef ContactManagerTmpl<MESH_TYPE,Contact::SerialMortarMethod,ComputeNormals,CoarseSearch,FineSearch>
                 SerialMortarContactSearch;
}
}

#endif
