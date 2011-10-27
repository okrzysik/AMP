#ifndef  included_AMP_FineSearch_H
#define  included_AMP_FineSearch_H

#include <boost/shared_ptr.hpp>
#include <map>

namespace AMP {
namespace Mesh {

  template <typename MESH , typename NORMAL_MAP , typename ITERATOR1 , typename ITERATOR2 = ITERATOR1>
  class  FineSearchParameters
  {
    public:
      typedef NORMAL_MAP                                                           NormalMap;
      typedef typename NormalMap::mapped_type                                      Vector;
      typedef boost::shared_ptr<FineSearchParameters<MESH,NORMAL_MAP,ITERATOR1,ITERATOR2> >   shared_ptr;

      typename MESH::shared_ptr        d_mesh1, 
                                       d_mesh2;


      NORMAL_MAP                       d_Normals;
      std::multimap<size_t , Vector>   d_Intersections;
  }
  ;

  template <typename MESH , typename NORMAL_MAP , typename ITERATOR1 , typename ITERATOR2 = ITERATOR1>
  class  FineSearch
  {
    public:
      typedef boost::shared_ptr < CoarseSearch <MESH,ITERATOR1,ITERATOR2> >   shared_ptr;
      typedef FineSearchParameters<MESH,NORMAL_MAP,ITERATOR1,ITERATOR2>       Parameters;

      typedef std::multimap<ITERATOR1 , ITERATOR2>        FirstToSecondMap;
      typedef std::multimap<ITERATOR2 , ITERATOR1>        SecondToFirstMap;

      typedef NORMAL_MAP                                  NormalMap;
      typedef typename NormalMap::mapped_type             Normal;

      typedef Normal                                      Vector;
      typedef std::multimap<size_t , Vector>              IntersectionMultiMap;

      typedef typename MESH::Node                         Node;
      typedef typename MESH::shared_ptr                   MeshSharedPtr;

    protected:
      class J
      {
        public:
          double mat[9];

          void  GEPP ( Vector &x );
      };

      typename Parameters::shared_ptr    d_Parameters;

      bool  isObject1InContactWithObject2 ( MeshNode &n , MeshElement &e , MeshSharedPtr );
      bool  isObject1InContactWithObject2 ( MeshElement &n , MeshNode &e , MeshSharedPtr );
      bool  isObject1InContactWithObject2 ( MeshElement &e1 , MeshElement &e2 , MeshSharedPtr );

      void  interpolatePosition ( MeshNode &n , MeshElement &e , Normal &norm , MeshSharedPtr , Vector &input , Vector &answer , J &deriv , Vector &interp );
      void  findIntersection ( MeshNode &n , MeshElement &e , Normal &norm , MeshSharedPtr , Vector & , Vector & );


      template <typename CONTAINER>
      void  refineMap ( CONTAINER & , CONTAINER & , typename MESH::shared_ptr );

      template <typename IN_CONTAINER , typename OUT_CONTAINER>
      void  buildReverseMap ( IN_CONTAINER & , OUT_CONTAINER & );

    public:
      FineSearch ( typename Parameters::shared_ptr p );

      virtual void  refineNeighbors ( FirstToSecondMap &inM1ToM2 ,
                                      SecondToFirstMap &inM2ToM1 ,
                                      FirstToSecondMap &outM1ToM2 ,
                                      SecondToFirstMap &outM2ToM1 );

  };

}
}

#include "FineSearch.tmpl.h"
#endif
