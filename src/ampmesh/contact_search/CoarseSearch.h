#ifndef  included_AMP_CoarseSearch_H
#define  included_AMP_CoarseSearch_H

#include <map>
#include "utils/AMP_MPI.h"
#include <ampmesh/MeshNode.h>
#include <ampmesh/MeshElement.h>
#include <boost/shared_ptr.hpp>

namespace AMP {
namespace Mesh {

  template <typename MESH , typename ITERATOR1 , typename ITERATOR2 = ITERATOR1>
  class CoarseSearchParameters
  {
    public:
      typedef boost::shared_ptr < CoarseSearchParameters <MESH,ITERATOR1,ITERATOR2> >   shared_ptr;
      ITERATOR1                    d_beginM1 , d_endM1;
      ITERATOR2                    d_beginM2 , d_endM2;
      typename MESH::shared_ptr    d_mesh1   , d_mesh2;
      AMP_MPI                      d_comm;

      // Stupid libMesh doesn't allow unassigned iterators!!!
      CoarseSearchParameters ( ITERATOR1 a , ITERATOR1 b , ITERATOR2 c , ITERATOR2 d ) :
        d_beginM1 ( a ),
        d_endM1 ( b ),
        d_beginM2 ( c ),
        d_endM2 ( d )
      {
      }
  };

  template <typename MESH , typename ITERATOR1 , typename ITERATOR2 = ITERATOR1>
  class CoarseSearch
  {
    public:
      typedef boost::shared_ptr < CoarseSearch <MESH,ITERATOR1,ITERATOR2> >   shared_ptr;
      typedef CoarseSearchParameters<MESH,ITERATOR1,ITERATOR2>                Parameters;

      typedef std::multimap<ITERATOR1 , ITERATOR2>        FirstToSecondMap;
      typedef std::multimap<ITERATOR2 , ITERATOR1>        SecondToFirstMap;

    protected:
      struct Point
      {
        double x , y , z;
      };

      typedef std::pair<struct Point,struct Point>    Box;


      class Index
      {
        public:
          int  x , y , z;

          Index  getNeighbor ( size_t i ) const
          {
            Index answer;
            answer.x = x - 1 + (i%3);
            answer.y = y - 1 + (i/3)%3;
            answer.z = z - 1 + (i/9)%3;
            return answer;
          }

          bool operator < ( const Index &rhs ) const
          {
            if ( x > rhs.x ) return false;
            if ( x < rhs.x ) return true;
            if ( y > rhs.y ) return false;
            if ( y < rhs.y ) return true;
            if ( z < rhs.z ) return true;
            return false;
          }

          bool operator == ( const Index &rhs ) const
          {
            if ( x != rhs.x ) return false;
            if ( y != rhs.y ) return false;
            if ( z != rhs.z ) return false;
            return true;
          }

          bool operator != ( const Index &rhs ) const
          {
            return !operator == ( rhs );
          }
      };

      Box                                 d_BoundingBox;
      double                              d_SmallestPartitionSide;
      typename Parameters::shared_ptr     d_Parameters;

      std::multimap<Index , ITERATOR1>         d_Mesh1Decomposition;
      std::multimap<Index , ITERATOR2>         d_Mesh2Decomposition;

      template <typename ITERATOR>
      double   growBoundingBox ( ITERATOR begin , ITERATOR end , typename MESH::shared_ptr );


      Index    computeIndex ( MeshNode & );

      template <typename ITERATOR>
      void     putInContainer ( MeshNode & , ITERATOR , std::multimap<Index,ITERATOR> & , typename MESH::shared_ptr );

      template <typename ITERATOR>
      void     putInContainer ( MeshElement & , ITERATOR , std::multimap<Index,ITERATOR> & , typename MESH::shared_ptr );

      template <typename ITERATOR>
      void     sortObjects ( ITERATOR begin , ITERATOR end , typename MESH::shared_ptr , std::multimap<Index,ITERATOR>  &decomp );

      template <typename ITERATOR , typename CONTAINER , typename OUTPUT>
      void     buildMap ( ITERATOR , ITERATOR , CONTAINER & , OUTPUT & );

      void   initBox ( MeshNode &n , typename MESH::shared_ptr , Box & );
      void   initBox ( MeshElement &e , typename MESH::shared_ptr , Box & );

      double   expandBox ( MeshNode &n , typename MESH::shared_ptr , Box & );
      double   expandBox ( MeshElement &e , typename MESH::shared_ptr , Box & );


      void     computeBoundingBox ();
      void     decomposeDomain ();
      void     buildMaps ( std::multimap<ITERATOR1 , ITERATOR2>  &m1ToM2 ,
                           std::multimap<ITERATOR2 , ITERATOR1>  &m2toM1 );

    public:
      CoarseSearch ( typename Parameters::shared_ptr  p );
      virtual void  findNeighbors ( FirstToSecondMap  &m1ToM2 , 
                                    SecondToFirstMap  &m2ToM1 );

  };

}
}

#include "CoarseSearch.tmpl.h"
#endif
