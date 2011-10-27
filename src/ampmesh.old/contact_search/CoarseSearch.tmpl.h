#include "mesh.h"
#include <algorithm>

namespace AMP {
namespace Mesh {

  template <typename MESH , typename ITERATOR1 , typename ITERATOR2>
  CoarseSearch<MESH , ITERATOR1,ITERATOR2>::CoarseSearch ( typename Parameters::shared_ptr p )
     : d_Parameters ( p )
  {
  }

  template <typename MESH , typename ITERATOR1 , typename ITERATOR2>
  void  CoarseSearch<MESH , ITERATOR1,ITERATOR2>::
                    findNeighbors ( std::multimap<ITERATOR1 , ITERATOR2>  &m1ToM2 ,
                                    std::multimap<ITERATOR2 , ITERATOR1>  &m2ToM1 )
  {
    if ( d_Parameters->d_beginM1 == d_Parameters->d_endM1 ) return;
    if ( d_Parameters->d_beginM2 == d_Parameters->d_endM2 ) return;
    computeBoundingBox ();
    decomposeDomain ();
    buildMaps ( m1ToM2 , m2ToM1 );
  }

  template <typename MESH , typename ITERATOR1 , typename ITERATOR2>
  void CoarseSearch<MESH,ITERATOR1,ITERATOR2>::initBox ( MeshNode &n , typename MESH::shared_ptr , Box &box )
  {
    box.first.x = box.second.x = n.x();
    box.first.y = box.second.y = n.y();
    box.first.z = box.second.z = n.z();
  }


  template <typename MESH , typename ITERATOR1 , typename ITERATOR2>
  void CoarseSearch<MESH,ITERATOR1,ITERATOR2>::initBox ( MeshElement &e , typename MESH::shared_ptr mesh , Box &box )
  {
    typename MESH::Node node = mesh->getNode ( e.getNodeID ( 0 ) );
    initBox ( node , mesh , box );
    for ( size_t i = 1 ; i < e.numNodes() ; i++ )
    {
      node = mesh->getNode ( e.getNodeID ( i ) );
      expandBox ( node , mesh , box );
    }
  }

  template <typename MESH , typename ITERATOR1 , typename ITERATOR2>
  double CoarseSearch<MESH,ITERATOR1,ITERATOR2>::expandBox ( MeshNode &n , typename MESH::shared_ptr , Box &box )
  {
    box.first.x = std::min ( n.x() , box.first.x );
    box.first.y = std::min ( n.y() , box.first.y );
    box.first.z = std::min ( n.z() , box.first.z );

    box.second.x = std::max ( n.x() , box.second.x );
    box.second.y = std::max ( n.y() , box.second.y );
    box.second.z = std::max ( n.z() , box.second.z );

    return 0;
  }

  template <typename MESH , typename ITERATOR1 , typename ITERATOR2>
  double CoarseSearch<MESH,ITERATOR1,ITERATOR2>::expandBox ( MeshElement &e , typename MESH::shared_ptr mesh , Box &box )
  {
    Box  localBox;
    typename MESH::Node node = mesh->getNode ( e.getNodeID ( 0 ) );
    initBox ( node , mesh , localBox );
    for ( size_t i = 0 ; i < e.numNodes() ; i++ )
    {
      node = mesh->getNode ( e.getNodeID ( i ) );
      expandBox ( node , mesh , box );
      expandBox ( node , mesh , localBox );
    }

    double answer = localBox.second.x - localBox.first.x;
    answer = std::max ( localBox.second.y - localBox.first.y , answer );
    answer = std::max ( localBox.second.z - localBox.first.z , answer );
    return answer;
  }

  template <typename MESH , typename ITERATOR1 , typename ITERATOR2>
  template <typename ITERATOR>
  double  CoarseSearch<MESH,ITERATOR1,ITERATOR2>::growBoundingBox ( ITERATOR begin , ITERATOR end , typename MESH::shared_ptr mesh )
  {
    double answer = 0;
    while ( begin != end )
    {
      answer = std::max ( answer , expandBox ( *begin , mesh , d_BoundingBox ) );
      begin++;
    }

    return answer;
  }

  template <typename MESH , typename ITERATOR1 , typename ITERATOR2>
  void  CoarseSearch<MESH,ITERATOR1,ITERATOR2>::computeBoundingBox ()
  {
    double mesh1MaxFeatureSize;
    double mesh2MaxFeatureSize;

    typename ITERATOR1::value_type  firstThing = *(d_Parameters->d_beginM1);
    initBox ( firstThing , d_Parameters->d_mesh1 , d_BoundingBox );

    mesh1MaxFeatureSize = growBoundingBox ( d_Parameters->d_beginM1 , 
                                            d_Parameters->d_endM1 , 
                                            d_Parameters->d_mesh1 );

    mesh2MaxFeatureSize = growBoundingBox ( d_Parameters->d_beginM2 , 
                                            d_Parameters->d_endM2 , 
                                            d_Parameters->d_mesh2 );

    d_SmallestPartitionSide = 1.01 * std::max ( mesh1MaxFeatureSize , mesh2MaxFeatureSize );
    AMP_ASSERT ( d_SmallestPartitionSide > 0 );
  }

  template <typename MESH , typename ITERATOR1 , typename ITERATOR2>
  void  CoarseSearch<MESH , ITERATOR1,ITERATOR2>::decomposeDomain ()
  {
    sortObjects ( d_Parameters->d_beginM1 ,
                  d_Parameters->d_endM1 ,
                  d_Parameters->d_mesh1 ,
                  d_Mesh1Decomposition );

    sortObjects ( d_Parameters->d_beginM2 ,
                  d_Parameters->d_endM2 ,
                  d_Parameters->d_mesh1 ,
                  d_Mesh2Decomposition );
  }

  template <typename MESH , typename ITERATOR1 , typename ITERATOR2>
  template <typename ITERATOR , typename CONTAINER , typename OUTPUT>
  void  CoarseSearch<MESH,ITERATOR1,ITERATOR2>::buildMap ( ITERATOR begin , ITERATOR end , CONTAINER &input , OUTPUT &output )
  {
    typedef std::pair < typename OUTPUT::key_type , typename OUTPUT::mapped_type > local_pair;

    while ( begin != end )
    {
      for ( size_t i = 0 ; i != 27 ; i++ )
      {
        Index curNdx = begin->first.getNeighbor (i);
        typename CONTAINER::iterator lb = input.lower_bound ( curNdx );
        typename CONTAINER::iterator ub = input.upper_bound ( curNdx );
        while ( lb != ub )
        {
          typename OUTPUT::iterator olb = output.lower_bound ( begin->second );
          typename OUTPUT::iterator oub = output.upper_bound ( begin->second );
          bool found = false;
          while ( olb != oub )
          {
            if ( olb->second == lb->second )
            {
              found = true;
              break;
            }
            olb++;
          }
          if ( !found )
          {
            local_pair t_pair ( begin->second , lb->second );
            output.insert ( t_pair );
          }
          lb++;
        }
      }
      begin++;
    }
  }


  template <typename MESH , typename ITERATOR1 , typename ITERATOR2>
  void   CoarseSearch<MESH , ITERATOR1,ITERATOR2>::buildMaps ( std::multimap<ITERATOR1 , ITERATOR2> &m1ToM2 , std::multimap<ITERATOR2,ITERATOR1> &m2ToM1 )
  {
    buildMap ( d_Mesh1Decomposition.begin() ,
               d_Mesh1Decomposition.end() ,
               d_Mesh2Decomposition ,
               m1ToM2 );

    buildMap ( d_Mesh2Decomposition.begin() ,
               d_Mesh2Decomposition.end() ,
               d_Mesh1Decomposition ,
               m2ToM1 );
  }

  template <typename MESH , typename ITERATOR1 , typename ITERATOR2>
  typename CoarseSearch<MESH,ITERATOR1,ITERATOR2>::Index 
         CoarseSearch<MESH,ITERATOR1,ITERATOR2>::computeIndex ( MeshNode &n )
  {
    Index answer;

    answer.x = n.x()/d_SmallestPartitionSide;
    answer.y = n.y()/d_SmallestPartitionSide;
    answer.z = n.z()/d_SmallestPartitionSide;
    if ( n.x() < 0 ) answer.x--;
    if ( n.y() < 0 ) answer.y--;
    if ( n.z() < 0 ) answer.z--;

    return answer;
  }

  template <typename MESH , typename ITERATOR1 , typename ITERATOR2>
  template <typename ITERATOR>
  void CoarseSearch<MESH,ITERATOR1,ITERATOR2>::putInContainer ( MeshNode &n , ITERATOR i , std::multimap<Index,ITERATOR> &container , typename MESH::shared_ptr )
  {
    size_t j = container.size();
    container.insert ( std::pair<Index,ITERATOR> ( computeIndex ( n ) , i ) );
    AMP_ASSERT ( j+1 == container.size() );
  }

  template <typename MESH , typename ITERATOR1 , typename ITERATOR2>
  template <typename ITERATOR>
  void CoarseSearch<MESH,ITERATOR1,ITERATOR2>::putInContainer ( MeshElement &e , ITERATOR iter , std::multimap<Index,ITERATOR> &container , typename MESH::shared_ptr mesh )
  {
    for ( size_t i = 0 ; i != e.numNodes() ; i++ )
    {
      typename MESH::Node node = mesh->getNode ( e.getNodeID ( i ) );
      putInContainer ( node , iter , container , mesh );
    }
  }

  template <typename MESH , typename ITERATOR1 , typename ITERATOR2>
  template <typename ITERATOR>
  void  CoarseSearch<MESH,ITERATOR1,ITERATOR2>::sortObjects ( ITERATOR begin , ITERATOR end , typename MESH::shared_ptr mesh , std::multimap<Index,ITERATOR>  &decomp )
  {
    while ( begin != end )
    {
      putInContainer ( *begin , begin , decomp , mesh );
      begin++;
    }
  }


}
}

