
#include "utils/FC.h"
#define  DGSEV  FC_GLOBAL_ ( dgesv , DGESV )

extern "C" 
{
  void DGSEV ( int* n, int* nrhs, double* a, int* lda, int* ipiv,
                double* b, int* ldb, int* info );
}


namespace AMP {
namespace Mesh {


  template <typename MESH , typename NORMAL_MAP , typename ITERATOR1 , typename ITERATOR2>
  FineSearch<MESH,NORMAL_MAP,ITERATOR1,ITERATOR2>::FineSearch ( typename Parameters::shared_ptr p ) 
    : d_Parameters ( p )
  {
  }

  template <typename MESH , typename NORMAL_MAP , typename ITERATOR1 , typename ITERATOR2>
  void FineSearch<MESH,NORMAL_MAP,ITERATOR1,ITERATOR2>::J::GEPP ( Vector &x )
  {
    int n = 3;
    int piv[3];
    int info;
    DGSEV ( &n , &n , mat , &n , piv , x.x , &n , &info );
  }

  template <typename MESH , typename NORMAL_MAP , typename ITERATOR1 , typename ITERATOR2>
  void FineSearch<MESH,NORMAL_MAP,ITERATOR1,ITERATOR2>::interpolatePosition ( MeshNode &n , MeshElement &e , Normal &norm , typename MESH::shared_ptr mesh , Vector &input , Vector &answer , J &deriv , Vector &interp )
  {
    double N[4];
    double dN[8];
    N[0] = .25*(1.-input.x[0])*(1.-input.x[1]);
    N[1] = .25*(1.-input.x[0])*(1.+input.x[1]);
    N[2] = .25*(1.+input.x[0])*(1.+input.x[1]);
    N[3] = .25*(1.+input.x[0])*(1.-input.x[1]);

    dN[0] = -.25 * (1.-input.x[1]);
    dN[1] = -.25 * (1.-input.x[0]);
    dN[2] = -.25 * (1.+input.x[1]);
    dN[3] =  .25 * (1.-input.x[0]);
    dN[4] =  .25 * (1.+input.x[1]);
    dN[5] =  .25 * (1.+input.x[0]);
    dN[6] =  .25 * (1.-input.x[1]);
    dN[7] = -.25 * (1.+input.x[0]);

    Vector project, Nxeta1, Nxeta2 , node, node2;
    interp = Nxeta1 = Nxeta2 = 0.;
    project = norm;
    project *= input.x[2];
    node.x[0] = n.x();
    node.x[1] = n.y();
    node.x[2] = n.z();

    for ( size_t i = 0 ; i != 4 ; i++ )
    {
      Node curNode = mesh->getNode ( e.getNodeID ( i ) );

      interp.x[0] += N[i] * curNode.x();
      interp.x[1] += N[i] * curNode.y();
      interp.x[2] += N[i] * curNode.z();

      Nxeta1.x[0] += dN[2*i]   * curNode.x();
      Nxeta1.x[1] += dN[2*i]   * curNode.y();
      Nxeta1.x[2] += dN[2*i]   * curNode.z();
      Nxeta2.x[0] += dN[2*i+1] * curNode.x();
      Nxeta2.x[1] += dN[2*i+1] * curNode.y();
      Nxeta2.x[2] += dN[2*i+1] * curNode.z();
    }

    answer = node;
    answer += project;
    answer -= interp;

    for ( size_t i = 0 ; i != 3 ; i++ )
    {
      deriv.mat[0 + i] = -Nxeta1.x[i];
      deriv.mat[3 + i] = -Nxeta2.x[i];
      deriv.mat[6 + i] = norm.x[i];
    }
  }

  template <typename MESH , typename NORMAL_MAP , typename ITERATOR1 , typename ITERATOR2>
  void FineSearch<MESH,NORMAL_MAP,ITERATOR1,ITERATOR2>::findIntersection ( MeshNode &n , MeshElement &e , Normal &norm , typename MESH::shared_ptr mesh , Vector &curGuess , Vector &intersection )
  {
    Vector curAnswer;
    curGuess = 0.;
    J      curJacob;
    interpolatePosition ( n , e , norm , mesh , curGuess , curAnswer , curJacob , intersection );
    int iter = 0;
    while ( curAnswer.norm() > 1.e-12 )
    {
      if (( fabs(curGuess.x[0]) > 5. ) ||
          ( fabs(curGuess.x[1]) > 5. ) ||
          ( fabs(curGuess.x[2]) > 5. )) // diverged
      {
        break;
      }
      if ( iter == 50 )
      {
        std::runtime_error ( "Stalled computation" );
      }
      iter++;
      curJacob.GEPP ( curAnswer );
      curGuess -= curAnswer;
      interpolatePosition ( n , e , norm , mesh , curGuess , curAnswer , curJacob , intersection );
    }
  }

  template <typename MESH , typename NORMAL_MAP , typename ITERATOR1 , typename ITERATOR2>
  bool FineSearch<MESH,NORMAL_MAP,ITERATOR1,ITERATOR2>::isObject1InContactWithObject2 ( MeshNode &n , MeshElement &e , typename MESH::shared_ptr mesh )
  {
    // Assuming quads at the moment...
    Vector local_projection;
    Vector global_projection;


    findIntersection ( n , e , d_Parameters->d_Normals[n.globalID()] , mesh , local_projection , global_projection );

    // The imprecision of the solve forces a little fudge factor
    if ( local_projection.x[0] < -1.0001 )
      return false;
    if ( local_projection.x[1] < -1.0001 )
      return false;
    if ( local_projection.x[0] > 1.0001 )
      return false;
    if ( local_projection.x[1] > 1.0001 )
      return false;

    d_Parameters->d_Intersections.insert ( std::pair<size_t,Vector> ( n.globalID() , global_projection ) );

    return true;
  }

  template <typename MESH , typename NORMAL_MAP , typename ITERATOR1 , typename ITERATOR2>
  bool FineSearch<MESH,NORMAL_MAP,ITERATOR1,ITERATOR2>::isObject1InContactWithObject2 ( MeshElement &e , MeshNode &n , MeshSharedPtr mesh )
  {
    return isObject1InContactWithObject2 ( n , e , mesh );
  }

  template <typename MESH , typename NORMAL_MAP , typename ITERATOR1 , typename ITERATOR2>
  bool FineSearch<MESH,NORMAL_MAP,ITERATOR1,ITERATOR2>::isObject1InContactWithObject2 ( MeshElement &e1 , MeshElement &e2 , typename MESH::shared_ptr obj1_mesh )
  {
    bool answer = false;
    for ( size_t i = 0 ; i != e1.numNodes() ; i++ )
    {
      typename MESH::Node node = obj1_mesh->getNode ( e1.getNodeID ( i ) );
      answer = answer || isObject1InContactWithObject2 ( node , e2 , obj1_mesh );
    }

    return answer;
  }

  template <typename MESH , typename NORMAL_MAP , typename ITERATOR1 , typename ITERATOR2>
  template <typename CONTAINER>
  void FineSearch<MESH,NORMAL_MAP,ITERATOR1,ITERATOR2>::refineMap ( CONTAINER &in , CONTAINER &out , typename MESH::shared_ptr obj1_mesh )
  {
    typedef std::pair<typename CONTAINER::key_type ,
                      typename CONTAINER::mapped_type >  local_pair;

    typename CONTAINER::iterator curPair = in.begin();
    while ( curPair != in.end() )
    {
      typename CONTAINER::key_type::value_type  obj1 = *(curPair->first);
      typename CONTAINER::mapped_type::value_type  obj2 = *(curPair->second);

      if ( isObject1InContactWithObject2 ( obj1 , obj2 , obj1_mesh ) )
      {
        out.insert ( local_pair ( curPair->first , curPair->second ) );
      }
      curPair++;
    }
  }

  template <typename MESH , typename NORMAL_MAP , typename ITERATOR1 , typename ITERATOR2>
  template <typename IN_CONTAINER , typename OUT_CONTAINER>
  void FineSearch<MESH,NORMAL_MAP,ITERATOR1,ITERATOR2>::buildReverseMap ( IN_CONTAINER &in , OUT_CONTAINER &out )
  {
    typedef std::pair<typename OUT_CONTAINER::key_type ,
                      typename OUT_CONTAINER::mapped_type >  local_pair;

    typename IN_CONTAINER::iterator  curPair = in.begin();
    while ( curPair != in.end() )
    {
      out.insert ( local_pair ( curPair->second , curPair->first ) );
      curPair++;
    }
  }

  template <typename MESH , typename NORMAL_MAP , typename ITERATOR1 , typename ITERATOR2>
  void FineSearch<MESH,NORMAL_MAP,ITERATOR1,ITERATOR2>::refineNeighbors ( 
                         FirstToSecondMap &inM1ToM2 , 
                         SecondToFirstMap &inM2ToM1 , 
                         FirstToSecondMap &outM1ToM2 , 
                         SecondToFirstMap &outM2ToM1 )
  {
    refineMap ( inM1ToM2 , outM1ToM2 , d_Parameters->d_mesh1 );
    buildReverseMap ( outM1ToM2 , outM2ToM1 );
  }

}
}

