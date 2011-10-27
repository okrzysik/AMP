
#include "utils/Utilities.h"

namespace AMP {
namespace Mesh {

  template <typename MESH , typename ITERATOR_FACE>
  void ComputeNormals<MESH,ITERATOR_FACE>::compute ( ITERATOR_FACE begin , ITERATOR_FACE end , NormalMap &ans , typename MESH::shared_ptr mesh , NormalType t )
  {
    std::map<size_t , double>   numNorms;
    while ( begin != end )
    {
      typename MESH::Node n1 = mesh->getNode ( begin->getNodeID(0) );
      typename MESH::Node n2 = mesh->getNode ( begin->getNodeID(1) );
      typename MESH::Node n3 = mesh->getNode ( begin->getNodeID(2) );

      double x[3];
      double y[3];
      x[0] = n2.x() - n1.x();
      x[1] = n2.y() - n1.y();
      x[2] = n2.z() - n1.z();
      y[0] = n3.x() - n1.x();
      y[1] = n3.y() - n1.y();
      y[2] = n3.z() - n1.z();
      Normal  thisNorm;
      thisNorm.x[0] = x[1]*y[2] - x[2]*y[1];
      thisNorm.x[1] = x[0]*y[2] - x[2]*y[0];
      thisNorm.x[2] = x[0]*y[1] - x[1]*y[0];
      if ( t == Facet )
      {
        numNorms[begin->globalID()] = 0.0;
        ans[begin->globalID()] = thisNorm;
      }
      else
      {
        AMP_ASSERT ( t == Node );
        for ( size_t i = 0 ; i != begin->numNodes() ; i++ )
        {
          if ( ans.find ( begin->getNodeID(i) ) == ans.end() )
          {
            numNorms[begin->getNodeID(i)] = 0.;
          }
          double &curDenom = numNorms[begin->getNodeID ( i )];
          Normal &curNormal = ans[begin->getNodeID ( i )];
          curNormal *= curDenom;
          curNormal += thisNorm;
          curDenom += 1.0;
          curNormal /= curDenom;
        }
      }
      begin++;
    }
  }

}
}

