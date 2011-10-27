#ifndef  included_AMP_MeshOperations
#define  included_AMP_MeshOperations

#include "MeshVariable.h"
#include "vectors/Vector.h"

namespace AMP { 
namespace Mesh {

  template <typename MESH>
  class  DeformMesh
  {
    protected:
      AMP::LinearAlgebra::Vector::shared_ptr          d_Vector;
      typename MESH::shared_ptr   d_Mesh;
      DOFMap::shared_ptr          d_DOFMap;
      double                      d_Direction;

    public:
      typedef typename MESH::NodeIterator   iterator;

      DeformMesh ()
      {
        AMP_ERROR( "Deform requires a vector" );
      }

      DeformMesh ( AMP::LinearAlgebra::Vector::shared_ptr vec , typename MESH::shared_ptr m , double scale )
          : d_Vector ( vec )
          , d_Mesh ( m )
          , d_Direction ( scale )
      {
        AMP_ASSERT ( d_Vector->getVariable()->isA<Nodal3VectorVariable> () );
        d_DOFMap = d_Mesh->getDOFMap ( d_Vector->getVariable() );
      }

      void operator () ( iterator curNode )
      {
        std::vector<unsigned int>  dofs , empty;
        d_DOFMap->getDOFs ( *curNode , dofs , empty );
        if ( d_Vector->containsGlobalElement (dofs[0]) ) 
        {
          curNode->x() += d_Direction*d_Vector->getValueByGlobalID ( dofs[0] );
          curNode->y() += d_Direction*d_Vector->getValueByGlobalID ( dofs[1] );
          curNode->z() += d_Direction*d_Vector->getValueByGlobalID ( dofs[2] );
        }
      }
  };



}
}

#endif
