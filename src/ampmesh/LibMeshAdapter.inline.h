#include "vectors/MultiVector.h"
#include "vectors/VectorSelector.h"

namespace AMP { 
namespace Mesh {


  inline
  LibMeshAdapter::~LibMeshAdapter ()
  {
  }

  inline
  void  LibMeshAdapter::clearData ()
  {
    d_mData.clear();
  }

  inline
  AMP::LinearAlgebra::Variable::shared_ptr LibMeshAdapter::getMeshSpecificVariable ( AMP::LinearAlgebra::Variable::shared_ptr in )
  {
    if ( d_UseMeshNameInVariable )
    {
      std::stringstream new_name;
      new_name << in->getName() << "_" << d_MeshName;
      return in->cloneVariable ( new_name.str() );
    }
    else
    {
      return in;
    }
  }

  inline
  std::string     LibMeshAdapter::getMeshName () 
  { 
    return d_MeshName; 
  }

  inline
  void LibMeshAdapter::setName ( const std::string &name , bool t )
  {
    d_UseMeshNameInVariable = t;
    d_MeshName = name;
  }

  inline
  size_t   LibMeshAdapter::geometricDimension () 
  { 
    return d_libMesh->mesh_dimension(); 
  }

  inline
  void LibMeshAdapter::appendMeshNameToVariable ( AMP::LinearAlgebra::Variable::shared_ptr in )
  {
    if ( d_UseMeshNameInVariable )
    {
      std::string name = in->getName();
      size_t  spot = name.rfind ( d_MeshName );
      size_t  end = name.size() - d_MeshName.size();
      bool  appendName = false;
      // Append the mesh name if
      //  1)  The mesh name is not found in the variable name
      appendName = appendName || ( spot == std::string::npos );
      //  2)  The mesh name is not at the end of the variable
      appendName = appendName || ( spot != end );
      //  3)  The beginning of the variable name is the name of the mesh
      appendName = appendName || ( spot == 0 );
      if ( appendName )
      {
        name += "_";
        name += d_MeshName;
        in->setName ( name );
      }
    }
  }

  inline
  void    LibMeshAdapter::registerVectorAsData ( AMP::LinearAlgebra::Vector::shared_ptr vec , const std::string &s )
  {
    std::string t = s;
    AMP::LinearAlgebra::Variable::shared_ptr var = vec->getVariable();
    if ( s == "" )
    {
      t = var->getName();
    }
    if ( vec->isA<AMP::LinearAlgebra::MultiVector> () )
    {
      AMP::LinearAlgebra::MultiVector &mv = vec->castTo<AMP::LinearAlgebra::MultiVector> ();
      for ( size_t i = 0 ; i != mv.getNumberOfSubvectors () ; i++ )
      {
        std::stringstream new_name;
        new_name << t << "_" << mv.getVector(i)->getVariable()->getName();
        registerVectorAsData ( mv.getVector(i) , new_name.str() );
      }
      return;
    }
    size_t n_dofs = vec->getVariable()->DOFsPerObject();
    if ( n_dofs > 3 )
    {
      for ( size_t i = 0 ; i != n_dofs ; i++ )
      {
        std::stringstream newname;
        newname << t << "_" << i;
        registerVectorAsData ( vec->select ( AMP::LinearAlgebra::VS_Stride ( newname.str() , i , n_dofs ) , newname.str() ) );
      }
      return;
    }

    if ( d_mData.find ( t ) != d_mData.end() )
    {
      AMP::LinearAlgebra::Vector::shared_ptr othVec = d_mData[t];
      AMP_INSIST ( othVec == vec , "Different vector with smae name already being written to Silo file." );
    }
    d_mData[t] = vec;
  }

  inline
  LibMeshAdapter::DataIterator      LibMeshAdapter::beginData ()
  {
    return d_mData.begin();
  }

  inline
  LibMeshAdapter::DataIterator      LibMeshAdapter::endData ()
  {
    return d_mData.end();
  }

  inline
  size_t  LibMeshAdapter::sizeData ()
  {
    return d_mData.size();
  }

  inline
  LibMeshAdapter::ElementIterator   LibMeshAdapter::beginElement()
  {
    return ElementIterator ( d_libMesh->local_elements_begin() );
  }

  inline
  LibMeshAdapter::ElementIterator   LibMeshAdapter::endElement()
  {
    return ElementIterator ( d_libMesh->local_elements_end() );
  }

  inline
  LibMeshAdapter::ConstElementIterator  LibMeshAdapter::beginElement () const
  {
    return ConstElementIterator ( d_libMesh->local_elements_begin() );
  }

  inline
  LibMeshAdapter::ConstElementIterator  LibMeshAdapter::endElement ()   const
  {
    return ConstElementIterator ( d_libMesh->local_elements_end() );
  }

  inline
  LibMeshAdapter::NodeElementIterator  LibMeshAdapter::beginElementForNode ( const Node &n )
  {
    return d_NodeElemMap.beginElement ( n.globalID() );
  }

  inline
  LibMeshAdapter::NodeElementIterator  LibMeshAdapter::endElementForNode ( const Node &n )
  {
    return d_NodeElemMap.endElement ( n.globalID() );
  }

  inline
  LibMeshAdapter::Element  LibMeshAdapter::getElement ( size_t i )
  {
    return Element ( d_libMesh->elem ( i ) );
  }

  inline
  LibMeshAdapter::Node     LibMeshAdapter::getNode    ( size_t i )
  {
    return Node ( d_libMesh->node_ptr ( i ) );
  }

  inline
  LibMeshAdapter::NodeIterator   LibMeshAdapter::beginNode()
  {
    return NodeIterator ( d_libMesh->local_nodes_begin() );
  }

  inline
  LibMeshAdapter::NodeIterator   LibMeshAdapter::endNode()
  {
    return NodeIterator ( d_libMesh->local_nodes_end() );
  }

  inline
  LibMeshAdapter::NodeIterator   LibMeshAdapter::beginUsedNode ()
  {
    return NodeIterator ( d_libMesh->nodes_begin() );
  }

  inline
  LibMeshAdapter::NodeIterator   LibMeshAdapter::endUsedNode ()
  {
    return NodeIterator ( d_libMesh->nodes_end() );
  }

  inline
  LibMeshAdapter::OwnedNodeIterator   LibMeshAdapter::beginOwnedNode()
  {
    return OwnedNodeIterator ( d_libMesh->local_nodes_begin() , d_libMesh->local_nodes_end() );
  }

  inline
  LibMeshAdapter::OwnedNodeIterator   LibMeshAdapter::endOwnedNode()
  {
    return OwnedNodeIterator ( d_libMesh->local_nodes_end() , d_libMesh->local_nodes_end() );
  }

  inline
  LibMeshAdapter::ConstNodeIterator  LibMeshAdapter::beginNode () const
  {
    return ConstNodeIterator ( d_libMesh->local_nodes_begin() );
  }

  inline
  LibMeshAdapter::ConstNodeIterator  LibMeshAdapter::endNode ()   const
  {
    return ConstNodeIterator ( d_libMesh->local_nodes_end() );
  }

  inline
  LibMeshAdapter::BoundarySideIterator  LibMeshAdapter::beginSideBoundary ( short int which )
  {
    return d_BoundarySet.beginSideBoundary ( which );
  }

  inline
  LibMeshAdapter::BoundarySideIterator  LibMeshAdapter::endSideBoundary ( short int which )
  {
    return d_BoundarySet.endSideBoundary ( which );
  }

  inline
  LibMeshAdapter::BoundaryNodeIterator  LibMeshAdapter::beginBoundary ( short int which )
  {
    return d_BoundarySet.beginBoundary ( which );
  }

  inline
  LibMeshAdapter::BoundaryNodeIterator  LibMeshAdapter::endBoundary ( short int which )
  {
    return d_BoundarySet.endBoundary ( which );
  }

  inline
  LibMeshAdapter::OwnedBoundaryNodeIterator  LibMeshAdapter::beginOwnedBoundary ( short int which )
  {
    return OwnedBoundaryNodeIterator ( d_BoundarySet.beginBoundary ( which ) , d_BoundarySet.endBoundary ( which ) );
  }

  inline
  LibMeshAdapter::OwnedBoundaryNodeIterator  LibMeshAdapter::endOwnedBoundary ( short int which )
  {
    return OwnedBoundaryNodeIterator ( d_BoundarySet.endBoundary ( which ) , d_BoundarySet.endBoundary ( which ) );
  }

  inline
  size_t  LibMeshAdapter::numLocalElements ()
  {
    return d_libMesh->n_local_elem();
  }

  inline
  size_t  LibMeshAdapter::numLocalNodes ()
  {
    return d_libMesh->n_local_nodes();
  }

  inline
  size_t  LibMeshAdapter::numUsedNodes ()
  {
    return d_libMesh->n_nodes();
  }

  inline
  size_t  LibMeshAdapter::numUsedElements ()
  {
    return d_libMesh->n_elem();
  }

  inline
  size_t  LibMeshAdapter::numTotalElements ()
  {
    return d_libMesh->parallel_n_elem();
  }

  inline
  size_t  LibMeshAdapter::numTotalNodes ()
  {
    return d_libMesh->parallel_n_nodes();
  }

  inline
  ::Mesh  &LibMeshAdapter::getMesh()
  {
    return *d_libMesh;
  }

  inline
  const ::Mesh  &LibMeshAdapter::getMesh() const
  {
    return *d_libMesh;
  }

  inline
  AMP_MPI  LibMeshAdapter::getComm()
  {
    return AMP_MPI(libMesh::COMM_WORLD);
  }

  inline
  int  LibMeshAdapter::dim ()
  {
    return d_libMesh->mesh_dimension ();
  }


}
}


