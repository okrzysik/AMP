namespace AMP { 
namespace Mesh {

  inline
  size_t LibMeshNode::procID () const
  {
    return d_Node->processor_id();
  }

  inline
  LibMeshNode::LibMeshNode () 
    : MeshNode () 
    , d_Node ( 0 ) 
  {
  }

  inline
  LibMeshNode::LibMeshNode ( ::Node *n )  
    : MeshNode ()
    , d_Node ( n ) 
  {
  }

  inline
  LibMeshNode::LibMeshNode ( const LibMeshNode &rhs ) 
    : MeshNode ()
    , d_Node ( rhs.d_Node  )
  {
  }

  inline
  ::Node &LibMeshNode::getNode()       
  { 
    return *d_Node; 
  }

  inline
  const ::Node &LibMeshNode::getNode() const 
  { 
    return *d_Node; 
  }

  inline
  size_t   LibMeshNode::globalID () const 
  { 
    return d_Node->id(); 
  }

  inline
  void LibMeshNode::translate ( double x , double y , double z )
  {
    (*d_Node)(0) += x;
    (*d_Node)(1) += y;
    (*d_Node)(2) += z;
  }

  inline
  double   &LibMeshNode::x ()       
  { 
    return (*d_Node)(0); 
  }

  inline
  double   &LibMeshNode::y ()       
  { 
    return (*d_Node)(1); 
  }

  inline
  double   &LibMeshNode::z ()       
  { 
    return (*d_Node)(2); 
  }

  inline
  double    LibMeshNode::x () const 
  { 
    return (*d_Node)(0); 
  }

  inline
  double    LibMeshNode::y () const 
  { 
    return (*d_Node)(1); 
  }

  inline
  double    LibMeshNode::z () const 
  { 
    return (*d_Node)(2); 
  }


  inline
  double  LibMeshNode::operator () ( size_t t ) const 
  { 
    return (*d_Node)(t); 
  }

  inline
  double &LibMeshNode::operator () ( size_t t )       
  { 
    return (*d_Node)(t); 
  }

  inline
  bool    LibMeshNode::isOwned () const
  {
    return libMesh::processor_id() == d_Node->processor_id();
  }

}
}

