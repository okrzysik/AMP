namespace AMP { 
namespace Mesh {

  template <typename PTR_TYPE>
  size_t LibMeshElementTmpl<PTR_TYPE>::procID () const
  {
    return d_Elem->processor_id();
  }

  template <typename PTR_TYPE>
  LibMeshElementTmpl<PTR_TYPE>::LibMeshElementTmpl () 
    : d_Elem ( 0 ) 
  {
  }

  template <typename PTR_TYPE>
  LibMeshElementTmpl<PTR_TYPE>::LibMeshElementTmpl ( PTR_TYPE e ) 
    : d_Elem ( e ) 
  {
  }

  template <typename PTR_TYPE>
  LibMeshElementTmpl<PTR_TYPE>::LibMeshElementTmpl ( const LibMeshElementTmpl &rhs ) 
  { 
    LibMeshElementTmpl *nonConst = const_cast<LibMeshElementTmpl *> ( &rhs );
    d_Elem = nonConst->d_Elem; 
  }

  template <typename PTR_TYPE>
  PTR_TYPE  LibMeshElementTmpl<PTR_TYPE>::getPtr () 
  { 
    return d_Elem; 
  }

  template <typename PTR_TYPE>
  void LibMeshElementTmpl<PTR_TYPE>::nonConstConstruct ( const LibMeshElementTmpl rhs ) 
  { 
    d_Elem = rhs.d_Elem; 
  }

  template <typename PTR_TYPE>
  size_t   LibMeshElementTmpl<PTR_TYPE>::globalID () const 
  { 
    return d_Elem->id(); 
  }

  template <typename PTR_TYPE>
  bool  LibMeshElementTmpl<PTR_TYPE>::operator < ( const MeshObject &rhs ) const 
  { 
    return globalID() < rhs.globalID(); 
  }

  template <typename PTR_TYPE>
  bool  LibMeshElementTmpl<PTR_TYPE>::operator == ( const LibMeshElementTmpl &rhs ) const 
  { 
    return globalID() == rhs.globalID(); 
  }

  template <typename PTR_TYPE>
  ::Elem &LibMeshElementTmpl<PTR_TYPE>::getElem ()       
  { 
    return *d_Elem; 
  }

  template <typename PTR_TYPE>
  const ::Elem &LibMeshElementTmpl<PTR_TYPE>::getElem () const 
  { 
    return *d_Elem; 
  }

  template <typename PTR_TYPE>
  size_t   LibMeshElementTmpl<PTR_TYPE>::numNodes () const 
  { 
    return d_Elem->n_nodes(); 
  }

  template <typename PTR_TYPE>
  size_t   LibMeshElementTmpl<PTR_TYPE>::numSides () const 
  { 
    return d_Elem->n_sides(); 
  }

  template <typename PTR_TYPE>
  bool     LibMeshElementTmpl<PTR_TYPE>::hasNeighbor ( size_t k ) const 
  { 
    return d_Elem->neighbor(k) != 0; 
  }

  template <typename PTR_TYPE>
  LibMeshElementTmpl< AutoPtr < ::Elem> >  LibMeshElementTmpl<PTR_TYPE>::getSide ( size_t k , bool proxy ) 
  { 
    return LibMeshElementTmpl< AutoPtr < ::Elem> > ( d_Elem->build_side ( k , proxy ) ); 
  }

  template <typename PTR_TYPE>
  LibMeshPoint  LibMeshElementTmpl<PTR_TYPE>::getPoint ( size_t k ) 
  { 
    return LibMeshPoint ( d_Elem->point(k) ); 
  }

  template <typename PTR_TYPE>
  size_t   LibMeshElementTmpl<PTR_TYPE>::getNodeID ( size_t k ) const 
  { 
    AMP_INSIST ( k < numNodes() , "getNodeID(i), i >= numNodes()" );
    return d_Elem->node (k); 
  }

  template <typename PTR_TYPE>
  ElementType   LibMeshElementTmpl<PTR_TYPE>::getType () const 
  { 
    return d_Elem->type(); 
  }

  template <typename PTR_TYPE>
  bool    LibMeshElementTmpl<PTR_TYPE>::isOwned () const
  {
    return libMesh::processor_id() == d_Elem->processor_id();
  }

  template <typename PTR_TYPE>
  double  LibMeshElementTmpl<PTR_TYPE>::volume () const
  {
    return d_Elem->volume();
  }

}
}


