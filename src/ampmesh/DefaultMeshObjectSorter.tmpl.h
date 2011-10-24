namespace AMP { 
namespace Mesh {

  template <typename ITERATOR>
  DefaultMeshObjectSorter<ITERATOR>::DefaultMeshObjectSorter ( ITERATOR start , ITERATOR end )
  {
    d_Size = 0;
    bool resetFirstObject = true;
    while ( start != end )
    {
      if ( start->isOwned() )
      {
        if ( resetFirstObject )
        {
          resetFirstObject = false;
          d_FirstObject = start->globalID();
        }
        else
        {
          d_FirstObject = std::min ( start->globalID() , d_FirstObject );
        }
        d_Size++;
      }
      start++;
    }
  }

  template <typename ITERATOR>
  DefaultMeshObjectSorter<ITERATOR>::~DefaultMeshObjectSorter ()
  {
  }

  template <typename ITERATOR>
  int DefaultMeshObjectSorter<ITERATOR>::operator [] ( int i ) const
  {
    return i;
  }

  template <typename ITERATOR>
  size_t DefaultMeshObjectSorter<ITERATOR>::size () const
  {
    return d_Size;
  }

}
}
