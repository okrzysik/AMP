namespace AMP {
namespace LinearAlgebra {

  template <bool AFFINE_MAP>
  unsigned int VectorEntryMap<AFFINE_MAP>::numDOFsPerObject ()
  {
    return d_Variable->DOFsPerObject ();
  }

  template <bool AFFINE_MAP>
  CommunicationList::shared_ptr & VectorEntryMap<AFFINE_MAP>::getCommunicationList ()
  {
    return d_CommList;
  }

  template <bool AFFINE_MAP>
  unsigned int   VectorEntryMap<AFFINE_MAP>::firstElement () 
  { 
    return d_FirstRow;
  }

  template <bool AFFINE_MAP>
  unsigned int   VectorEntryMap<AFFINE_MAP>::endElement () 
  { 
    return this->firstElement() + this->numLocalElements(); 
  }

  template <bool AFFINE_MAP>
  size_t  VectorEntryMap<AFFINE_MAP>::ComputeFirstDOF ()
  {
    size_t  retVal;
    if ( AFFINE_MAP )
    {
      retVal = d_SortedObjects->getFirstObject() * d_Variable->DOFsPerObject();
    }
    else
    {
      d_Comm.sumScan((int*)&d_NumRows,(int*)&retVal,1);
      retVal -= d_NumRows;
    }
    return retVal;
  }

  template <bool AFFINE_MAP>
  VectorEntryMap<AFFINE_MAP>::VectorEntryMap ( Parameters::shared_ptr ptr )
    : d_CommList ( ptr->d_CommList )
    , d_SortedObjects ( ptr->d_SortedObjects )
    , d_Variable ( ptr->d_Variable )
    , d_Comm ( ptr->d_Comm )
    , d_NumRows ( AFFINE_MAP ? ptr->d_CommList->numLocalRows() : 0 )
    , d_FirstRow ( AFFINE_MAP ? ptr->d_CommList->getStartGID() : 0 )
  {
    if ( AFFINE_MAP )
    {
      d_TotalSize = d_Comm.sumReduce(d_NumRows);
    }
  }

  template <bool AFFINE_MAP>
  VectorEntryMap<AFFINE_MAP>::~VectorEntryMap ()
  {
  }

  template <bool AFFINE_MAP>
  unsigned int  VectorEntryMap<AFFINE_MAP>::numLocalElements ()
  {
    return d_NumRows;
  }

  template <bool AFFINE_MAP>
  size_t VectorEntryMap<AFFINE_MAP>::getGlobalID ( size_t obj_id , size_t dof ) const
  {
    size_t  retVal;
    if ( AFFINE_MAP )
    {
      retVal = (*d_SortedObjects)[obj_id] * d_Variable->DOFsPerObject() + dof;
    }
    else
    {
      size_t  sort_obj_id = (*d_SortedObjects)[obj_id];
      retVal = d_NonAffineMap.find ( sort_obj_id )->second + dof;
    }
    return retVal;
  }

  template <bool AFFINE_MAP>
  void VectorEntryMap<AFFINE_MAP>::getLocalID ( size_t gid , size_t &lid , size_t &dof ) const
  {
    if ( AFFINE_MAP )
    {
      lid = gid / d_Variable->DOFsPerObject ();
      dof = gid % d_Variable->DOFsPerObject ();
    }
    else
    {
      AMP_ERROR( "Not yet implemented" );
    }
  }

  template <bool AFFINE_MAP>
  template <typename DOF_FIELD , typename ITERATOR>
  void VectorEntryMap<AFFINE_MAP>::computeNonAffineMap ( ITERATOR begin , ITERATOR end )
  {
    d_NumRows = 0;
    if ( !AFFINE_MAP )
    {
      while ( begin != end )
      {
        size_t numDofs = DOF_FIELD ( begin );
        d_NonAffineMap [ begin->globalID() ] = numDofs;
        d_NumRows += numDofs;
        begin++;
      }
      d_FirstRow = ComputeFirstDOF ();
      d_TotalSize = d_Comm.sumReduce(d_NumRows);
    }
  }

  template <bool AFFINE_MAP>
  unsigned int VectorEntryMap<AFFINE_MAP>::numGlobalElements()
  {
    return d_TotalSize;
  }

  template <bool AFFINE_MAP>
  void VectorEntryMap<AFFINE_MAP>::setCommunicationList ( CommunicationList::shared_ptr in )
  {
    d_CommList = in;
  }

}
}

