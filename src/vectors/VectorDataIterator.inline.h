
#include <exception>

#include "Vector.h"
#include "MultiVector.h"


namespace AMP {
namespace LinearAlgebra {

  inline
  VectorDataIterator VectorDataIterator::operator - ( int offset )
  {
    if ( d_Indexer )
    {
      AMP_ERROR( "Operator - not implemented for SubsetVector views" );
    }
    return operator + ( -offset );
  }

  inline
  VectorDataIterator &VectorDataIterator::operator -= ( int offset )
  {
    return operator += ( -offset );
  }


  inline
  size_t   VectorDataIterator::dbSize () const
  {
    return d_Vec->sizeOfDataBlock ( 0 );
  }

  inline
  double   & VectorDataIterator::operator * ()
  {
    return d_Block[d_CurOffset];
  }

  inline
  bool       VectorDataIterator::operator == ( const VectorDataIterator &rhs ) const
  {
    bool RetVal = true;
    RetVal &= d_Vec == rhs.d_Vec;
    RetVal &= d_CurBlock == rhs.d_CurBlock;
    RetVal &= d_CurOffset == rhs.d_CurOffset;
    return RetVal;
  }

  inline
  double & VectorDataIterator::operator [] ( int i )
  {
    VectorDataIterator t = *this + i;
    return *t;
  }


  inline
  bool       VectorDataIterator::operator != ( const VectorDataIterator &rhs ) const
  {
    return !operator == (rhs);
  }


  inline
  const double & ConstVectorDataIterator::operator [] ( int i )
  {
    ConstVectorDataIterator t = *this + i;
    return *t;
  }

  inline
  ConstVectorDataIterator ConstVectorDataIterator::operator - ( int offset )
  {
    return operator + ( -offset );
  }

  inline
  ConstVectorDataIterator &ConstVectorDataIterator::operator -= ( int offset )
  {
    return operator += ( -offset );
  }


  inline
  size_t   ConstVectorDataIterator::dbSize () const
  {
    return d_Vec->sizeOfDataBlock ( 0 );
  }

  inline
  const double    &ConstVectorDataIterator::operator * ()
  {
    return d_Block[d_CurOffset];
  }

  inline
  bool       ConstVectorDataIterator::operator == ( const ConstVectorDataIterator &rhs ) const
  {
    bool RetVal = true;
    RetVal &= d_Vec == rhs.d_Vec;
    RetVal &= d_CurBlock == rhs.d_CurBlock;
    RetVal &= d_CurOffset == rhs.d_CurOffset;
    return RetVal;
  }

  inline
  bool       ConstVectorDataIterator::operator != ( const ConstVectorDataIterator &rhs ) const
  {
    return !operator == (rhs);
  }

}
}

