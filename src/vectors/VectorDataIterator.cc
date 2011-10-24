#include "VectorDataIterator.h"

namespace AMP {
namespace LinearAlgebra {

  VectorDataIterator  VectorDataIterator::operator + ( int offset )
  {
    VectorDataIterator ans ( *this );
    if ( d_Indexer ) // Not constant time
    {
      if ( offset < 0 )
      {
        for ( int i = 0 ; i != offset ; i++ )
        {
         --ans;
        }
      }
      else
      {
        for ( int i = 0 ; i != offset ; i++ )
        {
          ++ans;
        }
      }
    }
    else
    {
      if ( offset > 0 )
      {
        advance ( offset , ans );
      }
      if ( offset < 0 )
      {
        recede ( offset , ans );
      }
    }
    return ans;
  }
  VectorDataIterator::VectorDataIterator ()
        : d_MultiVector ( 0 )
        , d_Vec ( 0 )
        , d_CurBlock ( 0 )
        , d_CurOffset ( 0 )
       , d_Indexer ( VectorIndexer::shared_ptr () )
  {
    d_Block = 0;
  }

  VectorDataIterator::VectorDataIterator ( const VectorDataIterator &rhs )
        : d_MultiVector ( rhs.d_MultiVector )
        , d_Vec ( rhs.d_Vec )
        , d_CurBlock ( rhs.d_CurBlock )
        , d_CurOffset ( rhs.d_CurOffset )
        , d_Indexer ( rhs.d_Indexer )
  {
    d_Block = d_Vec->getRawDataBlock<double> ( 0 );
  }


  VectorDataIterator::VectorDataIterator ( Vector *p , int block , int offset , VectorIndexer::shared_ptr  ndx )
        : d_MultiVector ( 0 )
        , d_Vec ( p )
        , d_CurBlock ( block )
        , d_CurOffset ( offset )
        , d_Indexer ( ndx )
  {
    d_Block = d_Vec->getRawDataBlock<double> ( 0 );
  }

  VectorDataIterator::VectorDataIterator ( VectorDataIterator i , MultiVector *p , size_t blockNum )
        : d_MultiVector ( p )
        , d_Vec ( i.d_Vec )
        , d_CurBlock ( i.d_CurBlock + blockNum )
        , d_CurOffset ( i.d_CurOffset )
        , d_Indexer ( i.d_Indexer )
  {
    d_Block = d_Vec->getRawDataBlock<double> ( 0 );
  }

  ConstVectorDataIterator::ConstVectorDataIterator ( ConstVectorDataIterator i , const MultiVector *const p , size_t blockNum )
        : d_MultiVector ( const_cast<MultiVector *> ( p ) )
        , d_Vec ( i.d_Vec )
        , d_CurBlock ( i.d_CurBlock + blockNum )
        , d_CurOffset ( i.d_CurOffset )
        , d_Indexer ( i.d_Indexer )
  {
    d_Block = d_Vec->getRawDataBlock<double> ( 0 );
  }


  ConstVectorDataIterator::ConstVectorDataIterator ()
       : d_MultiVector ( 0 )
       , d_Vec ( 0 )
       , d_CurBlock ( 0 )
       , d_CurOffset ( 0 )
       , d_Indexer ( VectorIndexer::shared_ptr () )
  {
    d_Block = 0;
  }

  ConstVectorDataIterator::ConstVectorDataIterator ( const VectorDataIterator &rhs )
       : d_MultiVector ( rhs.d_MultiVector )
       , d_Vec ( rhs.d_Vec )
       , d_CurBlock ( rhs.d_CurBlock )
       , d_CurOffset ( rhs.d_CurOffset )
       , d_Indexer ( rhs.d_Indexer )
  {
    d_Block = d_Vec->getRawDataBlock<double> ( 0 );
  }

  ConstVectorDataIterator::ConstVectorDataIterator ( const ConstVectorDataIterator &rhs )
       : d_MultiVector ( rhs.d_MultiVector )
       , d_Vec ( rhs.d_Vec )
       , d_CurBlock ( rhs.d_CurBlock )
       , d_CurOffset ( rhs.d_CurOffset )
       , d_Indexer ( rhs.d_Indexer )
  {
    d_Block = d_Vec->getRawDataBlock<double> ( 0 );
  }

  ConstVectorDataIterator  &ConstVectorDataIterator::operator++ ()
  {
    if ( d_Indexer )
    {
      advance ( d_Indexer->getIncrement ( d_CurOffset ) , *this );
      return *this;
    }
    AMP_ASSERT ( d_CurOffset < dbSize() );
    d_CurOffset++;
    if ( d_CurOffset == dbSize() )
    {
      d_CurOffset = 0;
      if ( d_MultiVector )
      {
        while ( d_CurBlock != d_MultiVector->numberOfDataBlocks () )
        {
          d_CurBlock++;
          if ( dbSize() )
          {
            ConstVectorDataIterator t = d_MultiVector->getIterator ( d_CurBlock );
            d_Vec = t.d_Vec;
            d_Indexer = t.d_Indexer;
            d_Block = d_Vec->getRawDataBlock<double> ( 0 );
            break;
          }
        }
      }
      else
      {
        d_CurBlock++;
      }
    }
    return *this;
  }

  ConstVectorDataIterator   ConstVectorDataIterator::operator++ ( int  )
  {
    if ( d_Indexer )
    {
      advance ( d_Indexer->getIncrement ( d_CurOffset ) , *this );
      return *this;
    }
    AMP_ASSERT ( d_CurOffset < dbSize() );
    d_CurOffset++;
    if ( d_CurOffset == dbSize() )
    {
      d_CurOffset = 0;
      if ( d_MultiVector )
      {
        while ( d_CurBlock != d_MultiVector->numberOfDataBlocks () )
        {
          d_CurBlock++;
          if ( dbSize() )
          {
            ConstVectorDataIterator t = d_MultiVector->getIterator ( d_CurBlock );
            d_Vec = t.d_Vec;
            d_Indexer = t.d_Indexer;
            d_Block = d_Vec->getRawDataBlock<double> ( 0 );
            break;
          }
        }
      }
      else
      {
        d_CurBlock++;
      }
    }
    return *this;
  }

  ConstVectorDataIterator  &ConstVectorDataIterator::operator-- ()
  {
    if ( d_Indexer )
    {
      recede ( d_Indexer->getIncrement ( d_CurOffset ) , *this );
      return *this;
    }
    while ( d_CurOffset == 0 )
    {
      AMP_ASSERT ( d_CurBlock != 0 );
      d_CurBlock--;
      if ( d_MultiVector )
      {
        VectorDataIterator t = d_MultiVector->getIterator ( d_CurBlock );
        d_Vec = t.d_Vec;
        d_Indexer = t.d_Indexer;
        d_Block = d_Vec->getRawDataBlock<double> ( 0 );
      }
      else
      {
        AMP_ASSERT ( d_CurBlock == 0 );
      }
      d_CurOffset = dbSize();
    }
    d_CurOffset--;
    return *this;
  }

  ConstVectorDataIterator   ConstVectorDataIterator::operator-- ( int  )
  {
    if ( d_Indexer )
    {
      recede ( d_Indexer->getIncrement ( d_CurOffset ) , *this );
      return *this;
    }
    while ( d_CurOffset == 0 )
    {
      AMP_ASSERT ( d_CurBlock != 0 );
      d_CurBlock--;
      if ( d_MultiVector )
      {
        VectorDataIterator t = d_MultiVector->getIterator ( d_CurBlock );
        d_Vec = t.d_Vec;
        d_Indexer = t.d_Indexer;
        d_Block = d_Vec->getRawDataBlock<double> ( 0 );
      }
      else
      {
        AMP_ASSERT ( d_CurBlock == 0 );
      }
      d_CurOffset = dbSize();
    }
    d_CurOffset--;
    return *this;
  }

  ConstVectorDataIterator::ConstVectorDataIterator ( const Vector *p , int block , int offset , VectorIndexer::shared_ptr ndx )
       : d_MultiVector ( 0 )
       , d_Vec ( p )
       , d_CurBlock ( block )
       , d_CurOffset ( offset )
       , d_Indexer ( ndx )
  {
    d_Block = d_Vec->getRawDataBlock<double> ( 0 );
  }


  VectorDataIterator   VectorDataIterator::operator-- ( int  )
  {
    if ( d_Indexer )
    {
      recede ( d_Indexer->getIncrement ( d_CurOffset ) , *this );
      return *this;
    }
    while ( d_CurOffset == 0 )
    {
      AMP_ASSERT ( d_CurBlock != 0 );
      d_CurBlock--;
      if ( d_MultiVector )
      {
        VectorDataIterator t = d_MultiVector->getIterator ( d_CurBlock );
        d_Vec = t.d_Vec;
        d_Indexer = t.d_Indexer;
        d_Block = d_Vec->getRawDataBlock<double> ( 0 );
      }
      else
      {
        AMP_ASSERT ( d_CurBlock == 0 );
      }
      d_CurOffset = dbSize();
    }
    d_CurOffset--;
    return *this;
  }
  VectorDataIterator  &VectorDataIterator::operator-- ()
  {
    if ( d_Indexer )
    {
      recede ( d_Indexer->getIncrement ( d_CurOffset ) , *this );
      return *this;
    }
    while ( d_CurOffset == 0 )
    {
      AMP_ASSERT ( d_CurBlock != 0 );
      d_CurBlock--;
      if ( d_MultiVector )
      {
        VectorDataIterator t = d_MultiVector->getIterator ( d_CurBlock );
        d_Vec = t.d_Vec;
        d_Indexer = t.d_Indexer;
        d_Block = d_Vec->getRawDataBlock<double> ( 0 );
      }
      else
      {
        AMP_ASSERT ( d_CurBlock == 0 );
      }
      d_CurOffset = dbSize();
    }
    d_CurOffset--;
    return *this;
  }

  VectorDataIterator   VectorDataIterator::operator++ ( int  )
  {
    if ( d_Indexer )
    {
      advance ( d_Indexer->getIncrement ( d_CurOffset ) , *this );
      return *this;
    }
    AMP_ASSERT ( d_CurOffset < dbSize() );
    d_CurOffset++;
    if ( d_CurOffset == dbSize() )
    {
      d_CurOffset = 0;
      if ( d_MultiVector )
      {
        while ( d_CurBlock != d_MultiVector->numberOfDataBlocks () )
        {
          d_CurBlock++;
          if ( dbSize() )
          {
            VectorDataIterator t = d_MultiVector->getIterator ( d_CurBlock );
            d_Vec = t.d_Vec;
            d_Indexer = t.d_Indexer;
            d_Block = d_Vec->getRawDataBlock<double> ( 0 );
            break;
          }
        }
      }
      else
      {
        d_CurBlock++;
      }
    }
    return *this;
  }

  VectorDataIterator  &VectorDataIterator::operator++ ()
  {
    if ( d_Indexer )
    {
      advance ( d_Indexer->getIncrement ( d_CurOffset ) , *this );
      return *this;
    }
    AMP_ASSERT ( d_CurOffset < dbSize() );
    d_CurOffset++;
    if ( d_CurOffset == dbSize() )
    {
      d_CurOffset = 0;
      if ( d_MultiVector )
      {
        while ( d_CurBlock != d_MultiVector->numberOfDataBlocks () )
        {
          d_CurBlock++;
          if ( dbSize() )
          {
            VectorDataIterator t = d_MultiVector->getIterator ( d_CurBlock );
            d_Vec = t.d_Vec;
            d_Indexer = t.d_Indexer;
            d_Block = d_Vec->getRawDataBlock<double> ( 0 );
            break;
          }
        }
      }
      else
      {
        d_CurBlock++;
      }
    }
    return *this;
  }
  VectorDataIterator  &VectorDataIterator::operator += ( int offset )
  {
    if ( d_Indexer ) // Not constant time
    {
      if ( offset < 0 )
      {
        for ( int i = 0 ; i != offset ; i++ )
        {
         operator --();
        }
      }
      else
      {
        for ( int i = 0 ; i != offset ; i++ )
        {
          operator ++();
        }
      }
    }
    else
    {
      if ( offset > 0 )
      {
        advance ( offset , *this );
      }
      if ( offset < 0 )
      {
        recede ( offset , *this );
      }
    }
    return *this;
  }

  void VectorDataIterator::advance ( int i , VectorDataIterator &d )
  {
    while ( 1 )
    {
      int toGo = d.dbSize() - d.d_CurOffset;
      if ( toGo > i )
      {
        break;
      }
      if ( d.d_CurBlock == d.d_Vec->numberOfDataBlocks() )
      {
        if ( d.d_Indexer )  // Check for a strided iterator
        {
          d.d_CurOffset = 0;
          return;
        }
        AMP_ERROR( "Attempt to advance past the end of a buffer" );
      }
      i -= toGo;
      d.d_CurOffset = 0;
      d.d_CurBlock++;
      if ( d.d_MultiVector )
      {
        if ( d.d_CurBlock == d.d_MultiVector->numberOfDataBlocks() )
        {
          return;
        }
        VectorDataIterator t = d.d_MultiVector->getIterator ( d.d_CurBlock );
        d.d_Vec = t.d_Vec;
        d.d_Indexer = t.d_Indexer;
        d.d_Block = d.d_Vec->getRawDataBlock<double> ( 0 );
      }
      else
      {
        return;
      }
    }
    d.d_CurOffset += i;
  }

  void VectorDataIterator::recede ( int i , VectorDataIterator &d )
  {
    while ( 1 )
    {
      int toGo = d.d_CurOffset;
      if ( toGo > i )
      {
        break;
      }
      if ( d.d_CurBlock == 0 )
      {
        AMP_ERROR( "Attempt to recede past the end of a buffer" );
      }
      i -= toGo;
      d.d_CurBlock--;
      d.d_CurOffset = d.dbSize();
      AMP_ASSERT ( d.d_MultiVector );
      VectorDataIterator t = d.d_MultiVector->getIterator ( d.d_CurBlock );
      d.d_Vec = t.d_Vec;
      d.d_Indexer = t.d_Indexer;
      d.d_Block = d.d_Vec->getRawDataBlock<double> ( 0 );
    }
    d.d_CurOffset -= i;
  }

  int VectorDataIterator::operator - ( const VectorDataIterator &rhs ) const
  {
    int ans;
    if ( !d_MultiVector )
    {
      AMP_ASSERT ( rhs.d_Vec == d_Vec );
    }
    else
    {
      AMP_ASSERT ( rhs.d_MultiVector == d_MultiVector );
    }

    if ( rhs.d_CurBlock < d_CurBlock )
    {
      ans = subtract ( *this , rhs );
    }
    else if ( rhs.d_CurBlock > d_CurBlock )
    {
      ans = -subtract ( rhs , *this );
    }
    else if ( rhs.d_CurOffset < d_CurOffset )
    {
      ans = subtract ( *this , rhs );
    }
    else
    {
      ans = -subtract ( rhs , *this );
    }
    return ans;
  }

  int VectorDataIterator::subtract ( const VectorDataIterator &bigger , const VectorDataIterator &smaller )
  {
    int answer = 0;
    if ( bigger.d_CurBlock != smaller.d_CurBlock )
    {
      for ( size_t i = smaller.d_CurBlock+1 ; i < bigger.d_CurBlock ; i++ )
      {
        AMP_ASSERT ( smaller.d_MultiVector );
        answer += smaller.d_MultiVector->sizeOfDataBlock ( i );
      }
      answer += smaller.dbSize() - smaller.d_CurOffset + bigger.d_CurOffset;
    }
    else
    {
      answer = bigger.d_CurOffset - smaller.d_CurOffset;
    }
    return answer;
  }



  ConstVectorDataIterator  ConstVectorDataIterator::operator + ( int offset )
  {
    ConstVectorDataIterator ans ( *this );
    if ( d_Indexer ) // Not constant time
    {
      if ( offset < 0 )
      {
        for ( int i = 0 ; i != offset ; i++ )
        {
         --ans;
        }
      }
      else
      {
        for ( int i = 0 ; i != offset ; i++ )
        {
          ++ans;
        }
      }
    }
    else
    {
      if ( offset > 0 )
      {
        advance ( offset , ans );
      }
      if ( offset < 0 )
      {
        recede ( offset , ans );
      }
    }
    return ans;
  }

  ConstVectorDataIterator  &ConstVectorDataIterator::operator += ( int offset )
  {
    if ( d_Indexer ) // Not constant time
    {
      if ( offset < 0 )
      {
        for ( int i = 0 ; i != offset ; i++ )
        {
         operator --();
        }
      }
      else
      {
        for ( int i = 0 ; i != offset ; i++ )
        {
          operator ++();
        }
      }
    }
    else
    {
      if ( offset > 0 )
      {
        advance ( offset , *this );
      }
      if ( offset < 0 )
      {
        recede ( offset , *this );
      }
    }
    return *this;
  }

  void ConstVectorDataIterator::advance ( int i , ConstVectorDataIterator &d )
  {
    while ( 1 )
    {
      int toGo = d.dbSize() - d.d_CurOffset;
      if ( toGo > i )
      {
        break;
      }
      if ( d.d_CurBlock == d.d_Vec->numberOfDataBlocks() )
      {
        if ( d.d_Indexer )  // Check for a strided iterator
        {
          d.d_CurOffset = 0;
          return;
        }
        AMP_ERROR( "Attempt to advance past the end of a buffer" );
      }
      i -= toGo;
      d.d_CurOffset = 0;
      d.d_CurBlock++;
      if ( d.d_MultiVector )
      {
        if ( d.d_CurBlock == d.d_MultiVector->numberOfDataBlocks() )
        {
          return;
        }
        VectorDataIterator t = d.d_MultiVector->getIterator ( d.d_CurBlock );
        d.d_Vec = t.d_Vec;
        d.d_Indexer = t.d_Indexer;
        d.d_Block = d.d_Vec->getRawDataBlock<double> ( 0 );
      }
      else
      {
        return;
      }
    }
    d.d_CurOffset += i;
  }

  void ConstVectorDataIterator::recede ( int i , ConstVectorDataIterator &d )
  {
    while ( 1 )
    {
      int toGo = d.d_CurOffset;
      if ( toGo > i )
      {
        break;
      }
      if ( d.d_CurBlock == 0 )
      {
        AMP_ERROR( "Attempt to recede past the end of a buffer" );
      }
      i -= toGo;
      d.d_CurBlock--;
      d.d_CurOffset = d.dbSize();
      AMP_ASSERT ( d.d_MultiVector );
      ConstVectorDataIterator t = d.d_MultiVector->getIterator ( d.d_CurBlock );
      d.d_Vec = t.d_Vec;
      d.d_Indexer = t.d_Indexer;
      d.d_Block = d.d_Vec->getRawDataBlock<double> ( 0 );
    }
    d.d_CurOffset -= i;
  }

  int ConstVectorDataIterator::operator - ( const ConstVectorDataIterator &rhs ) const
  {
    int ans;
    if ( !d_MultiVector )
    {
      AMP_ASSERT ( rhs.d_Vec == d_Vec );
    }
    else
    {
      AMP_ASSERT ( rhs.d_MultiVector == d_MultiVector );
    }
    if ( rhs.d_CurBlock < d_CurBlock )
    {
      ans = subtract ( *this , rhs );
    }
    else if ( rhs.d_CurBlock > d_CurBlock )
    {
      ans = -subtract ( rhs , *this );
    }
    else if ( rhs.d_CurOffset < d_CurOffset )
    {
      ans = subtract ( *this , rhs );
    }
    else
    {
      ans = -subtract ( rhs , *this );
    }
    return ans;
  }

  int ConstVectorDataIterator::subtract ( const ConstVectorDataIterator &bigger , const ConstVectorDataIterator &smaller )
  {
    int answer = 0;
    if ( bigger.d_CurBlock != smaller.d_CurBlock )
    {
      for ( size_t i = smaller.d_CurBlock+1 ; i < bigger.d_CurBlock ; i++ )
      {
        answer += smaller.d_Vec->sizeOfDataBlock ( i );
      }
      answer += smaller.dbSize() - smaller.d_CurOffset + bigger.d_CurOffset;
    }
    else
    {
      answer = bigger.d_CurOffset - smaller.d_CurOffset;
    }
    return answer;
  }

}
}

