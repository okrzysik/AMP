#include "utils/Utilities.h"
#include "ManagedVector.h"
#include <stdexcept>



namespace AMP {
namespace LinearAlgebra {


  Vector::shared_ptr  ManagedVector::subsetVectorForVariable ( const Variable::shared_ptr &name )
  {
    Vector::shared_ptr  retVal;
    if ( !d_vBuffer )
    {
      retVal = d_Engine->castTo<Vector>().subsetVectorForVariable ( name );
    }
    if ( !retVal )
    {
      retVal = Vector::subsetVectorForVariable ( name );
    }
    return retVal;
  }

  bool ManagedVector::isAnAliasOf ( Vector &rhs )
  {
    bool retVal = false;
    if ( rhs.isA<ManagedVector>() )
    {
      ManagedVector &other = rhs.castTo<ManagedVector> ();
      if ( d_vBuffer && ( other.d_vBuffer == d_vBuffer ) )
      {
        retVal = true;
      }
    }
    return retVal;
  }

  ManagedVector::ManagedVector ( shared_ptr  alias )
    : Vector ( boost::dynamic_pointer_cast<VectorParameters> ( alias->castTo<ManagedVector>().getParameters() ) ) ,
      d_vBuffer ( alias->castTo<ManagedVector>().d_vBuffer ) ,
      d_Engine ( alias->castTo<ManagedVector>().d_Engine->cloneEngine ( d_vBuffer ) )
  {
    d_Engine = alias->castTo<ManagedVector>().d_Engine;
    setVariable ( alias->getVariable() );
    d_pParameters = alias->castTo<ManagedVector>().d_pParameters;
    aliasGhostBuffer ( alias );
  }


  void ManagedVector::copyVector ( const Vector &other )
  {
    if ( other.getLocalSize() != getLocalSize() )
    {  // Another error condition
      AMP_ERROR( "Destination vector and source vector not the same size" );
    }
    fireDataChange();
    VectorDataIterator  cur1 = begin();
    VectorDataIterator  end1 = end();
    ConstVectorDataIterator  cur2 = other.begin();
    while ( cur1 != end1 )
    {
      *cur1 = *cur2;
      ++cur1;
      ++cur2;
    }
    copyGhostValues ( other );
  }


  void ManagedVector::swapVectors ( Vector &other )
  {
    ManagedVector &in = other.castTo<ManagedVector> ();
    d_vBuffer.swap ( in.d_vBuffer );
    std::swap ( d_pParameters , in.d_pParameters );
    d_Engine->swapEngines ( in.d_Engine );
  }


  void ManagedVector::aliasVector ( Vector &other )
  {
    ManagedVector &in = other.castTo<ManagedVector> ();
    d_pParameters = in.d_pParameters;
    d_vBuffer = in.d_vBuffer;
  }

  void ManagedVector::getLocalValuesByGlobalID ( int numVals , int *ndx , double *vals ) const
  {
    INCREMENT_COUNT("Virtual");
    if ( d_vBuffer )
    {
      for ( int i = 0 ; i != numVals ; i++ )
        vals[i] = (*d_vBuffer)[ndx[i] - d_CommList->getStartGID() ];
    }
    else
    {
      d_Engine->getLocalValuesByGlobalID ( numVals , ndx , vals );
    }
  }

}
}
