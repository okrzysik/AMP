#include "DualVariable.h"

namespace AMP {
namespace LinearAlgebra {

  size_t  DualVariable::variableID () const
  { 
    AMP_ERROR( "I was told this situation would never happen. ``Composing systems will not involve composing matrices,'' said Bobby.  -- Bill" ); 
    return 0;
  }

  bool   DualVariable::operator == ( const Variable &rhs ) const
  { 
    if ( !rhs.isA<DualVariable>() ) return false;
    return ( d_Var1 == rhs.castTo<DualVariable>().d_Var1 ) && 
           ( d_Var2 == rhs.castTo<DualVariable>().d_Var2 ) &&
           ( d_VariableName == rhs.castTo<DualVariable>().d_VariableName ); 
  }

  size_t   DualVariable::DOFsPerObject () const 
  {
    AMP_ERROR( "Not implemented for dual variables" );
    return 0;
  }

}
}

