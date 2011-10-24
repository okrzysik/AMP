
#include "NullVariable.h"

namespace AMP {
namespace LinearAlgebra {

  NullVariable::NullVariable ( const std::string &a ) : Variable ( a ) 
  {
  }

  size_t  NullVariable::variableID () const 
  { 
    return static_cast<unsigned int>(-1); 
  }

  size_t  NullVariable::DOFsPerObject () const 
  { 
    return 0; 
  }

  Variable::shared_ptr  NullVariable::cloneVariable ( const std::string &a ) const
  { 
    return Variable::shared_ptr ( new NullVariable ( a ) ); 
  }

}
}

