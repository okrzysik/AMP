
#include "Variable.h"

namespace AMP {
namespace LinearAlgebra {

  Variable::Variable ( const std::string &n ) : d_VariableName ( n ) 
  {
  }

  Variable::~Variable () 
  {
  }

  Variable::shared_ptr  Variable::cloneVariable ( const std::string &name ) const
  {
    return Variable::shared_ptr ( new Variable ( name ) );
  }

  size_t  Variable::DOFsPerObject () const
  {
    return 0;
  }

  size_t  Variable::variableID () const
  {
    AMP_ERROR( "Cannot use base class Variable for memory assignment" );
    return 0;
  }

  const std::string  &Variable::getName() const 
  { 
    return d_VariableName; 
  }

  bool Variable::operator == ( const Variable & rhs ) const
  {
    return d_VariableName == rhs.d_VariableName;
  }

  bool Variable::operator != ( const Variable & rhs ) const
  {
    return ( !( (*this) == rhs ) );
  }

  void Variable::setName ( const std::string &s ) 
  { 
    d_VariableName = s; 
  }

  Variable::shared_ptr  Variable::cloneVariable () const 
  { 
    return cloneVariable ( d_VariableName ); 
  }

  bool  Variable::isSameTypeAs ( Variable::shared_ptr rhs ) const
  {
    return (DOFsPerObject()-variableID()) == (rhs->DOFsPerObject()-rhs->variableID());
  }

}
}

