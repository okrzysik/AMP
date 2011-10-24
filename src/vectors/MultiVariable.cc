
#include "utils/Utilities.h"
#include "MultiVariable.h"
#include <map>

namespace AMP {
namespace LinearAlgebra {

  MultiVariable::MultiVariable ( const std::string &name ) : Variable ( name ) {}

  MultiVariable::~MultiVariable () {}

  Variable::shared_ptr  MultiVariable::getVariable ( size_t which )
  {
    AMP_ASSERT ( which < d_vVariables.size() );
    return d_vVariables[which];
  }

  size_t  MultiVariable::variableID ()  const
  { 
    AMP_ERROR( "I was told this situation would never happen. ``Composing systems will not involve composing matrices,'' said Bobby.  -- Bill" ); 
    return 0;
  }

  size_t  MultiVariable::numVariables ()
  {
    return d_vVariables.size();
  }

  void MultiVariable::setVariable ( size_t i , Variable::shared_ptr & p ) 
  { 
    AMP_ASSERT ( i < d_vVariables.size() );
    d_vVariables[i] = p; 
  }

  class MVSortByName
  {
    private:
      std::map<std::string,int> new_order;

    public:
      MVSortByName ( const std::vector<std::string> &in ) 
      {
        std::vector<std::string>::const_iterator cur = in.begin();
        int i = 0;
        while ( cur != in.end() )
        {
          new_order[*cur] = i++;
          cur++;
        }
      }

      bool operator () ( const Variable::shared_ptr left , const Variable::shared_ptr right ) 
      {
        std::string lname = left->getName();
        std::string rname = right->getName();
        return new_order[lname] < new_order[rname];
      }
  };

  void   MultiVariable::sortVariablesByName ( const std::vector<std::string> &order )
  {
    MVSortByName sorter ( order );
    std::sort ( beginVariable() , endVariable() , sorter );
  }

  void   MultiVariable::add ( Variable::shared_ptr newVar ) 
  {
    if ( (newVar.get() != NULL) && (newVar->isA<MultiVariable> ()) )
    {
      iterator curVar = newVar->castTo<MultiVariable>().beginVariable();
      while ( curVar != newVar->castTo<MultiVariable>().endVariable() )
      {
        add ( *curVar );
        curVar++;
      }
    }
    else
    {
      d_vVariables.push_back ( newVar ); 
    }
  }

  bool   MultiVariable::operator == ( const Variable &rhs ) const
  { 
    if ( !rhs.isA<MultiVariable>() ) return false;
    for ( size_t i = 0 ; i != d_vVariables.size() ; i++ )
    {
      if ( i == rhs.castTo<MultiVariable>().d_vVariables.size() )
      {
        return false;
      }
      if ( d_vVariables[i] != rhs.castTo<MultiVariable>().d_vVariables[i] )
      {
        return false;
      }
    }
    return d_VariableName == rhs.castTo<MultiVariable>().d_VariableName;
  }

  Variable::shared_ptr  MultiVariable::cloneVariable ( const std::string &name ) const
  {
    Variable::shared_ptr retVal ( new MultiVariable ( name ) );
    MultiVariable  &ans = retVal->castTo<MultiVariable> ();
    ans.d_vVariables.resize ( d_vVariables.size() );
    std::copy ( d_vVariables.begin() , d_vVariables.end() , ans.d_vVariables.begin() );
    return retVal;
  }
  
}
}

