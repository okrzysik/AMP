
#include "utils/Utilities.h"
#include "MultiVariable.h"
#include <map>

namespace AMP {
namespace LinearAlgebra {


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


MultiVariable::MultiVariable ( const std::string &name ) : Variable ( name ) {}


MultiVariable::~MultiVariable () {}


Variable::shared_ptr  MultiVariable::getVariable ( size_t which )
{
    AMP_ASSERT ( which < d_vVariables.size() );
    return d_vVariables[which];
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


void   MultiVariable::sortVariablesByName ( const std::vector<std::string> &order )
{
    MVSortByName sorter ( order );
    std::sort ( beginVariable() , endVariable() , sorter );
}


void   MultiVariable::add ( Variable::shared_ptr newVar ) 
{
    boost::shared_ptr<MultiVariable> multivariable = boost::dynamic_pointer_cast<MultiVariable>(newVar);
    if ( multivariable.get() != NULL ) {
        iterator curVar = multivariable->beginVariable();
        while ( curVar != multivariable->endVariable() ) {
            add ( *curVar );
            curVar++;
        }
    } else {
        d_vVariables.push_back ( newVar ); 
    }
}


bool   MultiVariable::operator == ( const Variable &rhs ) const
{ 
    const MultiVariable *multivariable = dynamic_cast<const MultiVariable*>( &rhs );
    if ( multivariable==NULL ) 
        return false;
    for (size_t i=0; i!=d_vVariables.size(); i++) {
        if ( i == multivariable->d_vVariables.size() )
            return false;
        if ( d_vVariables[i] != multivariable->d_vVariables[i] )
            return false;
    }
    return d_VariableName == multivariable->d_VariableName;
}


Variable::shared_ptr  MultiVariable::cloneVariable ( const std::string &name ) const
{
    boost::shared_ptr<MultiVariable> retVal( new MultiVariable ( name ) );
    retVal->d_vVariables.resize ( d_vVariables.size() );
    std::copy ( d_vVariables.begin() , d_vVariables.end() , retVal->d_vVariables.begin() );
    return retVal;
}

  
}
}

