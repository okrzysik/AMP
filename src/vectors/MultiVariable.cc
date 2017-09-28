#include "MultiVariable.h"
#include "utils/Utilities.h"
#include <algorithm>
#include <map>
#include <utility>


namespace AMP {
namespace LinearAlgebra {


class MVSortByName
{
private:
    std::map<std::string, int> new_order;

public:
    explicit MVSortByName( const std::vector<std::string> &in )
    {
        auto cur = in.begin();
        int i    = 0;
        while ( cur != in.end() ) {
            new_order[*cur] = i++;
            ++cur;
        }
    }

    bool operator()( const Variable::shared_ptr left, const Variable::shared_ptr right )
    {
        std::string lname = left->getName();
        std::string rname = right->getName();
        return new_order[lname] < new_order[rname];
    }
};


/****************************************************************
 * Constructors/Destructors                                      *
 ****************************************************************/
MultiVariable::MultiVariable( const std::string &name, const std::vector<Variable::shared_ptr>& vars )
    : Variable( name ), d_vVariables( vars )
{
}
MultiVariable::~MultiVariable() = default;


/****************************************************************
 * Get/set a variable                                            *
 ****************************************************************/
Variable::shared_ptr MultiVariable::getVariable( size_t which )
{
    AMP_ASSERT( which < d_vVariables.size() );
    return d_vVariables[which];
}
Variable::const_shared_ptr MultiVariable::getVariable( size_t which ) const
{
    AMP_ASSERT( which < d_vVariables.size() );
    return d_vVariables[which];
}
void MultiVariable::setVariable( size_t i, Variable::shared_ptr &p )
{
    AMP_ASSERT( i < d_vVariables.size() );
    d_vVariables[i] = p;
}
void MultiVariable::add( Variable::shared_ptr newVar )
{
    AMP::shared_ptr<MultiVariable> multivariable =
        AMP::dynamic_pointer_cast<MultiVariable>( newVar );
    if ( multivariable.get() != nullptr ) {
        auto curVar = multivariable->beginVariable();
        while ( curVar != multivariable->endVariable() ) {
            add( *curVar );
            ++curVar;
        }
    } else {
        d_vVariables.push_back( newVar );
    }
}


/****************************************************************
 * Misc                                                          *
 ****************************************************************/
size_t MultiVariable::numVariables() const { return d_vVariables.size(); }

void MultiVariable::sortVariablesByName( const std::vector<std::string> &order )
{
    MVSortByName sorter( order );
    std::sort( beginVariable(), endVariable(), sorter );
}


bool MultiVariable::operator==( const Variable &rhs ) const
{
    const auto *multivariable = dynamic_cast<const MultiVariable *>( &rhs );
    if ( multivariable == nullptr ) {
        // We are comparing a multi variable to another variable
        // The two variables match if the variable equals all sub-variable and
        // the names match
        if ( rhs.getName() != this->getName() )
            return false;
        for ( size_t i = 0; i != d_vVariables.size(); i++ ) {
            if ( *d_vVariables[i] != rhs )
                return false;
        }
    } else {
        // We are dealing with two multivariables, check that the internal variables match
        for ( size_t i = 0; i != d_vVariables.size(); i++ ) {
            if ( i == multivariable->d_vVariables.size() )
                return false;
            if ( ( *d_vVariables[i] ) != ( *( multivariable->d_vVariables[i] ) ) )
                return false;
        }
    }
    return true;
}


Variable::shared_ptr MultiVariable::cloneVariable( const std::string &name ) const
{
    AMP::shared_ptr<MultiVariable> retVal( new MultiVariable( name ) );
    retVal->d_vVariables.resize( d_vVariables.size() );
    std::copy( d_vVariables.begin(), d_vVariables.end(), retVal->d_vVariables.begin() );
    return retVal;
}


void MultiVariable::removeDuplicateVariables()
{
    // First remove any NULL pointers
    auto iterator = d_vVariables.begin();
    while ( iterator != d_vVariables.end() ) {
        if ( iterator->get() == nullptr )
            iterator = d_vVariables.erase( iterator );
        else
            ++iterator;
    }
    // Next remove any duplicate entries
    // Note: while it would be faster to sort, then remove duplicate entires,
    // this requires the < operator to be overloaded for the base class
    std::vector<Variable::shared_ptr> unique_list;
    unique_list.reserve( d_vVariables.size() );
    for ( auto &elem : d_vVariables ) {
        bool found = false;
        for ( auto &unique_list_j : unique_list ) {
            if ( elem->operator==( *( unique_list_j ) ) )
                found = true;
        }
        if ( !found )
            unique_list.push_back( elem );
    }
    d_vVariables = unique_list;
}
} // namespace LinearAlgebra
} // namespace AMP
