
#include "Variable.h"

#include <utility>


namespace AMP {
namespace LinearAlgebra {


Variable::Variable( const std::string &name ) : d_VariableName( name ) {}


Variable::~Variable() = default;


Variable::shared_ptr Variable::cloneVariable( const std::string &name ) const
{
    return AMP::make_shared<Variable>( name );
}


const std::string &Variable::getName() const { return d_VariableName; }


bool Variable::operator==( const Variable &rhs ) const
{
    return d_VariableName == rhs.d_VariableName;
}


bool Variable::operator!=( const Variable &rhs ) const { return ( !( ( *this ) == rhs ) ); }

void Variable::setUnits( const std::string &t ) { d_Units = t; }


const std::string &Variable::getUnits() const { return d_Units; }
} // namespace LinearAlgebra
} // namespace AMP
