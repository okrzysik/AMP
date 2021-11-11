#include "Variable.h"

#include <utility>


namespace AMP::LinearAlgebra {


Variable::Variable( const std::string &name ) : d_VariableName( name ) {}


Variable::~Variable() = default;


std::shared_ptr<Variable> Variable::cloneVariable( const std::string &name ) const
{
    return std::make_shared<Variable>( name );
}


const std::string &Variable::getName() const { return d_VariableName; }


bool Variable::operator==( const Variable &rhs ) const
{
    return d_VariableName == rhs.d_VariableName;
}


bool Variable::operator!=( const Variable &rhs ) const { return !( *this == rhs ); }

void Variable::setUnits( const Units &t ) { d_Units = t; }


const Units &Variable::getUnits() const { return d_Units; }

} // namespace AMP::LinearAlgebra
