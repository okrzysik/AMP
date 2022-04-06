#include "Variable.h"
#include "VectorSelector.h"


namespace AMP::LinearAlgebra {


Variable::Variable( const std::string &name ) : d_VariableName( name ) {}


Variable::~Variable() = default;


std::shared_ptr<Variable> Variable::cloneVariable( const std::string &name ) const
{
    return std::make_shared<Variable>( name );
}


bool Variable::operator==( const Variable &rhs ) const
{
    return d_VariableName == rhs.d_VariableName;
}


bool Variable::operator!=( const Variable &rhs ) const { return !( *this == rhs ); }


void Variable::setUnits( const Units &t ) { d_Units = t; }


const Units &Variable::getUnits() const { return d_Units; }


std::shared_ptr<VectorSelector> Variable::createVectorSelector() const
{
    return std::make_shared<VS_ByVariableName>( d_VariableName );
}


} // namespace AMP::LinearAlgebra
