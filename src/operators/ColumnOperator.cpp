#include "AMP/operators/ColumnOperator.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/MultiVariable.h"
#include "ProfilerApp.h"

namespace AMP {
namespace Operator {


/********************************************************
 * Constructors                                          *
 ********************************************************/
ColumnOperator::ColumnOperator( const std::shared_ptr<OperatorParameters> &params )
    : Operator( params )
{
}


/********************************************************
 * Apply/residual                                        *
 ********************************************************/
void ColumnOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                            AMP::LinearAlgebra::Vector::shared_ptr f )
{
    for ( auto &elem : d_operators ) {
        elem->apply( u, f );
    }
}
void ColumnOperator::residual( AMP::LinearAlgebra::Vector::const_shared_ptr f,
                               AMP::LinearAlgebra::Vector::const_shared_ptr u,
                               AMP::LinearAlgebra::Vector::shared_ptr r )
{
    for ( auto &elem : d_operators ) {
        AMP_INSIST( ( elem.get() != nullptr ), "ColumnOperator::operator component is NULL" );
        elem->residual( f, u, r );
    }
}


/********************************************************
 * getParameters                                         *
 ********************************************************/
std::shared_ptr<OperatorParameters>
ColumnOperator::getParameters( const std::string &type,
                               AMP::LinearAlgebra::Vector::const_shared_ptr u,
                               std::shared_ptr<OperatorParameters> params )
{
    std::shared_ptr<AMP::Database> db;
    auto opParameters = std::make_shared<ColumnOperatorParameters>( db );
    opParameters->d_OperatorParameters.resize( d_operators.size() );

    for ( unsigned int i = 0; i < d_operators.size(); i++ ) {
        opParameters->d_OperatorParameters[i] = d_operators[i]->getParameters( type, u, params );
    }
    return opParameters;
}


/********************************************************
 * reset                                                  *
 ********************************************************/
void ColumnOperator::reset( const std::shared_ptr<OperatorParameters> &params )
{
    auto columnParameters = std::dynamic_pointer_cast<ColumnOperatorParameters>( params );
    AMP_INSIST( ( columnParameters.get() != nullptr ),
                "ColumnOperator::reset parameter object is NULL" );
    AMP_INSIST( ( ( ( columnParameters->d_OperatorParameters ).size() ) == ( d_operators.size() ) ),
                " std::vector sizes do not match! " );
    for ( size_t i = 0; i < d_operators.size(); i++ ) {
        d_operators[i]->reset( ( columnParameters->d_OperatorParameters )[i] );
    }
}


/********************************************************
 * Add an operator                                       *
 ********************************************************/
void ColumnOperator::append( std::shared_ptr<Operator> op )
{
    AMP_INSIST( ( op.get() != nullptr ),
                "AMP::ColumnOperator::appendRow input argument is a NULL operator" );
    d_operators.push_back( op );
}


/********************************************************
 * Get the variables                                     *
 ********************************************************/
AMP::LinearAlgebra::Variable::shared_ptr ColumnOperator::getInputVariable()
{
    auto retVariable = std::make_shared<AMP::LinearAlgebra::MultiVariable>( "ColumnVariable" );
    for ( auto &elem : d_operators ) {
        auto opVar = elem->getInputVariable();
        if ( opVar.get() != nullptr ) {
            retVariable->add( opVar );
        }
    }
    retVariable->removeDuplicateVariables();
    return retVariable;
}
AMP::LinearAlgebra::Variable::shared_ptr ColumnOperator::getOutputVariable()
{
    auto retVariable = std::make_shared<AMP::LinearAlgebra::MultiVariable>( "ColumnVariable" );
    for ( auto &elem : d_operators ) {
        AMP::LinearAlgebra::Variable::shared_ptr opVar = elem->getOutputVariable();
        if ( opVar.get() != nullptr ) {
            retVariable->add( opVar );
        }
    }
    retVariable->removeDuplicateVariables();
    return retVariable;
}


/********************************************************
 * Check the input                                       *
 ********************************************************/
bool ColumnOperator::isValidInput( std::shared_ptr<AMP::LinearAlgebra::Vector> &u )
{
    bool bRetVal = true;
    for ( auto &elem : d_operators ) {
        bRetVal = bRetVal && elem->isValidInput( u );
    }
    return bRetVal;
}


/********************************************************
 * find                                                  *
 ********************************************************/
std::vector<std::shared_ptr<Operator>> ColumnOperator::find( const std::string &name )
{
    std::vector<std::shared_ptr<Operator>> ops;
    for ( auto op : d_operators ) {
        if ( op->type() == name )
            ops.push_back( op );
    };
    return ops;
}


} // namespace Operator
} // namespace AMP
