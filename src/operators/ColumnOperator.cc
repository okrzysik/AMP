#include "AMP/operators/ColumnOperator.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/MultiVariable.h"
#include "ProfilerApp.h"

namespace AMP {
namespace Operator {

void ColumnOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                            AMP::LinearAlgebra::Vector::shared_ptr f )
{
    for ( auto &elem : d_Operators ) {
        elem->apply( u, f );
    }
}

void ColumnOperator::residual( AMP::LinearAlgebra::Vector::const_shared_ptr f,
                               AMP::LinearAlgebra::Vector::const_shared_ptr u,
                               AMP::LinearAlgebra::Vector::shared_ptr r )
{
    for ( auto &elem : d_Operators ) {
        AMP_INSIST( ( elem.get() != nullptr ), "ColumnOperator::operator component is NULL" );
        elem->residual( f, u, r );
    }
}

std::shared_ptr<OperatorParameters>
ColumnOperator::getParameters( const std::string &type,
                               AMP::LinearAlgebra::Vector::const_shared_ptr u,
                               std::shared_ptr<OperatorParameters> params )
{
    std::shared_ptr<AMP::Database> db;
    std::shared_ptr<ColumnOperatorParameters> opParameters( new ColumnOperatorParameters( db ) );

    ( opParameters->d_OperatorParameters ).resize( d_Operators.size() );

    for ( unsigned int i = 0; i < d_Operators.size(); i++ ) {
        ( opParameters->d_OperatorParameters )[i] =
            ( d_Operators[i]->getParameters( type, u, params ) );
    }
    return opParameters;
}

void ColumnOperator::reset( const std::shared_ptr<OperatorParameters> &params )
{
    std::shared_ptr<ColumnOperatorParameters> columnParameters =
        std::dynamic_pointer_cast<ColumnOperatorParameters>( params );

    AMP_INSIST( ( columnParameters.get() != nullptr ),
                "ColumnOperator::reset parameter object is NULL" );

    AMP_INSIST( ( ( ( columnParameters->d_OperatorParameters ).size() ) == ( d_Operators.size() ) ),
                " std::vector sizes do not match! " );

    for ( unsigned int i = 0; i < d_Operators.size(); i++ ) {
        d_Operators[i]->reset( ( columnParameters->d_OperatorParameters )[i] );
    }
}

void ColumnOperator::append( std::shared_ptr<Operator> op )
{
    AMP_INSIST( ( op.get() != nullptr ),
                "AMP::ColumnOperator::appendRow input argument is a NULL operator" );

    d_Operators.push_back( op );
}

AMP::LinearAlgebra::Variable::shared_ptr ColumnOperator::getInputVariable()
{
    std::shared_ptr<AMP::LinearAlgebra::MultiVariable> retVariable(
        new AMP::LinearAlgebra::MultiVariable( "ColumnVariable" ) );

    for ( auto &elem : d_Operators ) {
        AMP::LinearAlgebra::Variable::shared_ptr opVar = elem->getInputVariable();
        if ( opVar.get() != nullptr ) {
            retVariable->add( opVar );
        }
    }
    retVariable->removeDuplicateVariables();

    return retVariable;
}

AMP::LinearAlgebra::Variable::shared_ptr ColumnOperator::getOutputVariable()
{
    std::shared_ptr<AMP::LinearAlgebra::MultiVariable> retVariable(
        new AMP::LinearAlgebra::MultiVariable( "ColumnVariable" ) );

    for ( auto &elem : d_Operators ) {
        AMP::LinearAlgebra::Variable::shared_ptr opVar = elem->getOutputVariable();
        if ( opVar.get() != nullptr ) {
            retVariable->add( opVar );
        }
    }
    retVariable->removeDuplicateVariables();

    return retVariable;
}

bool ColumnOperator::isValidInput( std::shared_ptr<AMP::LinearAlgebra::Vector> &u )
{
    bool bRetVal = true;

    for ( auto &elem : d_Operators ) {
        bRetVal = bRetVal && elem->isValidInput( u );
    }

    return bRetVal;
}
} // namespace Operator
} // namespace AMP
