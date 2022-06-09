#include "AMP/operators/OperatorFactory.h"
#include "AMP/operators/Operator.h"


namespace AMP::Operator {


// Create the operator
std::unique_ptr<Operator> OperatorFactory::create( std::shared_ptr<OperatorParameters> parameters )
{
    AMP_ASSERT( parameters != nullptr );
    auto inputDatabase = parameters->d_db;
    AMP_ASSERT( inputDatabase );
    auto objectName = inputDatabase->getString( "name" );
    return FactoryStrategy<Operator, std::shared_ptr<OperatorParameters>>::create( objectName,
                                                                                   parameters );
}


// Register operators
void registerOperatorFactories( void ) {}


} // namespace AMP::Operator
