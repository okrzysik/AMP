#include "AMP/operators/OperatorFactory.h"
#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/CoupledOperator.h"
#include "AMP/operators/Operator.h"


namespace AMP::Operator {


// Macro to register an operator
#define REGISTER_OPERATOR( NAME )                                         \
    static struct NAME##_INIT {                                           \
        NAME##_INIT()                                                     \
        {                                                                 \
            auto fun = []( std::shared_ptr<OperatorParameters> params ) { \
                return std::make_unique<NAME>( params );                  \
            };                                                            \
            OperatorFactory::registerFactory( #NAME, fun );               \
        }                                                                 \
    } NAME##_INIT


// Create the operator
std::unique_ptr<Operator> OperatorFactory::create( std::shared_ptr<OperatorParameters> parameters )
{
    AMP_ASSERT( parameters != nullptr );
    auto inputDatabase = parameters->d_db;
    AMP_ASSERT( inputDatabase );
    auto objectName = inputDatabase->getString( "name" );
    return create( objectName, parameters );
}


// Register operators
REGISTER_OPERATOR( CoupledOperator );
REGISTER_OPERATOR( ColumnOperator );


} // namespace AMP::Operator
