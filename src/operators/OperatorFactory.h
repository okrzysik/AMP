#ifndef included_AMP_Operator_OperatorFactory
#define included_AMP_Operator_OperatorFactory

#include "AMP/utils/FactoryStrategy.hpp"

#include <memory>


namespace AMP::Operator {

class Operator;
class OperatorParameters;


//! Operator factory class
class OperatorFactory : public FactoryStrategy<Operator, std::shared_ptr<OperatorParameters>>
{
public:
    static std::unique_ptr<Operator> create( std::shared_ptr<OperatorParameters> parameters );
};


//! Register operators known by AMP
void registerOperatorFactories();


} // namespace AMP::Operator

#endif
