#ifndef included_AMP_OperatorFactory_H_
#define included_AMP_OperatorFactory_H_

/* AMP files */

#include "AMP/utils/FactoryStrategy.hpp"

namespace AMP::Operator {

class Operator;
class OperatorParameters;

using OperatorFactory = FactoryStrategy<Operator, OperatorParameters>;

// free function to preregister operators known by AMP
// none at present
void registerOperatorFactories();

} // namespace AMP::Operator

#endif
