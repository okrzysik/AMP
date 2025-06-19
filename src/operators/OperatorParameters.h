#ifndef included_AMP_OperatorParameters
#define included_AMP_OperatorParameters

#include "AMP/utils/ParameterBase.h"
#include "AMP/utils/Utilities.h"

#include <memory>


namespace AMP::Mesh {
class Mesh;
}


namespace AMP::Operator {

class Operator;

/**\class OperatorParameters
 *
 * OperatorParameters encapsulates parameters used to initialize or reset
 * operators. It is an abstract base class.
 */

class OperatorParameters : public ParameterBase
{
public:
    /**
     * Construct and initialize a parameter list according to input
     * data.  Guess what the required and optional keywords are.
     */
    explicit OperatorParameters( std::shared_ptr<AMP::Database> db,
                                 std::shared_ptr<AMP::Mesh::Mesh> mesh = nullptr );

    // Get the memory location from a string
    static AMP::Utilities::MemoryType memoryLocationFromString( const std::string &name );

    //! Destructor
    virtual ~OperatorParameters() {}

    //! Optional mesh for the operator
    std::shared_ptr<AMP::Mesh::Mesh> d_Mesh;

    //! Allow for the case that a fully constructed operator is returned
    std::shared_ptr<AMP::Operator::Operator> d_pOperator;

    /**
     * Location (host/managed/device) where internally created
     * vectors and matrices should live. Host memory should always
     * work. More specialized parameter classes can overwrite this
     * if supported.
     */
    AMP::Utilities::MemoryType d_memory_location = AMP::Utilities::MemoryType::host;
};


} // namespace AMP::Operator

#endif
