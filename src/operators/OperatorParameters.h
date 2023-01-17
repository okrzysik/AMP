
#ifndef included_AMP_OperatorParameters
#define included_AMP_OperatorParameters

#include <memory>

#include "AMP/mesh/Mesh.h"
#include "AMP/utils/ParameterBase.h"


namespace AMP::Operator {

/**\class OperatorParameters
 *
 * OperatorParameters encapsulates parameters used to initialize or reset
 * operators. It is an abstract base class.
 */

class OperatorParameters : public ParameterBase
{
public:
    typedef std::shared_ptr<AMP::Operator::OperatorParameters> shared_ptr;

    /**
     * Construct and initialize a parameter list according to input
     * data.  Guess what the required and optional keywords are.
     */
    explicit OperatorParameters( std::shared_ptr<AMP::Database> db ) : ParameterBase( db ) {}

    /**
     * Destructor.
     */
    virtual ~OperatorParameters() {}

    AMP::Mesh::Mesh::shared_ptr d_Mesh;

protected:
private:
};
} // namespace AMP::Operator

#endif
