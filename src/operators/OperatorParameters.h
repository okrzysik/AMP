
#ifndef included_AMP_OperatorParameters
#define included_AMP_OperatorParameters

#include <memory>

#include "AMP/mesh/Mesh.h"
#include "AMP/utils/Database.h"
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
    explicit OperatorParameters( std::shared_ptr<AMP::Database> db ) : d_db( db ) {}

    /**
     * Destructor.
     */
    virtual ~OperatorParameters() {}

    /**
     *  Database object which needs to be initialized specific to the solver.
     *  Documentation for parameters required by each solver can be found in the
     *  documentation for the solver.
     */
    std::shared_ptr<AMP::Database> d_db;

    AMP::Mesh::Mesh::shared_ptr d_Mesh;

protected:
private:
};
} // namespace AMP::Operator

#endif
