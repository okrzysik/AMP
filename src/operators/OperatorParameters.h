
#ifndef included_AMP_OperatorParameters
#define included_AMP_OperatorParameters

#include "AMP/utils/shared_ptr.h"

#include "AMP/ampmesh/Mesh.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/ParameterBase.h"


namespace AMP {
namespace Operator {

/**\class OperatorParameters
 *
 * OperatorParameters encapsulates parameters used to initialize or reset
 * operators. It is an abstract base class.
 */

class OperatorParameters : public ParameterBase
{
public:
    typedef AMP::shared_ptr<AMP::Operator::OperatorParameters> shared_ptr;

    /**
     * Construct and initialize a parameter list according to input
     * data.  Guess what the required and optional keywords are.
     */
    explicit OperatorParameters( const AMP::shared_ptr<AMP::Database> &db ) : d_db( db ) {}

    /**
     * Destructor.
     */
    virtual ~OperatorParameters() {}

    /**
     *  Database object which needs to be initialized specific to the solver.
     *  Documentation for parameters required by each solver can be found in the
     *  documentation for the solver.
     */
    AMP::shared_ptr<AMP::Database> d_db;

    AMP::Mesh::Mesh::shared_ptr d_Mesh;

protected:
private:
};
} // namespace Operator
} // namespace AMP

#endif
