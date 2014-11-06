#ifndef included_AMP_AndersonStatusTest
#define included_AMP_AndersonStatusTest

// AMP includes
#include "utils/Database.h"

// Trilinos includes
#include "NOX_StatusTest_Generic.H"
#include "NOX_Thyra.H"
#include "NOX_Solver_Generic.H"


namespace AMP {
namespace Solver {


/**
  * The AndersonStatusTest implements a AMP-specific stopping criteria
  * for use with the NOX Anderson acceleration solver.
  */
class AndersonStatusTest : public NOX::StatusTest::Generic {
public:

    /**
     * main constructor
     @param [in] db The database object containing teh following fields:

     1. type: vector<string>, name: AndersonConvergenceVariables, default value: "",
     acceptable values (any variable name on solution vector)

     2. type: vector<double>, name: AndersonConvergenceTolerances, default value: "",
     acceptable values (>0.0)

    */
     AndersonStatusTest(AMP::shared_ptr<AMP::Database> db);

     /**
      * Default destructor.
      */
    virtual ~AndersonStatusTest();

    /**
    * Test the stopping criterion
    * @param solver Instance of solver for which convergence is being evaluated.
    *        checkType Type of convergence check
    */
    NOX::StatusTest::StatusType checkStatus( const NOX::Solver::Generic &solver, NOX::StatusTest::CheckType checkType);

    /**
    * Return the result of the most recent checkStatus call.
    */
    NOX::StatusTest::StatusType getStatus() const;

    /**
    * Output formatted description of stopping test to output stream.
    */
    std::ostream & print (std::ostream &stream, int indent=0) const;

protected:

    // Current status
    NOX::StatusTest::StatusType d_status;

    // List of variables to evaluate convergence
    std::vector<std::string> d_variableNames;
    std::vector<double>      d_tolerances;
    std::vector<double>      d_relativeResiduals;
};

}
}

#endif
