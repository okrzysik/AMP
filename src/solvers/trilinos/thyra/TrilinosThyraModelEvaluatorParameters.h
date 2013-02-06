#ifndef included_AMP_TrilinosThyraModelEvaluatorParameters
#define included_AMP_TrilinosThyraModelEvaluatorParameters


#include "discretization/DOF_Manager.h"


namespace AMP {
namespace Solver {


/**
  * The TrilinosThyraModelEvaluator is a wrapper for a Thyra ModelEvaluator to 
  * wrap AMP::Operators for use with Trilinos NOX solvers.
  */
class TrilinosThyraModelEvaluatorParameters
{
public:
    
    AMP::Discretization::DOFManager::shared_ptr d_dofs;

};


}
}

#endif

