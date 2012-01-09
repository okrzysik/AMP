#ifndef included_SolverStrategy
#include "SolverStrategy.h"
#endif
#include "utils/Utilities.h"
#include "utils/Utilities.h"


namespace AMP {
  namespace Solver {

    int SolverStrategy::d_iInstanceId=0;

    SolverStrategy::SolverStrategy()
    {
      d_iNumberIterations      = -1;
      d_dResidualNorm          = -1;
      d_bUseZeroInitialGuess = true;
      d_iDebugPrintInfoLevel = 0;
    }

    SolverStrategy::SolverStrategy(boost::shared_ptr<SolverStrategyParameters> parameters)
    {
      AMP_INSIST(parameters.get()!=NULL,"NULL SolverStrategyParameters object");

      d_iObjectId              = SolverStrategy::d_iInstanceId;
      d_iNumberIterations      = -1;
      d_dResidualNorm          = -1;
      d_dMaxRhs                = 1.0;
      d_bUseZeroInitialGuess = true;

      d_pOperator               = parameters->d_pOperator;

      SolverStrategy::d_iInstanceId++;
      SolverStrategy::getFromInput(parameters->d_db);   
    }

    SolverStrategy::~SolverStrategy()
    {
    }

    void 
      SolverStrategy::getFromInput(const boost::shared_ptr<AMP::Database>& db)
      {
        AMP_INSIST(db.get()!=NULL,"InputDatabase object must be non-NULL");

        d_iMaxIterations = db->getIntegerWithDefault("max_iterations", 1);

        d_dMaxError = db->getDoubleWithDefault("max_error", 1.0e-12);

        d_iDebugPrintInfoLevel = db->getIntegerWithDefault("print_info_level", 0);

        d_bUseZeroInitialGuess = db->getBoolWithDefault("zero_initial_guess", true);
      }

    void
      SolverStrategy::initialize(boost::shared_ptr< SolverStrategyParameters >  const parameters)
      {
        AMP_INSIST(parameters.get()!=NULL, "SolverStrategyParameters object cannot be NULL");
      }

    void
      SolverStrategy::resetOperator(const boost::shared_ptr<AMP::Operator::OperatorParameters> params)
      {
        if(d_pOperator.get()!=NULL)
        {
          d_pOperator->reset(params);
        }
      }

    void SolverStrategy::setConvergenceTolerance(
        const int max_iterations,
        const double max_error)
    {

      AMP_INSIST(max_iterations >= 0, "max_iterations must be non-negative");
      AMP_INSIST(max_error      >= 0.0, "max_eror must be non-negative");

      d_iMaxIterations = max_iterations;
      d_dMaxError      = max_error;
    }


  }
}

