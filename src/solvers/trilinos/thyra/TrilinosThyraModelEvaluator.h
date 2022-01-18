#ifndef included_AMP_TrilinosThyraModelEvaluator
#define included_AMP_TrilinosThyraModelEvaluator


#include "AMP/discretization/DOF_Manager.h"
#include "AMP/solvers/trilinos/thyra/TrilinosLinearOP.h"
#include "AMP/solvers/trilinos/thyra/TrilinosThyraModelEvaluatorParameters.h"
#include <memory>


// Trilinos includes
DISABLE_WARNINGS
#include "Thyra_StateFuncModelEvaluatorBase.hpp"
ENABLE_WARNINGS


namespace AMP::Solver {


/**
 * The TrilinosThyraModelEvaluator is a wrapper for a Thyra ModelEvaluator to
 * wrap AMP::Operators for use with Trilinos NOX solvers.
 */
class TrilinosThyraModelEvaluator : public ::Thyra::StateFuncModelEvaluatorBase<double>
{
public:
    //! Default constructor
    explicit TrilinosThyraModelEvaluator(
        std::shared_ptr<TrilinosThyraModelEvaluatorParameters> params );

    //! Destructor
    virtual ~TrilinosThyraModelEvaluator();

    //! Copy constructor
    TrilinosThyraModelEvaluator( const TrilinosThyraModelEvaluator & ) = delete;

    //! Assignment operator
    TrilinosThyraModelEvaluator &operator=( const TrilinosThyraModelEvaluator & ) = delete;

    //! Function to set the rhs vector
    void setRhs( AMP::LinearAlgebra::Vector::const_shared_ptr rhs );

    // Functions derived from Thyra::StateFuncModelEvaluatorBase<double>
    virtual Teuchos::RCP<const ::Thyra::VectorSpaceBase<double>> get_x_space() const;
    virtual Teuchos::RCP<const ::Thyra::VectorSpaceBase<double>> get_f_space() const;
    virtual ::Thyra::ModelEvaluatorBase::InArgs<double> getNominalValues() const;
    virtual Teuchos::RCP<::Thyra::LinearOpBase<double>> create_W_op() const;
    virtual void set_W_factory(
        const Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<double>> &W_factory );
    virtual Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<double>> get_W_factory() const;
    virtual ::Thyra::ModelEvaluatorBase::InArgs<double> createInArgs() const;

    // Functions derived from Thyra::ModelEvaluator
    virtual Teuchos::RCP<::Thyra::PreconditionerBase<double>> create_W_prec() const;

protected:
    // Return TrilinosLinearOP from Thyra::LinearOpBase<double>
    static std::shared_ptr<AMP::Solver::TrilinosLinearOP>
    view( Teuchos::RCP<Thyra::LinearOpBase<double>> op );

    // Functions derived from Thyra::StateFuncModelEvaluatorBase<double>
    virtual ::Thyra::ModelEvaluatorBase::OutArgs<double> createOutArgsImpl() const;
    virtual void evalModelImpl( const ::Thyra::ModelEvaluatorBase::InArgs<double> &inArgs,
                                const ::Thyra::ModelEvaluatorBase::OutArgs<double> &outArgs ) const;

private:
    //! Empty constructor
    TrilinosThyraModelEvaluator() = default;

    // Data members
    Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<double>> d_W_factory;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_icVec;
    std::shared_ptr<const AMP::LinearAlgebra::Vector> d_rhs;
    std::shared_ptr<AMP::Operator::Operator> d_nonlinearOp;
    std::shared_ptr<AMP::Operator::Operator> d_linearOp;
    std::shared_ptr<AMP::Solver::SolverStrategy> d_preconditioner;
    std::shared_ptr<AMP::Solver::PrePostOperator> d_prePostOperator;
};
} // namespace AMP::Solver

#endif
