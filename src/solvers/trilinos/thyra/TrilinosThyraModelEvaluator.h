#ifndef included_AMP_TrilinosThyraModelEvaluator
#define included_AMP_TrilinosThyraModelEvaluator


#include "boost/shared_ptr.hpp"
#include "solvers/trilinos/thyra/TrilinosLinearOP.h"
#include "solvers/trilinos/thyra/TrilinosThyraModelEvaluatorParameters.h"
#include "discretization/DOF_Manager.h"


// Trilinos includes
#include "Thyra_StateFuncModelEvaluatorBase.hpp"


namespace AMP {
namespace Solver {


/**
  * The TrilinosThyraModelEvaluator is a wrapper for a Thyra ModelEvaluator to 
  * wrap AMP::Operators for use with Trilinos NOX solvers.
  */
class TrilinosThyraModelEvaluator: public ::Thyra::StateFuncModelEvaluatorBase<double>
{
public:
    
    //! Default constructor
    TrilinosThyraModelEvaluator( boost::shared_ptr<TrilinosThyraModelEvaluatorParameters> params );

    //! Destructor
    virtual ~TrilinosThyraModelEvaluator();

    //! Function to set the rhs vector
    void setRhs( AMP::LinearAlgebra::Vector::const_shared_ptr rhs );

    // Functions derived from Thyra::StateFuncModelEvaluatorBase<double>
    virtual Teuchos::RCP<const ::Thyra::VectorSpaceBase<double> > get_x_space() const;
    virtual Teuchos::RCP<const ::Thyra::VectorSpaceBase<double> > get_f_space() const;
    virtual ::Thyra::ModelEvaluatorBase::InArgs<double> getNominalValues() const;
    virtual Teuchos::RCP< ::Thyra::LinearOpBase<double> > create_W_op() const;
    virtual void set_W_factory(const Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<double> >& W_factory);
    virtual Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<double> > get_W_factory() const;
    virtual ::Thyra::ModelEvaluatorBase::InArgs<double> createInArgs() const;

protected:

    // Return TrilinosLinearOP from Thyra::LinearOpBase<double>
    static boost::shared_ptr<AMP::Solver::TrilinosLinearOP> view( Teuchos::RCP< Thyra::LinearOpBase<double> > op );

    // Functions derived from Thyra::StateFuncModelEvaluatorBase<double>
    virtual ::Thyra::ModelEvaluatorBase::OutArgs<double> createOutArgsImpl() const;
    virtual void evalModelImpl( const ::Thyra::ModelEvaluatorBase::InArgs<double> &inArgs,
        const ::Thyra::ModelEvaluatorBase::OutArgs<double> &outArgs ) const;

private:

    //! Empty constructor
    TrilinosThyraModelEvaluator() {}

    // Data members
    Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<double> > d_W_factory;
    AMP::LinearAlgebra::Vector::shared_ptr d_icVec;
    AMP::LinearAlgebra::Vector::const_shared_ptr d_rhs;
    AMP::Operator::Operator::shared_ptr d_nonlinearOp;
    AMP::Operator::Operator::shared_ptr d_linearOp;

};


}
}

#endif

