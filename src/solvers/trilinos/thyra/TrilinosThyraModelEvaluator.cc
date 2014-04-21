#include "solvers/trilinos/thyra/TrilinosThyraModelEvaluator.h"
#include "solvers/trilinos/thyra/TrilinosLinearOP.h"
#include "vectors/trilinos/thyra/ThyraVectorSpaceWrapper.h"
#include "vectors/trilinos/thyra/ThyraVectorWrapper.h"
#include "vectors/trilinos/thyra/ThyraVector.h"
#include "operators/Operator.h"
#include "utils/Utilities.h"

#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_LinearOpBase.hpp"

namespace AMP {
namespace Solver {


/****************************************************************
*  Constructors                                                 *
****************************************************************/
TrilinosThyraModelEvaluator::TrilinosThyraModelEvaluator( boost::shared_ptr<TrilinosThyraModelEvaluatorParameters> params )
{
    AMP_ASSERT(params->d_nonlinearOp!=NULL);
    d_nonlinearOp = params->d_nonlinearOp;
    d_linearOp = params->d_linearOp;
    d_icVec = params->d_icVec;
    d_preconditioner = params->d_preconditioner;
    d_prePostOperator = params->d_prePostOperator;
}


/****************************************************************
*  Destructor                                                   *
****************************************************************/
TrilinosThyraModelEvaluator::~TrilinosThyraModelEvaluator()
{
}


/****************************************************************
*  Set the rhs vector                                           *
****************************************************************/
void TrilinosThyraModelEvaluator::setRhs( AMP::LinearAlgebra::Vector::const_shared_ptr rhs )
{
    d_rhs = rhs;
}


/****************************************************************
*  Evaluate the model                                           *
****************************************************************/
void TrilinosThyraModelEvaluator::evalModelImpl( const ::Thyra::ModelEvaluatorBase::InArgs<double> &inArgs,
    const ::Thyra::ModelEvaluatorBase::OutArgs<double> &outArgs ) const
{
         
    AMP::LinearAlgebra::Vector::const_shared_ptr x = AMP::LinearAlgebra::ThyraVector::constView( inArgs.get_x().get() );
    AMP::LinearAlgebra::Vector::shared_ptr   f_out = AMP::LinearAlgebra::ThyraVector::view( outArgs.get_f().get() );
    //const Thyra::ConstDetachedVectorView<double> x(inArgs.get_x());
    //const Teuchos::RCP< Thyra::VectorBase<double> > f_out = outArgs.get_f();
    //const Teuchos::RCP< Thyra::LinearOpBase<double> > W_out = outArgs.get_W_op();
    AMP_ASSERT(x!=NULL);

    // Temporary workaround to ensure x and rhs are consistent
    const_cast<AMP::LinearAlgebra::Vector*>(x.get())->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
    if ( d_rhs != NULL ) {
        const_cast<AMP::LinearAlgebra::Vector*>(d_rhs.get())->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
    }
    AMP_ASSERT(x->getUpdateStatus()==AMP::LinearAlgebra::Vector::UNCHANGED);
    AMP_ASSERT(d_rhs->getUpdateStatus()==AMP::LinearAlgebra::Vector::UNCHANGED);

    const Teuchos::RCP<Thyra::PreconditionerBase<double> > W_prec_out = outArgs.get_W_prec();
    if ( nonnull(W_prec_out) ) {
        // Reset the preconditioner
        AMP::LinearAlgebra::Vector::shared_ptr x2 = boost::const_pointer_cast<AMP::LinearAlgebra::Vector>(x);
        boost::shared_ptr<AMP::Operator::OperatorParameters> op_params = d_nonlinearOp->getJacobianParameters(x2);
        d_preconditioner->resetOperator(op_params);
    }

    if ( f_out != NULL ) {
        // Evaluate the residual:  r = A(u) - rhs
        f_out->zero();
        ::Thyra::ModelEvaluatorBase::Evaluation< ::Thyra::VectorBase<double> > eval = outArgs.get_f();
        bool exact = true;
        if (eval.getType() == ::Thyra::ModelEvaluatorBase::EVAL_TYPE_EXACT) {
            exact = true;
        } else if (eval.getType() == ::Thyra::ModelEvaluatorBase::EVAL_TYPE_APPROX_DERIV) {
            exact = false;
        } else if (eval.getType() == ::Thyra::ModelEvaluatorBase::EVAL_TYPE_VERY_APPROX_DERIV) {
            exact = false;
        }
        if ( d_prePostOperator != NULL )
            d_prePostOperator->runPreApply( x, f_out, exact );
        d_nonlinearOp->apply( d_rhs, x, f_out, 1.0, -1.0 );
        if ( d_prePostOperator != NULL )
            d_prePostOperator->runPostApply( x, f_out, exact );
    }

    if ( outArgs.supports(::Thyra::ModelEvaluatorBase::OUT_ARG_W_op) ) {
        Teuchos::RCP<Thyra::LinearOpBase<double> >  W_out = outArgs.get_W_op();
        if ( W_out.get() != NULL ) {
            // Get the jacobian
            AMP_ERROR("Not finished");
        }
    }
}


/****************************************************************
* Functions derived from Thyra::StateFuncModelEvaluatorBase     *
****************************************************************/
Teuchos::RCP<const ::Thyra::VectorSpaceBase<double> > TrilinosThyraModelEvaluator::get_x_space() const
{
    boost::shared_ptr<LinearAlgebra::ThyraVectorWrapper> vec( 
        new LinearAlgebra::ThyraVectorWrapper(std::vector<LinearAlgebra::Vector::shared_ptr>(1,d_icVec) ) );    
    Teuchos::RCP<LinearAlgebra::ThyraVectorSpaceWrapper> vector_space(new LinearAlgebra::ThyraVectorSpaceWrapper(vec));
    return vector_space;
}
Teuchos::RCP<const ::Thyra::VectorSpaceBase<double> > TrilinosThyraModelEvaluator::get_f_space() const
{
    boost::shared_ptr<LinearAlgebra::ThyraVectorWrapper> vec( 
        new LinearAlgebra::ThyraVectorWrapper(std::vector<LinearAlgebra::Vector::shared_ptr>(1,d_icVec) ) );    
    Teuchos::RCP<LinearAlgebra::ThyraVectorSpaceWrapper> vector_space(new LinearAlgebra::ThyraVectorSpaceWrapper(vec));
    return vector_space;
}
::Thyra::ModelEvaluatorBase::InArgs<double> TrilinosThyraModelEvaluator::getNominalValues() const
{
    return ::Thyra::ModelEvaluatorBase::InArgs<double>();
}
Teuchos::RCP< ::Thyra::LinearOpBase<double> > TrilinosThyraModelEvaluator::create_W_op() const
{
    return Teuchos::RCP<TrilinosLinearOP>( new TrilinosLinearOP(d_linearOp) );
}
void TrilinosThyraModelEvaluator::set_W_factory(const Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<double> >& W_factory)
{
    d_W_factory = W_factory;
}
Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<double> > TrilinosThyraModelEvaluator::get_W_factory() const
{
    return d_W_factory;
}
::Thyra::ModelEvaluatorBase::InArgs<double> TrilinosThyraModelEvaluator::createInArgs() const
{
    ::Thyra::ModelEvaluatorBase::InArgsSetup<double> inArgs;
    inArgs.setModelEvalDescription(this->description());
    inArgs.setSupports(::Thyra::ModelEvaluatorBase::IN_ARG_x);
    return inArgs;
}
::Thyra::ModelEvaluatorBase::OutArgs<double> TrilinosThyraModelEvaluator::createOutArgsImpl() const
{
    ::Thyra::ModelEvaluatorBase::OutArgsSetup<double> outArgs;
    outArgs.setModelEvalDescription(this->description());
    outArgs.setSupports(::Thyra::ModelEvaluatorBase::OUT_ARG_f);
    outArgs.setSupports(::Thyra::ModelEvaluatorBase::OUT_ARG_W_op);
    outArgs.setSupports(::Thyra::ModelEvaluatorBase::OUT_ARG_W_prec);
    return outArgs;
}


/****************************************************************
* Function to create the preconditioner                         *
****************************************************************/
Teuchos::RCP< ::Thyra::PreconditionerBase<double> > TrilinosThyraModelEvaluator::create_W_prec() const
{
    Teuchos::RCP<Thyra::DefaultPreconditioner<double> > preconditioner;
    if ( d_preconditioner != NULL ) {
        Teuchos::RCP<Thyra::LinearOpBase<double> > leftPrec;
        Teuchos::RCP<Thyra::LinearOpBase<double> > rightPrec( new TrilinosLinearOP(d_preconditioner) );
        preconditioner.reset( new Thyra::DefaultPreconditioner<double>( leftPrec, rightPrec ) );
    }
    return preconditioner;
}


/****************************************************************
* Return the views to the TrilinosLinearOP                      *
****************************************************************/
template<class T> static void nullDeleter( T* ) {};
boost::shared_ptr<AMP::Solver::TrilinosLinearOP> TrilinosThyraModelEvaluator::view( Teuchos::RCP< Thyra::LinearOpBase<double> > op )
{
    if ( op.is_null() )
        return boost::shared_ptr<AMP::Solver::TrilinosLinearOP>();
    AMP::Solver::TrilinosLinearOP* tmp = dynamic_cast<AMP::Solver::TrilinosLinearOP*>(op.get());
    AMP_ASSERT(tmp!=NULL);
    return boost::shared_ptr<AMP::Solver::TrilinosLinearOP>( tmp, nullDeleter<AMP::Solver::TrilinosLinearOP> );
}



}
}

