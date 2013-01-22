#ifndef included_AMP_TrilinosThyraModelEvaluator
#define included_AMP_TrilinosThyraModelEvaluator


#include "boost/shared_ptr.hpp"
#include "solvers/TrilinosLinearOP.h"


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
    
    //! Empty constructor
    TrilinosThyraModelEvaluator();

    //! Destructor
    virtual ~TrilinosThyraModelEvaluator();

    // Functions derived from Thyra::StateFuncModelEvaluatorBase<double>
    virtual Teuchos::RCP<const ::Thyra::VectorSpaceBase<double> > get_x_space() const;
    virtual Teuchos::RCP<const ::Thyra::VectorSpaceBase<double> > get_f_space() const;
    virtual ::Thyra::ModelEvaluatorBase::InArgs<double> getNominalValues() const;
    virtual Teuchos::RCP< ::Thyra::LinearOpBase<double> > create_W_op() const;
    virtual void set_W_factory(const Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<double> >& W_factory);
    virtual Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<double> > get_W_factory() const;
    virtual ::Thyra::ModelEvaluatorBase::InArgs<double> createInArgs() const;

private:

    // Functions derived from Thyra::StateFuncModelEvaluatorBase<double>
    virtual ::Thyra::ModelEvaluatorBase::OutArgs<double> createOutArgsImpl() const;
    virtual void evalModelImpl( const ::Thyra::ModelEvaluatorBase::InArgs<double> &inArgs,
        const ::Thyra::ModelEvaluatorBase::OutArgs<double> &outArgs ) const;

    // Data members
    Teuchos::RCP<TrilinosLinearOP> d_linearOP;
    Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<double> > d_W_factory;

};


}
}

#endif

