#ifndef included_IDATimeOperator
#define included_IDATimeOperator


#include "boost/shared_ptr.hpp"
#include "operators/Operator.h"
#include "operators/OperatorParameters.h"
#include "operators/libmesh/MassLinearFEOperator.h"
#include "vectors/Vector.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"
#include "time_integrators/TimeOperatorParameters.h"
#include "time_integrators/TimeOperator.h"
#include "operators/OperatorBuilder.h"

// BP : the next include is probably unnecessary
#include "operators/libmesh/VolumeIntegralOperator.h"

namespace AMP{
namespace TimeIntegrator{
  
  typedef  TimeOperatorParameters IDATimeOperatorParameters;
  
  /*!
    @brief operator class associated with IDATimeIntegrator
    
    Class IDATimeOperator is derived from TimeOperator. It
    is the operator class associated with a IDATimeIntegrator.
    
    @see IDATimeIntegrator
    @see TimeOperator
  */
  
  class IDATimeOperator: public TimeOperator
  {
  public:

    /**
     * Main constructor.
     @param [in] params: shared pointer to TimeOperatorParameters object.
     */
    IDATimeOperator(boost::shared_ptr<AMP::Operator::OperatorParameters > params);

    /**
     * virtual destructor
     */
    virtual ~IDATimeOperator();
    
    //virtual void reset(const boost::shared_ptr<AMP::Operator::OperatorParameters>& params);
    
      /**
        The function that computes the residual.
       * @param f: rhs vector for A(u)=f, this may be a null pointer if f=0. 
       * @param u: multivector of the state.
       * @param r: specific power in Watts per gram 
       * @param a: constnt multiplier applied to return of operator
       * @param b: constant multiplier applied to return of rhs vector
       The result of apply is
       * r = b*f+a*A(u)
       */
    void apply(AMP::LinearAlgebra::Vector::const_shared_ptr f, 
           AMP::LinearAlgebra::Vector::const_shared_ptr u,
           AMP::LinearAlgebra::Vector::shared_ptr r,
           const double a = -1.0, const double b=1.0);

    /**
     * registers the time derivative vector provided by IDA with this operator
     @param [in] vec   shared pointer to time derivative computed by IDA
     */
    void registerIDATimeDerivative(boost::shared_ptr<AMP::LinearAlgebra::Vector> vec) {d_pIDATimeDerivative = vec; }

    /**
     * registers a source term if any
     @param [in] vec   shared pointer to vector for source term
     */
    void registerSourceTerm(boost::shared_ptr<AMP::LinearAlgebra::Vector> vec) {d_pSourceTerm = vec; }

    /**
     * sets the current time for the operator
     @param [in] currentTime   the current time
     */
    void registerCurrentTime( double currentTime ) {d_current_time = currentTime;}        
    
  protected:

    IDATimeOperator();
    
    boost::shared_ptr<AMP::LinearAlgebra::Vector> d_pIDATimeDerivative;
    
    bool d_cloningHappened;
    
    //JL
    //The test we want to run has a source term which depends on time
    //The time comes from TimeIntegrator
    double d_current_time; 
    double d_beta;
    
  private:
    
    
  };
  
}
}

#endif

