#ifndef included_TimeIntegratorParameters
#define included_TimeIntegratorParameters


#include "boost/shared_ptr.hpp"
#include "utils/Database.h"
#include "utils/ParameterBase.h"
#include "vectors/Vector.h"
#include "operators/Operator.h"

#include <string>


namespace AMP{
namespace TimeIntegrator{

/*!
  @brief TimeIntegratorParameters is a base class for providing
  parameters for the TimeIntegrator's. The Database object contained
  must contain the following entries:

  Required input keys and data types:
  @param initial_time double value for the initial simulation time.
  @param final_time double value for the final simulation time. 
  @param max_integrator_steps integer value for the maximum number
  of timesteps allowed.
  
  All input data items described above, except for initial_time,
  may be overwritten by new input values when continuing from restart.

 */

class TimeIntegratorParameters: public ParameterBase
{
public:
   //! Convience typedef
   typedef boost::shared_ptr<AMP::TimeIntegrator::TimeIntegratorParameters>  shared_ptr;

   TimeIntegratorParameters(const boost::shared_ptr<AMP::Database> db);

   virtual ~TimeIntegratorParameters();
   /**
   *  Database object which needs to be initialized specific to the time integrator.
   *  Documentation for parameters required by each integrator can be found in the
   *  documentation for the integrator.
   */
   boost::shared_ptr<AMP::Database> d_db;

   /**
    * String used to identify specific class instantiation of the time integrator
    */
   std::string d_object_name; 

   /**
    * Initial conditions vector
    */
   boost::shared_ptr<AMP::LinearAlgebra::Vector> d_ic_vector;

   /**
    * source term for time integration, can also include boundary conditions for IBVP problems
    */
   boost::shared_ptr<AMP::LinearAlgebra::Vector> d_pSourceTerm;

   /**
    * The operator is the right hand side operator for an explicit integrator when the time integration problem is : u_t = f(u)
    * but in the case of implicit time integrators the operator represents u_t-f(u) 
    */
   boost::shared_ptr< AMP::Operator::Operator > d_operator;

   /**
    * The operator is the left hand side mass operator (for FEM formulations)
    */
   boost::shared_ptr< AMP::Operator::Operator > d_pMassOperator;

   /**
    * algebraic variable
    */
   boost::shared_ptr<AMP::LinearAlgebra::Variable> d_pAlgebraicVariable;

protected:

private:
   // not implemented
   TimeIntegratorParameters(){}
   TimeIntegratorParameters(const TimeIntegratorParameters&);
   void operator=(const TimeIntegratorParameters&);
};

}
}

#endif

