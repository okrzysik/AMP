
#ifndef included_AMP_RobinVectorCorrection
#define included_AMP_RobinVectorCorrection

#include "BoundaryOperator.h"
#include "NeumannVectorCorrection.h"
#include "NeumannVectorCorrectionParameters.h"

/* Libmesh files */
#include "fe_type.h"
#include "fe_base.h"
#include "elem.h"
#include "quadrature.h"

#include "enum_order.h"
#include "enum_fe_family.h"
#include "enum_quadrature_type.h"
#include "auto_ptr.h"
#include "string_to_enum.h"

#include <string>

namespace AMP {
namespace Operator {

  typedef NeumannVectorCorrectionParameters RobinVectorCorrectionParameters;

  /**
    A class to impose Robin Boundary conditions for a nonlinear operator. This can 
    be written as \f$\alpha k(u)*\frac{\partial u}{\partial n} + \beta h*u = \gamma*c \f$.
    Imposing this condition would involve evaluating the expression and adding the 
    contribution to the residual vector. This class is derived from NeumannVectorCorrection
    as it implements similar functionality.
    */

  class RobinVectorCorrection : public NeumannVectorCorrection 
  {
    public :

      /**
        Constructor. This function reads all the parameters required for Robin boundary 
        conditions. Since it is derived from NeumannVectorCorrection, its constructor
        will be called to read the required parameters.
        */
      RobinVectorCorrection(const boost::shared_ptr<NeumannVectorCorrectionParameters> & params)
        : NeumannVectorCorrection (params)
      {
          reset(params);
          d_InstanceID = d_iInstance_id;
      }

      ~RobinVectorCorrection() { }

      /**
        Sets Robin values into the appropriate locations of the output vector (r). 
        */
      void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
              AMP::LinearAlgebra::Vector::shared_ptr  &r, const double a = -1.0, const double b = 1.0);

      /**
        This function can be used to change the Robin boundary conditions i.e., change the
        RHS flux values.
        */
      void reset(const boost::shared_ptr<OperatorParameters>& params);

      /**
        This function returns a parameter object that can be used to reset the corresponding
        RobinMatrixCorrection operator.
        */
      boost::shared_ptr<OperatorParameters> getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& );

    protected :

      // input variable for the unkown rhs
      AMP::LinearAlgebra::Variable::shared_ptr d_srcVariable;

      double d_hef;  //Convective Coefficient

      double d_alpha; 
      std::vector<double> d_beta;
      std::vector<double> d_gamma;
      double d_beta2; 
      double d_gamma2; 

      bool d_skipParams; 

      std::vector<AMP::LinearAlgebra::Vector::shared_ptr> d_elementInputVec;


    private :

      int d_InstanceID;
 
  };

}
}

#endif

