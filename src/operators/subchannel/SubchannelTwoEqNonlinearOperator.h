
#ifndef included_AMP_SubchannelTwoEqNonlinearOperator
#define included_AMP_SubchannelTwoEqNonlinearOperator

#include "operators/Operator.h"
#include "SubchannelOperatorParameters.h"

#include "ampmesh/MeshElementVectorIterator.h"

namespace AMP {
namespace Operator {

  /**
    Nonlinear operator class for the 2-equation {enthalpy, pressure} formulation of the subchannel equations:
    see /AMPFuel-Docs/technicalInfo/flow/SubchannelFlow.pdf for details
    */
  class SubchannelTwoEqNonlinearOperator : public Operator
  {
    public :

      /**
        Constructor
        */
      SubchannelTwoEqNonlinearOperator(const boost::shared_ptr<SubchannelOperatorParameters> & params)
        : Operator (params)
      {
        AMP_INSIST( params->d_db->keyExists("InputVariable"), "Key 'InputVariable' does not exist");
        std::string inpVar = params->d_db->getString("InputVariable");
        d_inpVariable.reset(new AMP::LinearAlgebra::Variable(inpVar));

        AMP_INSIST( params->d_db->keyExists("OutputVariable"), "Key 'OutputVariable' does not exist");
        std::string outVar = params->d_db->getString("OutputVariable");
        d_outVariable.reset(new AMP::LinearAlgebra::Variable(outVar));

        reset(params);
      }

      /**
        Destructor
        */
      ~SubchannelTwoEqNonlinearOperator() { }

      /**
        For this operator we have an in-place apply.
        @param [in]  f auxillary/rhs vector. 
        @param [in]  u input vector. 
        @param [out] r residual/output vector. 
        @param [in]  a first constant used in the expression: r = a*A(u) + b*f. The default value is -1.
        @param [in]  b second constant used in the expression: r = a*A(u) + b*f. The default value is 1.
        */
      void apply(AMP::LinearAlgebra::Vector::const_shared_ptr f, AMP::LinearAlgebra::Vector::const_shared_ptr u,
          AMP::LinearAlgebra::Vector::shared_ptr r, const double a = -1.0, const double b = 1.0);

      void reset(const boost::shared_ptr<OperatorParameters>& params);

      AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() {
        return d_inpVariable;
      }

      AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
        return d_outVariable;
      }

      virtual AMP::LinearAlgebra::Vector::shared_ptr subsetOutputVector(AMP::LinearAlgebra::Vector::shared_ptr vec);
      virtual AMP::LinearAlgebra::Vector::const_shared_ptr subsetOutputVector(AMP::LinearAlgebra::Vector::const_shared_ptr vec);

      virtual AMP::LinearAlgebra::Vector::shared_ptr subsetInputVector(AMP::LinearAlgebra::Vector::shared_ptr vec);
      virtual AMP::LinearAlgebra::Vector::const_shared_ptr subsetInputVector(AMP::LinearAlgebra::Vector::const_shared_ptr vec);

      /**
        Gets parameters from nonlinear operator for use in linear operator
        */
      boost::shared_ptr<OperatorParameters> getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& );

    protected:

      boost::shared_ptr<SubchannelOperatorParameters> d_params;

      boost::shared_ptr<SubchannelPhysicsModel> d_subchannelPhysicsModel;

    private :

      /**
        Function used in reset to get double parameter or use default if missing
        */
      double getDoubleParameter(boost::shared_ptr<SubchannelOperatorParameters>, std::string, double);

      /**
        Function used in reset to get double parameter or use default if missing
        */
      std::string getStringParameter(boost::shared_ptr<SubchannelOperatorParameters>, std::string, std::string);

      boost::shared_ptr<AMP::LinearAlgebra::Variable> d_inpVariable;

      boost::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariable;
      
      boost::shared_ptr<AMP::LinearAlgebra::Vector> d_cladTemperature;

      std::vector<double> d_channelFractions;     // channel fraction to calculate effective radius 

      double d_Pout;     // exit pressure [Pa]
      double d_Tin;      // inlet temperature [K]
      double d_m;        // inlet mass flow rate [kg/s]
      double d_gamma;    // fission heating coefficient
      double d_theta;    // channel angle [rad]
      double d_friction; // friction factor
      double d_pitch;    // lattice pitch [m]
      double d_diameter; // fuel rod diameter [m]
      double d_K;        // form loss coefficient
      double d_Q;        // rod power

      double d_channelDia;
      double d_reynolds  ;
      double d_prandtl   ;

      std::string d_source; // heat source type
      std::string d_heatShape; // heat shape used if heat source type is "totalHeatGeneration"

      std::vector<double> d_x, d_y, d_z;
      std::vector<bool> d_ownSubChannel;                      // Which subchannels do I own (multple procs my own a subchannel)
      int getSubchannelIndex( double x, double y );
      void fillSubchannelGrid(AMP::Mesh::Mesh::shared_ptr);   // Function to fill the subchannel data for all processors
      int d_numSubchannels; 
  };

}
}

#endif

