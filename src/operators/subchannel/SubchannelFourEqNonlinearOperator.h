
#ifndef included_AMP_SubchannelFourEqNonlinearOperator
#define included_AMP_SubchannelFourEqNonlinearOperator

#include "operators/Operator.h"
#include "operators/subchannel/SubchannelOperatorParameters.h"

#include "ampmesh/MeshElementVectorIterator.h"

namespace AMP {
namespace Operator {

  /**
    Nonlinear operator class for the 2-equation {enthalpy, pressure} formulation of the subchannel equations:
    see /AMPFuel-Docs/technicalInfo/flow/SubchannelFlow.pdf for details
    */
  class SubchannelFourEqNonlinearOperator : public Operator
  {
    public :

      /**
        Constructor
        */
      SubchannelFourEqNonlinearOperator(const boost::shared_ptr<SubchannelOperatorParameters> & params)
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
      virtual ~SubchannelFourEqNonlinearOperator() { }

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

      /**
        Makes map of lateral gaps to their centroids
        */
      std::map<std::vector<double>,AMP::Mesh::MeshElement> getLateralFaces(AMP::Mesh::Mesh::shared_ptr);

      void getAxialFaces(AMP::Mesh::MeshElement,AMP::Mesh::MeshElement&,AMP::Mesh::MeshElement&);

      void fillSubchannelGrid(AMP::Mesh::Mesh::shared_ptr); // function to fill the subchannel data for all processors

      int getSubchannelIndex( double x, double y ); // function to give unique index for each subchannel

    protected:

      boost::shared_ptr<SubchannelPhysicsModel> d_subchannelPhysicsModel;

    private :

      /**
        Function used in reset to get double parameter or use default if missing
        */
      double getDoubleParameter(boost::shared_ptr<SubchannelOperatorParameters>, std::string, double);

      /**
        Function used in reset to get integer parameter or use default if missing
        */
      int getIntegerParameter(boost::shared_ptr<SubchannelOperatorParameters>, std::string, int);

      /**
        Function used in reset to get double parameter or use default if missing
        */
      std::string getStringParameter(boost::shared_ptr<SubchannelOperatorParameters>, std::string, std::string);

      boost::shared_ptr<AMP::LinearAlgebra::Variable> d_inpVariable;

      boost::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariable;

      double d_Pout;     // exit pressure [Pa]
      double d_Tin;      // inlet temperature [K]
      double d_min;      // inlet mass flow rate [kg/s]
      double d_win;      // inlet mass flow rate [kg/s]
      double d_gamma;    // fission heating coefficient
      double d_theta;    // channel angle [rad]
      double d_pitch;    // lattice pitch [m]
      double d_diameter; // fuel rod diameter [m]
      double d_turbulenceCoef; // proportionality constant relating turbulent momentum to turbulent energy transport

      std::string d_frictionModel; // friction model
      double d_friction; // friction factor
      double d_roughness; // surface roughness [m]

      size_t d_NGrid;                 // number of grid spacers
      std::vector<double> d_zMinGrid; // z min positions of each grid spacer
      std::vector<double> d_zMaxGrid; // z max positions of each grid spacer
      std::vector<double> d_lossGrid; // loss coefficients for each grid spacer

      std::string d_source; // heat source type
      std::string d_heatShape; // heat shape used if heat source type is "totalHeatGeneration"
      double d_Q;        // (sum of rod powers)/4 for each subchannel
      std::vector<double> d_QFraction; // fraction of max rod power in each subchannel

      std::vector<double> d_x, d_y, d_z;
      std::vector<bool> d_ownSubchannel; // does this processor own this subchannel (multiple processors may own a subchannel)?
      size_t d_numSubchannels; // number of subchannels

      double Volume(double,double);              // evaluates specific volume
      double Temperature(double,double);         // evaluates temperature
      double ThermalConductivity(double,double); // evaluates thermal conductivity
      double DynamicViscosity(double,double);    // evaluates dynamic viscosity
      double Enthalpy(double,double);            // evaluates specific enthalpy

      AMP::Mesh::MeshElement getAxiallyAdjacentLateralFace(AMP::Mesh::MeshElement*,AMP::Mesh::MeshElement,
         std::map<std::vector<double>,AMP::Mesh::MeshElement>);
  };

}
}

#endif

