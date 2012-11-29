
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

      //! Constructor
      SubchannelFourEqNonlinearOperator(const boost::shared_ptr<SubchannelOperatorParameters> & params);

      //! Destructor
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

      void setVector(AMP::LinearAlgebra::Vector::shared_ptr &frozenVec) {
          d_cladTemperature = frozenVec;
      }
      
      //! Gets parameters from nonlinear operator for use in linear operator
      boost::shared_ptr<OperatorParameters> getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& );

      //! Get the element physics model
      boost::shared_ptr<SubchannelPhysicsModel> getSubchannelPhysicsModel() { return d_subchannelPhysicsModel; }
    
      //! Get the Inlet Temperature [K]
      double getInletTemperature() { return d_Tin; }

      //! Get the Outlet Pressure [Pa]
      double getOutletPressure() { return d_Pout; }

      //! Get the current operator parameters
      boost::shared_ptr<SubchannelOperatorParameters> getParams() { return d_params; }

      //! Makes map of lateral gaps to their centroids
      std::map<std::vector<double>,AMP::Mesh::MeshElement> getLateralFaces(AMP::Mesh::Mesh::shared_ptr);

      //! Makes map of gap widths to their xy positions
      std::map<std::vector<double>,double> getGapWidths(AMP::Mesh::Mesh::shared_ptr,
         const std::vector<double>&, const std::vector<double>&, const std::vector<double>&);

      void getAxialFaces(AMP::Mesh::MeshElement,AMP::Mesh::MeshElement&,AMP::Mesh::MeshElement&);

      void fillSubchannelGrid(AMP::Mesh::Mesh::shared_ptr); // function to fill the subchannel data for all processors

      int getSubchannelIndex( double x, double y ); // function to give unique index for each subchannel

    protected:

      boost::shared_ptr<SubchannelOperatorParameters> d_params;

      boost::shared_ptr<SubchannelPhysicsModel> d_subchannelPhysicsModel;

    private :

      bool d_initialized;

      // Function used in reset to get double parameter or use default if missing
      double getDoubleParameter(boost::shared_ptr<SubchannelOperatorParameters>, std::string, double);

      // Function used in reset to get integer parameter or use default if missing
      int getIntegerParameter(boost::shared_ptr<SubchannelOperatorParameters>, std::string, int);

      // Function used in reset to get double parameter or use default if missing
      std::string getStringParameter(boost::shared_ptr<SubchannelOperatorParameters>, std::string, std::string);

      boost::shared_ptr<AMP::LinearAlgebra::Variable> d_inpVariable;
      boost::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariable;
      boost::shared_ptr<AMP::LinearAlgebra::Vector> d_cladTemperature;

      double d_Pout;     // exit pressure [Pa]
      double d_Tin;      // inlet temperature [K]
      double d_mass;     // inlet global mass flow rate [kg/s]
      double d_win;      // inlet lateral mass flow rate [kg/s]
      double d_gamma;    // fission heating coefficient
      double d_theta;    // channel angle [rad]
      double d_turbulenceCoef; // proportionality constant relating turbulent momentum to turbulent energy transport
      double d_reynolds;  // reynolds number
      double d_prandtl;   // prandtl number

      std::string d_frictionModel; // friction model
      double d_friction; // friction factor
      double d_roughness; // surface roughness [m]

      std::vector<double> d_channelDiamFric;  // Channel hydraulic diameter using the clad perimeter
      std::vector<double> d_channelDiamHeat;  // Channel hydraulic diameter using the total perimeter
      std::vector<double> d_channelArea;  // Channel flow area
      std::vector<double> d_rodDiameter;  // Average rod diameter for each subchannel
      std::vector<double> d_rodFraction;  // Fraction of a rod in each subchannel
      std::vector<double> d_channelMass;  // Mass flow rate for each subchannel [kg/s]

      size_t d_NGrid;                 // number of grid spacers
      std::vector<double> d_zMinGrid; // z min positions of each grid spacer
      std::vector<double> d_zMaxGrid; // z max positions of each grid spacer
      std::vector<double> d_lossGrid; // loss coefficients for each grid spacer

      std::string d_source; // heat source type
      std::string d_heatShape; // heat shape used if heat source type is "totalHeatGeneration"
      double d_Q;        // (sum of rod powers)/4 for each subchannel
      std::vector<double> d_QFraction; // fraction of max rod power in each subchannel

      std::vector<double> d_x, d_y, d_z;
      std::vector<bool> d_ownSubChannel; // does this processor own this subchannel (multiple processors may own a subchannel)?
      std::vector< std::vector<AMP::Mesh::MeshElement> >  d_subchannelElem;   // List of elements in each subchannel
      std::vector< std::vector<AMP::Mesh::MeshElement> >  d_subchannelFace;   // List of z-face elements in each subchannel
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

