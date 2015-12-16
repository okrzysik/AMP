
#ifndef included_AMP_SubchannelTwoEqNonlinearOperator
#define included_AMP_SubchannelTwoEqNonlinearOperator

#include "operators/Operator.h"
#include "operators/subchannel/SubchannelOperatorParameters.h"
#include "operators/subchannel/SubchannelPhysicsModel.h"

#include "ampmesh/MeshElementVectorIterator.h"

namespace AMP {
namespace Operator {

/**
  Nonlinear operator class for the 2-equation {enthalpy, pressure} formulation of the subchannel
  equations:
  see /AMPFuel-Docs/technicalInfo/flow/SubchannelFlow.pdf for details
  */
class SubchannelTwoEqNonlinearOperator : public Operator
{
public:
    //! Constructor
    explicit SubchannelTwoEqNonlinearOperator(
        const AMP::shared_ptr<SubchannelOperatorParameters> &params );

    //! Destructor
    virtual ~SubchannelTwoEqNonlinearOperator() {}

    /**
      For this operator we have an in-place apply.
      @param [in]  u input vector.
      @param [out] r residual/output vector.
     */
    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    void reset( const AMP::shared_ptr<OperatorParameters> &params ) override;

    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() override { return d_inpVariable; }

    AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() override { return d_outVariable; }

    void setVector( AMP::LinearAlgebra::Vector::shared_ptr frozenVec )
    {
        d_cladTemperature = frozenVec;
    }

    //! Get the element physics model
    AMP::shared_ptr<SubchannelPhysicsModel> getSubchannelPhysicsModel()
    {
        return d_subchannelPhysicsModel;
    }

    //! Get the Inlet Temperature [K]
    double getInletTemperature() { return d_Tin; }

    //! Get the Outlet Pressure [Pa]
    double getOutletPressure() { return d_Pout; }

    //! Get the current operator parameters
    AMP::shared_ptr<SubchannelOperatorParameters> getParams() { return d_params; }

protected:
    //! Gets parameters from nonlinear operator for use in linear operator
    AMP::shared_ptr<OperatorParameters>
    getJacobianParameters( AMP::LinearAlgebra::Vector::const_shared_ptr u ) override;

    AMP::shared_ptr<SubchannelOperatorParameters> d_params;

    AMP::shared_ptr<SubchannelPhysicsModel> d_subchannelPhysicsModel;

private:
    bool d_initialized;

    // Function used in reset to get double parameter or use default if missing
    double getDoubleParameter( AMP::shared_ptr<SubchannelOperatorParameters>, std::string, double );

    // Function used in reset to get integer parameter or use default if missing
    int getIntegerParameter( AMP::shared_ptr<SubchannelOperatorParameters>, std::string, int );

    // Function used in reset to get double parameter or use default if missing
    std::string getStringParameter( AMP::shared_ptr<SubchannelOperatorParameters>,
                                    std::string,
                                    std::string );

    AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_inpVariable;
    AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariable;
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_cladTemperature;

    double d_Pout;     // exit pressure [Pa]
    double d_Tin;      // inlet temperature [K]
    double d_mass;     // inlet global mass flow rate [kg/s]
    double d_gamma;    // fission heating coefficient
    double d_theta;    // channel angle [rad]
    double d_Q;        // rod power
    double d_reynolds; // reynolds number
    double d_prandtl;  // prandtl number

    std::string d_frictionModel; // friction model
    double d_friction;           // friction factor
    double d_roughness;          // surface roughness [m]

    std::vector<double> d_channelDiam; // Channel hydraulic diameter using the wetted perimeter
    std::vector<double> d_channelArea; // Channel flow area
    std::vector<double> d_rodDiameter; // Average rod diameter for each subchannel
    std::vector<double> d_rodFraction; // Fraction of a rod in each subchannel
    std::vector<double> d_channelMass; // Mass flow rate for each subchannel [kg/s]

    size_t d_NGrid;                 // number of grid spacers
    std::vector<double> d_zMinGrid; // z min positions of each grid spacer
    std::vector<double> d_zMaxGrid; // z max positions of each grid spacer
    std::vector<double> d_lossGrid; // loss coefficients for each grid spacer

    std::string d_source;    // heat source type
    std::string d_heatShape; // heat shape used if heat source type is "totalHeatGeneration"

    std::vector<double> d_x, d_y, d_z;
    std::vector<bool> d_ownSubChannel; // Which subchannels do I own
    std::vector<std::vector<AMP::Mesh::MeshElement>>
        d_subchannelElem; // List of elements in each subchannel
    std::vector<std::vector<AMP::Mesh::MeshElement>>
        d_subchannelFace; // List of z-face elements in each subchannel
    int getSubchannelIndex( double x, double y );
    size_t d_numSubchannels;
};
}
}

#endif
