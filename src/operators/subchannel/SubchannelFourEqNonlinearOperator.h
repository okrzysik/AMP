
#ifndef included_AMP_SubchannelFourEqNonlinearOperator
#define included_AMP_SubchannelFourEqNonlinearOperator

#include "AMP/operators/Operator.h"
#include "AMP/operators/subchannel/SubchannelOperatorParameters.h"

#include "AMP/ampmesh/MeshElementVectorIterator.h"

namespace AMP {
namespace Operator {

/**
  Nonlinear operator class for the 2-equation {enthalpy, pressure} formulation of the subchannel
  equations:
  see /AMPFuel-Docs/technicalInfo/flow/SubchannelFlow.pdf for details
  */
class SubchannelFourEqNonlinearOperator : public Operator
{
public:
    //! Constructor
    explicit SubchannelFourEqNonlinearOperator(
        std::shared_ptr<const SubchannelOperatorParameters> params );

    //! Destructor
    virtual ~SubchannelFourEqNonlinearOperator() {}

    //! Return the name of the operator
    std::string type() const override { return "SubchannelFourEqNonlinearOperator"; }

    /**
      For this operator we have an in-place apply.
      @param [in]  u input vector.
      @param [out] f output vector.
      */
    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    void reset( std::shared_ptr<const OperatorParameters> params ) override;

    std::shared_ptr<AMP::LinearAlgebra::Variable> getInputVariable() override
    {
        return d_inpVariable;
    }

    std::shared_ptr<AMP::LinearAlgebra::Variable> getOutputVariable() override
    {
        return d_outVariable;
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr
    subsetOutputVector( AMP::LinearAlgebra::Vector::shared_ptr vec ) override;
    virtual AMP::LinearAlgebra::Vector::const_shared_ptr
    subsetOutputVector( AMP::LinearAlgebra::Vector::const_shared_ptr vec ) override;

    virtual AMP::LinearAlgebra::Vector::shared_ptr
    subsetInputVector( AMP::LinearAlgebra::Vector::shared_ptr vec ) override;
    virtual AMP::LinearAlgebra::Vector::const_shared_ptr
    subsetInputVector( AMP::LinearAlgebra::Vector::const_shared_ptr vec ) override;

    void setVector( AMP::LinearAlgebra::Vector::shared_ptr frozenVec )
    {
        d_cladTemperature = frozenVec;
    }

    //! Get the element physics model
    std::shared_ptr<SubchannelPhysicsModel> getSubchannelPhysicsModel()
    {
        return d_subchannelPhysicsModel;
    }

    //! Get the Inlet Temperature [K]
    double getInletTemperature() { return d_Tin; }

    //! Get the Outlet Pressure [Pa]
    double getOutletPressure() { return d_Pout; }

    //! Get the current operator parameters
    auto getParams() { return d_params; }

    //! Makes map of lateral gaps to their centroids
    void getLateralFaces( AMP::Mesh::Mesh::shared_ptr,
                          std::map<AMP::Mesh::Point, AMP::Mesh::MeshElement> &,
                          std::map<AMP::Mesh::Point, AMP::Mesh::MeshElement> & );

    //! Makes map of gap widths to their xy positions
    std::map<AMP::Mesh::Point, double> getGapWidths( AMP::Mesh::Mesh::shared_ptr,
                                                     const std::vector<double> &,
                                                     const std::vector<double> &,
                                                     const std::vector<double> & );

    void getAxialFaces( const AMP::Mesh::MeshElement &,
                        AMP::Mesh::MeshElement &,
                        AMP::Mesh::MeshElement & );

    void fillSubchannelGrid(
        AMP::Mesh::Mesh::shared_ptr ); // function to fill the subchannel data for all processors

    int getSubchannelIndex( double x,
                            double y ); // function to give unique index for each subchannel

protected:
    //! Gets parameters from nonlinear operator for use in linear operator
    std::shared_ptr<OperatorParameters>
    getJacobianParameters( AMP::LinearAlgebra::Vector::const_shared_ptr u ) override;

    std::shared_ptr<const SubchannelOperatorParameters> d_params;

    std::shared_ptr<SubchannelPhysicsModel> d_subchannelPhysicsModel;

private:
    bool d_initialized;

    // Function used in reset to get double parameter or use default if missing
    double
    getDoubleParameter( std::shared_ptr<const SubchannelOperatorParameters>, std::string, double );

    // Function used in reset to get integer parameter or use default if missing
    int
    getIntegerParameter( std::shared_ptr<const SubchannelOperatorParameters>, std::string, int );

    // Function used in reset to get string parameter or use default if missing
    std::string getStringParameter( std::shared_ptr<const SubchannelOperatorParameters>,
                                    std::string,
                                    std::string );

    // Function used in reset to get bool parameter or use default if missing
    bool getBoolParameter( std::shared_ptr<const SubchannelOperatorParameters>, std::string, bool );

    std::shared_ptr<AMP::LinearAlgebra::Variable> d_inpVariable;
    std::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariable;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_cladTemperature;

    bool d_forceNoConduction; // option to force conduction terms to zero; used for testing
    bool d_forceNoTurbulence; // option to force turbulence terms to zero; used for testing
    bool d_forceNoHeatSource; // option to force heat source terms to zero; used for testing
    bool
        d_forceNoFriction; // option to force friction and form loss terms to zero; used for testing

    double d_Pout;           // exit pressure [Pa]
    double d_Tin;            // inlet temperature [K]
    double d_mass;           // inlet global mass flow rate [kg/s]
    double d_win;            // inlet lateral mass flow rate [kg/s]
    double d_gamma;          // fission heating coefficient
    double d_theta;          // channel angle [rad]
    double d_turbulenceCoef; // proportionality constant relating turbulent momentum to turbulent
                             // energy transport
    double d_reynolds;       // reynolds number
    double d_prandtl;        // prandtl number
    double d_KG;             // lateral form loss coefficient

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

    std::string d_source;            // heat source type
    std::string d_heatShape;         // heat shape used if heat source type is "totalHeatGeneration"
    double d_Q;                      // (sum of rod powers)/4 for each subchannel
    std::vector<double> d_QFraction; // fraction of max rod power in each subchannel

    std::vector<double> d_x, d_y, d_z;
    std::vector<bool> d_ownSubChannel; // does this processor own this subchannel (multiple
                                       // processors may own a subchannel)?
    std::vector<std::vector<AMP::Mesh::MeshElement>>
        d_subchannelElem; // List of elements in each subchannel
    std::vector<std::vector<AMP::Mesh::MeshElement>>
        d_subchannelFace;    // List of z-face elements in each subchannel
    size_t d_numSubchannels; // number of subchannels

    double Volume( double, double );              // evaluates specific volume
    double Temperature( double, double );         // evaluates temperature
    double ThermalConductivity( double, double ); // evaluates thermal conductivity
    double DynamicViscosity( double, double );    // evaluates dynamic viscosity
    double Enthalpy( double, double );            // evaluates specific enthalpy

    AMP::Mesh::MeshElement
    getAxiallyAdjacentLateralFace( AMP::Mesh::MeshElement *,
                                   const AMP::Mesh::MeshElement &,
                                   const std::map<AMP::Mesh::Point, AMP::Mesh::MeshElement> & );
};


} // namespace Operator
} // namespace AMP

#endif
