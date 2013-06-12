// This file containts some useful constants for subchannel
#ifndef included_AMP_SubchannelConstants
#define included_AMP_SubchannelConstants

namespace AMP {
namespace Operator {
namespace Subchannel {

const double scaleAxialMassFlowRate = 6.2e-4;      // Scale the axial mass flow rate by this constant in the vector (controls the norm)
const double scaleEnthalpy = 2.6e3;               // Scale the enthalapy by this constant in the vector (controls the norm)
const double scalePressure = 3.1e4;               // Scale the pressure by this constant in the vector (controls the norm)
const double scaleLateralMassFlowRate = 1.0e-7;    // Scale the lateral mass flow rate by this constant in the vector (controls the norm)

}
}
}

#endif

