// This file containts some useful constants for subchannel
#ifndef included_AMP_SubchannelConstants
#define included_AMP_SubchannelConstants

namespace AMP {
namespace Operator {
namespace Subchannel {

const double scaleAxialMassFlowRate = 1.0;      // Scale the axial mass flow rate by this constant in the vector (controls the norm)
const double scaleEnthalpy = 1.0;               // Scale the enthalapy by this constant in the vector (controls the norm)
const double scalePressure = 1.0;               // Scale the pressure by this constant in the vector (controls the norm)
const double scaleLateralMassFlowRate = 1.0;    // Scale the lateral mass flow rate by this constant in the vector (controls the norm)

}
}
}

#endif

