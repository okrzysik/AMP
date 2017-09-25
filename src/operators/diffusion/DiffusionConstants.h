#ifndef included_AMP_DiffusionConstants
#define included_AMP_DiffusionConstants

#include <string>

namespace AMP {
namespace Operator {

/**
 * Namespace for all of the diffusion operator code.
 */
namespace Diffusion {
//// These constants define indices in an array of standard material parameters used as an argument
/// to a material
/// property evaluation function
//@{
const unsigned int TEMPERATURE      = 0;
const unsigned int CONCENTRATION    = 1;
const unsigned int BURNUP           = 2;
const unsigned int NUMBER_VARIABLES = 3; //< Total number of material parameters.
//@}

/// ASCI names of the standard material parameters appearing as arguments to material property
/// evaluation functions.
const std::string names[NUMBER_VARIABLES] = { "temperature", "concentration", "burnup" };
} // namespace Diffusion
} // namespace Operator
} // namespace AMP

#endif
