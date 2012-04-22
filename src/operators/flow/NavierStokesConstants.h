
#ifndef included_AMP_NavierStokesConstants
#define included_AMP_NavierStokesConstants

namespace AMP {
namespace Operator {
  namespace NavierStokes {
    const unsigned int VELOCITY = 0; 
    const unsigned int PRESSURE = 1; 
    const unsigned int TEMPERATURE = 2; 
    const unsigned int PRINCIPALSTRESS = 3; 
    const unsigned int SHEARSTRESS = 4; 
    const unsigned int TOTAL_NUMBER_OF_VARIABLES = 5;
  }
}
}

#endif

