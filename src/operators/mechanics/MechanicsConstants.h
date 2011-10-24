
#ifndef included_AMP_MechanicsConstants
#define included_AMP_MechanicsConstants

namespace AMP {
  namespace Operator {
    namespace Mechanics {
      const unsigned int DISPLACEMENT = 0; /**< Global constant used to identify displacement variables within mechanics operators. */
      const unsigned int TEMPERATURE = 1; /**<  Global constant used to identify temperature variables within mechanics operators. */
      const unsigned int BURNUP = 2; /**< Global constant used to indentify burnup variables within mechanics operators. */
      const unsigned int OXYGEN_CONCENTRATION = 3; /**< Global constant used to indentify oxygen concentration variables within mechanics operators. */
      const unsigned int LHGR = 4; /**< Global constant used to identify Linear Heat Generation Rate variable within mechanics operators (needed for relocation). */
      const unsigned int TOTAL_NUMBER_OF_VARIABLES = 5; /**< Global constant that equals the number of 
                                                          different types of variables used within mechanics operators. */
    }
  }
}

#endif

