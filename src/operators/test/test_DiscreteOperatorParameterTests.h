#include "../../utils/InputDatabase.h"
#include "../OperatorParameters.h"

namespace AMP {
namespace unit_test {

class InstantiateOperatorParameter
{
public:
    static const char *get_test_name() { return "instantiate OperatorParameters"; }

    template <typename UTILS>
    static void run_test( UTILS *utils )
    {
        AMP::shared_ptr<AMP::InputDatabase> new_db( new AMP::InputDatabase( "Dummy db" ) );
        AMP::Operator::OperatorParameters params( new_db );
        utils->passes( "instantiate OperatorParameters" );
    }
};
}
}
