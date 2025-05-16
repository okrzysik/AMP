#ifndef included_AMP_OperatorParameters
#define included_AMP_OperatorParameters

#include "AMP/utils/ParameterBase.h"
#include "AMP/utils/Utilities.h"

#include <memory>


namespace AMP::Mesh {
class Mesh;
}


namespace AMP::Operator {

class Operator;

/**\class OperatorParameters
 *
 * OperatorParameters encapsulates parameters used to initialize or reset
 * operators. It is an abstract base class.
 */

class OperatorParameters : public ParameterBase
{
public:
    typedef std::shared_ptr<AMP::Operator::OperatorParameters> shared_ptr;

    /**
     * Construct and initialize a parameter list according to input
     * data.  Guess what the required and optional keywords are.
     */
    explicit OperatorParameters( std::shared_ptr<AMP::Database> db ) : ParameterBase( db )
    {
        if ( db ) {
            auto memLoc       = db->getWithDefault<std::string>( "MemoryLocation", "host" );
            d_memory_location = memoryLocationFromString( memLoc );
        } else if ( d_memory_location == AMP::Utilities::MemoryType::none ) {
            d_memory_location = AMP::Utilities::MemoryType::host;
        }
    }

    static AMP::Utilities::MemoryType memoryLocationFromString( const std::string &name )
    {
#ifdef USE_DEVICE
        if ( name == "managed" || name == "Managed" ) {
            return AMP::Utilities::MemoryType::managed;
        } else if ( name == "device" || name == "Device" ) {
            return AMP::Utilities::MemoryType::device;
        }
#endif
        (void) name;
        return AMP::Utilities::MemoryType::host;
    }

    /**
     * Destructor.
     */
    virtual ~OperatorParameters() {}

    std::shared_ptr<AMP::Mesh::Mesh> d_Mesh;
    /**
     * Allow for the case that a fully constructed operator is returned
     */
    std::shared_ptr<AMP::Operator::Operator> d_pOperator;

    /**
     * Location (host/managed/device) where internally created
     * vectors and matrices should live. Host memory should always
     * work. More specialized parameter classes can overwrite this
     * if supported.
     */
    AMP::Utilities::MemoryType d_memory_location = AMP::Utilities::MemoryType::none;
};


} // namespace AMP::Operator

#endif
