#ifndef included_AMP_OperatorParameters
#define included_AMP_OperatorParameters

#include "AMP/utils/ParameterBase.h"
#include "AMP/utils/Utilities.h"

#include <memory>


namespace AMP::Mesh {
class Mesh;
}


namespace AMP::Operator {


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
        if ( name != "host" && name != "Host" ) {
            AMP_WARNING( "Unrecognized memory location, returning host" );
        }
        return AMP::Utilities::MemoryType::host;
    }

    /**
     * Destructor.
     */
    virtual ~OperatorParameters() {}

    std::shared_ptr<AMP::Mesh::Mesh> d_Mesh;

    /**
     * Location (host/managed/device) where internally created
     * vectors and matrices should live. Host memory should always
     * work. More specialized parameter classes can overwrite this
     * if supported.
     */
    AMP::Utilities::MemoryType d_memory_location = AMP::Utilities::MemoryType::host;
};


} // namespace AMP::Operator

#endif
