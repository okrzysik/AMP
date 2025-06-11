#include "AMP/operators/OperatorParameters.h"


namespace AMP::Operator {


OperatorParameters::OperatorParameters( std::shared_ptr<AMP::Database> db,
                                        std::shared_ptr<AMP::Mesh::Mesh> mesh )
    : ParameterBase( db ), d_Mesh( mesh )
{
    if ( db ) {
        auto memLoc       = db->getWithDefault<std::string>( "MemoryLocation", "host" );
        d_memory_location = memoryLocationFromString( memLoc );
    }
}

AMP::Utilities::MemoryType OperatorParameters::memoryLocationFromString( const std::string &name )
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


} // namespace AMP::Operator
