#include "Material.h"

// Include all Material headers (they are responsible for registering themselves with the factory)
#include "CylindricallySymmetric.h"
#include "Dr_nonlinear.h"
#include "FixedClad.h"
#include "FixedFuel.h"
#include "Independent.h"
#include "Ox_MSRZC_09.h"
#include "Steel316_MSRZC_09.h"
#include "UO2_MSRZC_09.h"
#include "WaterLibrary.h"


namespace AMP {
namespace Materials {


// check if a property exists in the material
bool Material::hasProperty( std::string type )
{
    return d_propertyMap->find( type ) != d_propertyMap->end();
}


// get a pointer to a specific property through its name
std::shared_ptr<Property<double>> Material::property( std::string type )
{
    auto it = d_propertyMap->find( type );
    AMP_INSIST( it != d_propertyMap->end(), std::string( "property " ) + type + " is not defined" );
    return it->second;
}


// return a list of all properties in this material
// note: this list is in error if property has an embedded underscore
std::vector<std::string> Material::list()
{
    std::vector<std::string> result;
    for ( auto it : *d_propertyMap ) {
        std::string name        = it.second->get_name();
        size_t usIndex          = name.rfind( "_" );
        std::string nameReduced = name.substr( usIndex + 1 );
        result.push_back( nameReduced );
    }
    return result;
}
} // namespace Materials
} // namespace AMP
