#include "AMP/materials/Material.h"
#include "AMP/materials/MaterialList.h"


namespace AMP::Materials {


// check if a property exists in the material
bool Material::hasProperty( std::string type )
{
    return d_propertyMap.find( type ) != d_propertyMap.end();
}


// get a pointer to a specific property through its name
std::shared_ptr<Property> Material::property( std::string type )
{
    auto it = d_propertyMap.find( type );
    AMP_INSIST( it != d_propertyMap.end(), std::string( "property " ) + type + " is not defined" );
    return it->second;
}


// return a list of all properties in this material
// note: this list is in error if property has an embedded underscore
std::vector<std::string> Material::list()
{
    std::vector<std::string> result;
    for ( auto it : d_propertyMap )
        result.push_back( it.first );
    return result;
}


// Return the material
std::shared_ptr<Material> getMaterial( const std::string &name )
{
    return AMP::voodoo::Factory<AMP::Materials::Material>::instance().create( name );
}


// Return a list of available materials
std::vector<std::string> getMaterialList()
{
    return AMP::voodoo::Factory<AMP::Materials::Material>::instance().getKeys();
}


} // namespace AMP::Materials
