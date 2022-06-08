#include "AMP/materials/Material.h"
#include "AMP/materials/MaterialList.h"
#include "AMP/materials/ScalarProperty.h"


namespace AMP::Materials {


// check if a property exists in the material
bool Material::hasProperty( const std::string &type ) const
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
std::vector<std::string> Material::list() const
{
    std::vector<std::string> result;
    for ( auto it : d_propertyMap )
        result.push_back( it.first );
    return result;
}


// Add a property
void Material::addScalarProperty( std::string name,
                                  double value,
                                  const AMP::Units &unit,
                                  std::string source )
{
    addProperty<ScalarProperty>( std::move( name ), value, unit, std::move( source ) );
}
void Material::addScalarProperty( std::string name,
                                  AMP::Array<double> value,
                                  const AMP::Units &unit,
                                  std::string source )
{
    addProperty<ScalarProperty>( std::move( name ), std::move( value ), unit, std::move( source ) );
}
void Material::addPolynomialProperty( std::string name,
                                      std::string source,
                                      const AMP::Units &unit,
                                      std::vector<double> params,
                                      std::vector<std::string> args,
                                      std::vector<std::array<double, 2>> ranges,
                                      std::vector<AMP::Units> argUnits )
{
    addProperty<PolynomialProperty>( std::move( name ),
                                     std::move( source ),
                                     unit,
                                     std::move( params ),
                                     std::move( args ),
                                     std::move( ranges ),
                                     std::move( argUnits ) );
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
