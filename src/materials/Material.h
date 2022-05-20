#ifndef MATERIAL_H
#define MATERIAL_H

#include "AMP/materials/Property.h"
#include "AMP/utils/Factory.h"
#include "AMP/utils/Utilities.h"

#include <map>
#include <memory>
#include <string>
#include <vector>


namespace AMP::Materials {


/**
 * Material base class.
 * Loose organizer to collect a group of properties.
 */
class Material
{
public:
    Material() {}

    virtual ~Material() {}

public:
    //! check if a property exists in the material
    bool hasProperty( std::string type );

    //! get a pointer to a specific scalar property
    std::shared_ptr<Property> property( std::string type );

    //! return a list of all properties in this material
    std::vector<std::string> list();

protected:
    /// database of scalar properties
    std::map<std::string, std::shared_ptr<Property>> d_propertyMap;
};


//! Get a material
std::shared_ptr<Material> getMaterial( const std::string &name );


//! Macro to register a material
#define registerMaterial( CLASS, NAME )                                                    \
    static struct CLASS##_INIT {                                                           \
        CLASS##_INIT()                                                                     \
        {                                                                                  \
            static AMP::voodoo::Registration<AMP::Materials::Material, CLASS> reg( NAME ); \
        }                                                                                  \
    } CLASS##_init


//! Macro to register a scalar property
#define registerScalarProperty( PROPERTY, VALUE, UNITS )        \
    d_propertyMap[PROPERTY] = std::make_shared<ScalarProperty>( \
        AMP::Utilities::demangle( typeid( *this ).name() ) + "::" + PROPERTY, VALUE, UNITS );


//! Scalar property class
class ScalarProperty final : public Property
{
public:
    ScalarProperty( std::string name, double value, AMP::Units unit )
        : Property( name ), d_value( value ), d_unit( unit )
    {
    }
    double eval( const std::vector<double> & ) override { return d_value; }

private:
    double d_value;
    AMP::Units d_unit;
};


} // namespace AMP::Materials


#endif
