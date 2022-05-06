#ifndef MATERIAL_H
#define MATERIAL_H

#include "AMP/utils/Factory.h"
#include "Property.h"
#include <memory>

#include <map>
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


/*// This macro is to be placed after each material class (UO2, Pu, etc.)
// It will register the material with the factory
#define REGISTER_MATERIAL(name)											\
namespace																\
{																		\
    AMP::voodoo::Registration<Material,name> reg(#name);				\
}*/


} // namespace AMP::Materials


#endif
