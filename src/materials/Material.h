#ifndef included_AMP_Material
#define included_AMP_Material

#include "AMP/materials/Property.h"
#include "AMP/utils/Factory.h"
#include "AMP/utils/Utilities.h"

#include <map>
#include <memory>
#include <string>
#include <vector>


// This macro is to be placed after each material class (UO2, Pu, etc.)
// It will register the material with the factory
#define REGISTER_MATERIAL( NAME )                                                          \
    static struct NAME##_INIT {                                                            \
        NAME##_INIT()                                                                      \
        {                                                                                  \
            static AMP::voodoo::Registration<AMP::Materials::Material, NAME> reg( #NAME ); \
        }                                                                                  \
    } NAME##_init


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

//! Get the list of materials available
std::vector<std::string> getMaterialList();


} // namespace AMP::Materials


#endif
