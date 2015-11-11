#ifndef MATERIAL_H
#define MATERIAL_H

#include <limits>
#include <map>
#include <string>
#include <vector>
#include "utils/shared_ptr.h"

#include "utils/Utilities.h"
#include "utils/Factory.h"

#include "Property.h"

// do not use property name with an embedded underscore
#define INSERT_PROPERTY_IN_MAP(name,space) \
	d_propertyMap->insert(std::make_pair( #name, AMP::shared_ptr<space::name##Prop>(new space::name##Prop)));

namespace AMP {
namespace Materials {

typedef AMP::shared_ptr<Property<double> > PropertyPtr;

/**
 * Material base class.
 * Loose organizer to collect a group of properties.
 */
class Material
{
public:
	Material(): d_propertyMap(NULL)
	{
	}

	virtual ~Material()
	{
		delete d_propertyMap;
	}

	/// specific shared pointer for this class
	typedef AMP::shared_ptr<Material> shared_ptr;

	static size_t counter;

public:

	/// check if a property exists in the material
	bool hasProperty(std::string type);

	/// get a pointer to a specific scalar property through its name
	PropertyPtr property(std::string type);

	/// return a list of all properties in this material
	std::vector<std::string> list();

protected:

	/// database of scalar properties
	std::map<std::string, PropertyPtr > *d_propertyMap;
};

/*// This macro is to be placed after each material class (UO2, Pu, etc.)
// It will register the material with the factory
#define REGISTER_MATERIAL(name)											\
namespace																\
{																		\
	AMP::voodoo::Registration<Material,name> reg(#name);				\
}*/

}
}


// Include all Material headers (they are responsible for registering themselves with the factory)
#include "Dr_nonlinear.h"
#include "Independent.h"
#include "Ox_MSRZC_09.h"
#include "Steel316_MSRZC_09.h"
#include "UO2_MSRZC_09.h"
#include "CylindricallySymmetric.h"
#include "WaterLibrary.h"
#include "FixedFuel.h"
#include "FixedClad.h"



#endif
