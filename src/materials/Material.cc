/*
 * This is a dummy file to pacify cmake.
 * Material.cc
 *
 *  Created on: Aug 25, 2011
 *      Author: gad
 */

#include "Material.h"

namespace AMP
{
namespace Materials
{

size_t Material::counter=0;


// check if a property exists in the material
bool Material::hasProperty(std::string type) 
{
	return d_propertyMap->find(type) != d_propertyMap->end();
}


// get a pointer to a specific property through its name
PropertyPtr Material::property(std::string type)
{
	std::map<std::string, PropertyPtr>::iterator it;
	it = d_propertyMap->find(type);
	AMP_INSIST(it != d_propertyMap->end(),
			std::string("property ") + type + " is not defined");
	return it->second;
}


// return a list of all properties in this material
// note: this list is in error if property has an embedded underscore
std::vector<std::string> Material::list()
{
	std::vector<std::string> result;
	for (std::map<std::string, PropertyPtr>::iterator it=d_propertyMap->begin();
			it != d_propertyMap->end(); ++it)
	{
		std::string name = it->second->get_name();
		size_t usIndex = name.rfind("_"); 
		std::string nameReduced = name.substr(usIndex+1);

		result.push_back(nameReduced);
	}
	return result;
}

}
}
