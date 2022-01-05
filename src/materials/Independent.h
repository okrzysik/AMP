#ifndef Independent_H
#define Independent_H

#include "Material.h"


// Define the material
namespace AMP::Materials {
class Independent : public Material
{
public:
    Independent();
};
} // namespace AMP::Materials


// Add static initialize to force symbols to be included
// It will register the material with the factory
static struct Independent_INIT {
    Independent_INIT()
    {
        static AMP::voodoo::Registration<AMP::Materials::Material, AMP::Materials::Independent> reg(
            "Independent" );
    }
} Independent_init;


#endif
