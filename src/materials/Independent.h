#ifndef Independent_H
#define Independent_H

#include "Material.h"
#include "utils/Factory.h"


// Define the material
namespace AMP {
namespace Materials {
class Independent : public Material {
public:
    Independent();
};
}
}


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
