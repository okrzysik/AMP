#ifndef UO2_MSRZC_09_H
#define UO2_MSRZC_09_H

#include "Material.h"
#include "utils/Factory.h"


// Define the material
namespace AMP {
namespace Materials {
class UO2_MSRZC_09 : public Material
{
public:
    UO2_MSRZC_09();
};
}
}


// Add static initialize to force symbols to be included
// It will register the material with the factory
static struct UO2_MSRZC_09_INIT {
    UO2_MSRZC_09_INIT()
    {
        static AMP::voodoo::Registration<AMP::Materials::Material, AMP::Materials::UO2_MSRZC_09>
            reg( "UO2_MSRZC_09" );
    }
} UO2_MSRZC_09_init;


#endif
