#ifndef Ox_MSRZC_09_H
#define Ox_MSRZC_09_H

#include "Material.h"


// Define the material
namespace AMP {
namespace Materials {
class Ox_MSRZC_09 : public Material
{
public:
    Ox_MSRZC_09();
};
}
}


// Add static initialize to force symbols to be included
// It will register the material with the factory
static struct Ox_MSRZC_09_INIT {
    Ox_MSRZC_09_INIT()
    {
        static AMP::voodoo::Registration<AMP::Materials::Material, AMP::Materials::Ox_MSRZC_09> reg(
            "Ox_MSRZC_09" );
    }
} Ox_MSRZC_09_init;


#endif
