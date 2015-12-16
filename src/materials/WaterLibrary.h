#ifndef WaterLibrary_H
#define WaterLibrary_H

#include "Material.h"
#include "utils/Factory.h"


// Define the material
namespace AMP {
namespace Materials {
class WaterLibrary : public Material {
public:
    WaterLibrary();
};
}
}


// Add static initialize to force symbols to be included
// It will register the material with the factory
static struct WaterLibrary {
    WaterLibrary()
    {
        static AMP::voodoo::Registration<AMP::Materials::Material, AMP::Materials::WaterLibrary>
            reg( "WaterLibrary" );
    }
} WaterLibrary_init;


#endif
