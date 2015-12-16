#ifndef Steel316_MSRZC_09_H
#define Steel316_MSRZC_09_H

#include "Material.h"
#include "utils/Factory.h"


// Define the material
namespace AMP {
namespace Materials {
class Steel316_MSRZC_09 : public Material
{
public:
    Steel316_MSRZC_09();
};
}
}


// Add static initialize to force symbols to be included
// It will register the material with the factory
static struct Steel316_MSRZC_09 {
    Steel316_MSRZC_09()
    {
        static AMP::voodoo::Registration<AMP::Materials::Material,
                                         AMP::Materials::Steel316_MSRZC_09>
            reg( "Steel316_MSRZC_09" );
    }
} Steel316_MSRZC_09_init;


#endif
