#ifndef Dr_nonlinear_H
#define Dr_nonlinear_H

#include "Material.h"


// Define the material
namespace AMP {
namespace Materials {
class Dr_nonlinear : public Material
{
public:
    Dr_nonlinear();
};
} // namespace Materials
} // namespace AMP


// Add static initialize to force symbols to be included
// It will register the material with the factory
static struct Dr_nonlinear_INIT {
    Dr_nonlinear_INIT()
    {
        static AMP::voodoo::Registration<AMP::Materials::Material, AMP::Materials::Dr_nonlinear>
            reg( "Dr_nonlinear" );
    }
} Dr_nonlinear_init;


#endif
