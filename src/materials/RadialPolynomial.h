#ifndef RadialPolynomial_H
#define RadialPolynomial_H

#include "materials/Material.h"
#include "utils/Factory.h"


// Define the material 
namespace AMP {
namespace Materials {
    class RadialPolynomial : public AMP::Materials::Material { public:
    	RadialPolynomial();
    };
}
}


// Add static initialize to force symbols to be included 
// It will register the material with the factory
static struct RadialPolynomial_INIT {
    RadialPolynomial_INIT() {
        static AMP::voodoo::Registration<AMP::Materials::Material,AMP::Materials::RadialPolynomial> reg("RadialPolynomial");
     }
} RadialPolynomial_init;


#endif
