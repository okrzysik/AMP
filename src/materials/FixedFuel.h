//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   materials/FixedFuel.h
 * \author Aaron Phillippe
 * \brief  Header file for constant fuel properties  
 */
//---------------------------------------------------------------------------//
#ifndef FixedFuel_H
#define FixedFuel_H

#include "Material.h"
#include "utils/Factory.h"


// Define the material 
namespace AMP {
namespace Materials {
    class FixedFuel : public Material { public:
    	FixedFuel();
    };
}
}


// Add static initialize to force symbols to be included 
// It will register the material with the factory
static struct FixedFuel_INIT {
    FixedFuel_INIT() {
        static AMP::voodoo::Registration<AMP::Materials::Material,AMP::Materials::FixedFuel> reg("FixedFuel");
     }
} FixedFuel_init;


#endif
