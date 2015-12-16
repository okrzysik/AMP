//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   materials/FixedClad.h
 * \author Aaron Phillippe
 * \brief  Header file for constant cladding properties
 */
//---------------------------------------------------------------------------//
#ifndef FixedClad_H
#define FixedClad_H

#include "Material.h"
#include "utils/Factory.h"


// Define the material
namespace AMP {
namespace Materials {
class FixedClad : public Material {
public:
    FixedClad();
};
}
}


// Add static initialize to force symbols to be included
// It will register the material with the factory
static struct FixedClad_INIT {
    FixedClad_INIT()
    {
        static AMP::voodoo::Registration<AMP::Materials::Material, AMP::Materials::FixedClad> reg(
            "FixedClad" );
    }
} FixedClad_init;


#endif
