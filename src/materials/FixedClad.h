//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   materials/FixedClad.h
 * \author Aaron Phillippe
 * \brief  Header file for constant cladding properties
 */
//---------------------------------------------------------------------------//
#ifndef FixedClad_H
#define FixedClad_H

#include "AMP/utils/Factory.h"
#include "Material.h"


// Define the material
namespace AMP::Materials {
class FixedClad : public Material
{
public:
    FixedClad();
};
} // namespace AMP::Materials


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
