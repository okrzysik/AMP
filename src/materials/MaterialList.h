// Define the materials in AMP
#ifndef included_AMP_MaterialList
#define included_AMP_MaterialList

#include "AMP/materials/Material.h"


namespace AMP::Materials {


class CylindricallySymmetric : public AMP::Materials::Material
{
public:
    CylindricallySymmetric();
};

class Dr_nonlinear : public Material
{
public:
    Dr_nonlinear();
};

class FixedClad : public Material
{
public:
    FixedClad();
};

class FixedFuel : public Material
{
public:
    FixedFuel();
};

class Independent : public Material
{
public:
    Independent();
};

class Ox_MSRZC_09 : public Material
{
public:
    Ox_MSRZC_09();
};

class Steel316_MSRZC_09 : public Material
{
public:
    Steel316_MSRZC_09();
};

class UO2_MSRZC_09 : public Material
{
public:
    UO2_MSRZC_09();
};

class WaterLibrary : public Material
{
public:
    WaterLibrary();
};


// Register the materials
REGISTER_MATERIAL( CylindricallySymmetric );
REGISTER_MATERIAL( Dr_nonlinear );
REGISTER_MATERIAL( FixedClad );
REGISTER_MATERIAL( FixedFuel );
REGISTER_MATERIAL( Independent );
REGISTER_MATERIAL( Ox_MSRZC_09 );
REGISTER_MATERIAL( Steel316_MSRZC_09 );
REGISTER_MATERIAL( UO2_MSRZC_09 );
REGISTER_MATERIAL( WaterLibrary );


} // namespace AMP::Materials

#endif
