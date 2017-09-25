/*
 * \file materials/WaterLibrary.h
 * \brief Material library that contains all the properties of water as a coolant for subchannel
 * flow.
 */

#include "materials/WaterLibrary.h"
#include "materials/Material.h"
#include "materials/Property.h"

#include <iomanip>
#include <string>
#include <valarray>

namespace AMP {
namespace Materials {
namespace WaterLibrary_NS {

//=================== Constants =====================================================

static const std::string name_base( "WaterLibrary" );
static const std::string source( "\
M. P. Paulsen, et. al, \
RETRAN-3D, \
Electric Power Research Institute, Technical Document NP-7450, GeomType::Volume 1, \
September 1998" );

// Temperature as a function of enthalpy and pressure
static const unsigned int TempNumArgs         = 2;
static const unsigned int TempNumParams       = 2 * 4 + 5 * 5 + 5 * 5 + 5 * 5;
static const double TempParams[TempNumParams] = {
    .3276275552e2,    .9763617000e0,    .1857226027e-3,   -.4682674330e-6, // CT1
    .3360880214e-2,   -.5595281760e-4,  .1618595991e-6,   -.1180204381e-9,  .6390801208e3,
    -.3055217235e1,   .8713231868e-2,   -.6269403683e-5,  -.9844700000e-17, // CT2
    -.4302857237e0,   .2673303442e-2,   -.5198380474e-5,  .3037825558e-8,   .3309704045e-12,
    .1174524584e-3,   -.6839200986e-6,  .1168011772e-8,   -.4260074181e-12, -.2732087163e-15,
    -.1473977290e-7,  .8018858166e-10,  -.1164901355e-12, .4559400289e-17,  .5825511142e-19,
    .7104327342e-12,  -.3649500626e-14, .4457387575e-17,  .1678398723e-20,  -.3756804091e-23,
    -.1179100862e5,   .2829274345e2,    -.2678181564e-1,  .1218742752e-4,   -.2092033147e-8, // CT3
    .1256160907e3,    -.3333448850,     .3326901268e-3,   -.1477890326e-6,  .2463258371e-10,
    -.1083713369,     .2928177730e-3,   -.2972436458e-6,  .1342639113e-9,   -.2275585718e-13,
    .3278071846e-4,   -.8970959364e-7,  .9246248312e-10,  -.4249155515e-13, .7338316751e-17,
    -.3425564927e-8,  .9527692453e-11,  -.1001409043e-13, .4703914404e-17,  -.8315044742e-21,
    .3795256853e4,    -.6347031007e1,   .2867228326e-2,   .5953599813e-8,   .4798207438e-10, // CT4
    -.3910086240e1,   .1222747819e-1,   -.1404664699e-4,  .7505679464e-8,   -.1608693653e-11,
    .3410500159e-4,   .7010900113e-9,   -.1030201866e-9,  .5731099333e-14,  .3720795449e-16,
    .1527377542e-6,   -.5356866315e-9,  .6823225984e-12,  -.3668096142e-15, .6946004624e-19,
    -.1437179752e-10, .5006731336e-13,  -.6365519546e-16, .3473711350e-19,  -.6842306083e-23
};
static const std::string TempArgs[TempNumArgs] = { "enthalpy", "pressure" };
static const double TempHminVal                = 0;        // minimum enthalpy [J/kg]
static const double TempHmaxVal                = 4.0705e6; // maximum enthalpy of correlation [J/kg]
static const double TempPminVal                = 689.4757; // minimum pressure [Pa]
static const double TempPmaxVal =
    4.1369e7; // This is limiting value of correlation (6000 psi), not the critical pressure [Pa]
static const double TempRanges[2][2] = { { TempHminVal, TempHmaxVal },
                                         { TempPminVal, TempPmaxVal } };

// Saturated liquid enthalpy as a function of pressure
static const unsigned int HfSatNumArgs          = 1;
static const unsigned int HfSatNumParams        = 3 * 9;
static const double HfSatParams[HfSatNumParams] = {
    .6970887859e2,   .3337529994e2,  .2318240735e1,   .1840599513e0,   -.5245502284e-2,
    .2878007027e-2,  .1753652324e-2, -.4334859620e-3, .3325699282e-4,  .8408618802e6,
    .3637413208e6,   -.4634506669e6, .1130306339e6,   -.4350217298e3,  -.3898988188e4,
    .6697399434e3,   -.4730726377e2, .1265125057e1,   .9060030436e3,   -.1426813520e2,
    .1522233257e1,   -.6973992961e0, .1743091663e0,   -.2319717696e-1, .1694019149e-2,
    -.6454771710e-4, .1003003098e-5
};
static const std::string HfSatArgs[HfSatNumArgs] = { "pressure" };
static const double HfSatPminVal                 = 689.4757; // minimum pressure [Pa]
static const double HfSatPmaxVal      = 22119759.4074; // critical pressure; maximum pressure [Pa]
static const double HfSatRanges[1][2] = { { HfSatPminVal, HfSatPmaxVal } };

// Saturated vapor enthalpy as a function of pressure
static const unsigned int HgSatNumArgs           = 1;
static const unsigned int HgSatNumParams         = 12 + 9 + 7;
static const double HgSatParams[HgSatNumParams]  = { 0.1105836875e4,
                                                    0.1436943768e2,
                                                    0.8018288621,
                                                    0.1617232913e-1,
                                                    -0.1501147505e-2,
                                                    0.0,
                                                    0.0,
                                                    0.0,
                                                    0.0,
                                                    -0.1237675562e-4,
                                                    0.3004773304e-5,
                                                    -0.2062390734e-6,
                                                    -0.2234264997e7,
                                                    0.1231247634e7,
                                                    -0.1978847871e6,
                                                    0.1859988044e2,
                                                    -0.2765701318e1,
                                                    0.1036033878e4,
                                                    -0.2143423131e3,
                                                    0.1690507762e2,
                                                    -0.4864322134,
                                                    0.9059978254e3,
                                                    0.5561957539e1,
                                                    0.3434189609e1,
                                                    -0.6406390628,
                                                    0.5918579484e-1,
                                                    -0.2725378570e-2,
                                                    0.5006336938e-4 };
static const std::string HgSatArgs[HgSatNumArgs] = { "pressure" };
static const double HgSatPminVal                 = 689.4757; // minimum pressure [Pa]
static const double HgSatPmaxVal      = 22119759.4074; // critical pressure; maximum pressure [Pa]
static const double HgSatRanges[1][2] = { { HgSatPminVal, HgSatPmaxVal } };


// Specific volume as a function of enthalpy and pressure
static const unsigned int VolNumArgs        = 2;
static const unsigned int VolNumParams      = 3 * 3 + 3 * 5 + 4 * 3 + 4 * 4 + 3 + 4 + 3 + 4 + 1;
static const double VolParams[VolNumParams] = {
    1.1956901695e-9,   3.7591804374e-11,
    -2.4473796276e-13, // CN0
    1.6173258743e-13,  -2.1383283863e-14, 9.3054844544e-17,  7.4927085737e-17,  4.2176141427e-18,
    -1.1512516405e-20, -.4117961750e1,    -.3811294543e-3,   .4308265942e-5,    -.9160120130e-8,
    .8017924673e-11, // CN1
    -.4816067020e-5,   .7744786733e-7,    -.6988467605e-9,   .1916720525e-11,   -.1760288590e-14,
    -.1820625039e-8,   .1440785930e-10,   -.2082170753e-13,  -.3603625114e-16,  .7407124321e-19,
    -0.1403086182e4,   0.1802594763e1,
    -0.2097279215e-3, // CN2
    0.3817195017,      -0.5394444747e-3,  0.1855203702e-6,   -0.6449501159e-4,  0.8437637660e-7,
    -0.2713755001e-10, 0.7823817858e-8,   -0.1053834646e-10, 0.3629590764e-14,  2.6481168e11,
    -2.652416617e11,   2.0712464649e13,
    1.7871898061e13, // CN3
    -8.219232e8,       8.219232e8,        -6.23137536e10,    -5.6473168896e10,  8.44992e5,
    -8.44992e5,        6.2332416e7,       5.9346432e7,       -2.88e2,           2.88e2,
    -2.0736e4,         -2.0736e4,         -2.3368086721,     -2.6889781574e-4,
    1.7249452574e-8,                                                           // CP
    2.584132857e2,     1.983901141e-2,    -5.818965046e-6,   7.6017789231e-10, // CX
    8.7195773415e-3,   -1.7052703162e-6,
    1.0827981153e-10,                                                           // CT
    1.36221661279,     -1.4985169728e-4,  2.7387521579e-8,   -2.9162058556e-12, // CJ
    -2.3256803936e-9
};
static const std::string VolArgs[VolNumArgs] = { "enthalpy", "pressure" };
static const double VolHminVal               = 0; // minimum enthalpy [J/kg]
static const double VolHmaxVal = 4.0705e6;        // maximum enthalpy; critical enthalpy [J/kg]
static const double VolPminVal = 689.4757;        // minimum pressure [Pa]
static const double VolPmaxVal =
    4.1369e7; // This is limiting value of correlation (6000 psi), not the critical pressure [Pa]
static const double VolRanges[2][2] = { { VolHminVal, VolHmaxVal }, { VolPminVal, VolPmaxVal } };

// Thermal conductivity as a function of temperature and density
static const unsigned int CondNumArgs          = 2;
static const unsigned int CondNumParams        = 0;
static const double CondParams[1]              = { 0 };
static const std::string CondArgs[CondNumArgs] = { "temperature", "density" };
static const double CondTminVal                = 200.0; // minimum temperature [K]
static const double CondTmaxVal =
    1500; // maximum temperature [K] (arbitrary "very high" temperature)
static const double CondRhominVal = 318; // minimum density [kg/m3]
static const double CondRhomaxVal =
    1200.; // maximum density [kg/m3] (arbitrary "very high" density)
static const double CondRanges[2][2] = { { CondTminVal, CondTmaxVal },
                                         { CondRhominVal, CondRhomaxVal } };

// Convective Heat Coefficient as a function of thermal Conductivity,
// Diameter, Reynolds Number and Prandtl Number
static const unsigned int ConvNumArgs          = 5;
static const unsigned int ConvNumParams        = 0;
static const double ConvParams[1]              = { 0 };
static const std::string ConvArgs[ConvNumArgs] = {
    "temperature", "density", "diameter", "reynolds", "prandtl"
};
static const double ConvTminVal = 200.0; // minimum temperature [K]
static const double ConvTmaxVal =
    1500; // maximum temperature [K] (arbitrary "very high" temperature)
static const double ConvRhominVal = 318; // minimum density [kg/m3]
static const double ConvRhomaxVal =
    1200.; // maximum density [kg/m3] (arbitrary "very high" density)
static const double ConvDminVal      = 0.0; // minimum diameter [m]
static const double ConvDmaxVal      = 1.0; // maximum diameter [m] (arbitrary "very high" diameter)
static const double ConvReyminVal    = 1e3; // minimum reynolds # []
static const double ConvReymaxVal    = 1e6; // maximum reynolds # []
static const double ConvPrtminVal    = 0.87; // minimum Prandtl # []
static const double ConvPrtmaxVal    = 14.;  // maximum Prandtl # []
static const double ConvRanges[5][2] = { { ConvTminVal, ConvTmaxVal },
                                         { ConvRhominVal, ConvRhomaxVal },
                                         { ConvDminVal, ConvDmaxVal },
                                         { ConvReyminVal, ConvReymaxVal },
                                         { ConvPrtminVal, ConvPrtmaxVal } };


// dynamic viscosity as a function of temperature and density
static const unsigned int ViscNumArgs         = 2;
static const unsigned int ViscNumParams       = 4 + 5 * 6;
static const double ViscParams[ViscNumParams] = {
    0.0181583,  0.0177624,  0.0105287,  -0.0036744, 0.501938,  0.162888,  -0.130356,
    0.907919,   -0.551119,  0.146543,   0.235622,   0.789393,  0.673665,  1.207552,
    0.0670665,  -0.0843370, -0.274637,  -0.743539,  -0.959456, -0.687343, -0.497089,
    0.195286,   0.145831,   0.263129,   0.347247,   0.213486,  0.100754,  -0.032932,
    -0.0270448, -0.0253093, -0.0267758, -0.0822904, 0.0602253, -0.0202595
};
static const std::string ViscArgs[ViscNumArgs] = { "temperature", "density" };
static const double ViscTminVal                = 0.0; // minimum temperature [K]
static const double ViscTmaxVal =
    1.5e3; // maximum temperature [K] (arbitrary "very high" temperature)
static const double ViscRhominVal = 318; // minimum density [kg/m3]
static const double ViscRhomaxVal =
    1200.; // maximum density [kg/m3] (arbitrary "very high" density)
static const double ViscRanges[2][2] = { { ViscTminVal, ViscTmaxVal },
                                         { ViscRhominVal, ViscRhomaxVal } };

// enthalpy as a function of temperature and pressure
static const unsigned int EnthalpyNumArgs              = 2;
static const unsigned int EnthalpyNumParams            = 6;
static const double EnthalpyParams[EnthalpyNumParams]  = { -256638.942,    -203.118982,
                                                          0.760349801,    -3848757.66,
                                                          -0.00106377488, 0.0000006177396046 };
static const std::string EnthalpyArgs[EnthalpyNumArgs] = { "temperature", "pressure" };
static const double EnthalpyTminVal                    = 0.0; // minimum temperature [K]
static const double EnthalpyTmaxVal =
    1.5e3; // maximum temperature [K] (arbitrary "very high" temperature)
static const double EnthalpyPminVal = 689.4757; // minimum pressure [Pa]
// static const double       EnthalpyPmaxVal = 22119759.4074;    // critical pressure; maximum
// pressure [Pa]
static const double EnthalpyPmaxVal =
    4.1369e7; // This is limiting value of correlation (6000 psi), not the critical pressure [Pa]
static const double EnthalpyRanges[2][2] = { { EnthalpyTminVal, EnthalpyTmaxVal },
                                             { EnthalpyPminVal, EnthalpyPmaxVal } };

//=================== Classes =======================================================

class TemperatureProp : public Property<double>
{
public:
    TemperatureProp()
        : Property<double>( name_base + "_" + "Temperature", // Name string
                            source,                          // Reference source
                            TempParams,                      // Property parameters
                            TempNumParams,                   // Number of parameters
                            TempArgs,    // Names of arguments to the eval function
                            TempNumArgs, // Number of arguments
                            TempRanges )
    {
    } // Range of variables

    double eval( std::vector<double> &args ) override;
};

class SaturatedLiquidEnthalpyProp : public Property<double>
{
public:
    SaturatedLiquidEnthalpyProp()
        : Property<double>( name_base + "_" + "SaturatedLiquidEnthalpy", // Name string
                            source,                                      // Reference source
                            HfSatParams,                                 // Property parameters
                            HfSatNumParams,                              // Number of parameters
                            HfSatArgs,    // Names of arguments to the eval function
                            HfSatNumArgs, // Number of arguments
                            HfSatRanges )
    {
    } // Range of variables

    double eval( std::vector<double> &args ) override;
};

class SaturatedVaporEnthalpyProp : public Property<double>
{
public:
    SaturatedVaporEnthalpyProp()
        : Property<double>( name_base + "_" + "SaturatedVaporEnthalpy", // Name string
                            source,                                     // Reference source
                            HgSatParams,                                // Property parameters
                            HgSatNumParams,                             // Number of parameters
                            HgSatArgs,    // Names of arguments to the eval function
                            HgSatNumArgs, // Number of arguments
                            HgSatRanges )
    {
    } // Range of variables

    double eval( std::vector<double> &args ) override;
};

class SpecificVolumeProp : public Property<double>
{
public:
    SpecificVolumeProp()
        : Property<double>( name_base + "_" + "SpecificVolume", // Name string
                            source,                             // Reference source
                            VolParams,                          // Property parameters
                            VolNumParams,                       // Number of parameters
                            VolArgs,    // Names of arguments to the eval function
                            VolNumArgs, // Number of arguments
                            VolRanges )
    {
    } // Range of variables

    double eval( std::vector<double> &args ) override;
};

class ThermalConductivityProp : public Property<double>
{
public:
    ThermalConductivityProp()
        : Property<double>( name_base + "_" + "ThermalConductivity", // Name string
                            source,                                  // Reference source
                            CondParams,                              // Property parameters
                            CondNumParams,                           // Number of parameters
                            CondArgs,    // Names of arguments to the eval function
                            CondNumArgs, // Number of arguments
                            CondRanges )
    {
    } // Range of variables

    double eval( std::vector<double> &args ) override;
};

class ConvectiveHeatProp : public Property<double>
{
public:
    ConvectiveHeatProp()
        : Property<double>( name_base + "_" + "ConvectiveHeat", // Name string
                            source,                             // Reference source
                            ConvParams,                         // Property parameters
                            ConvNumParams,                      // Number of parameters
                            ConvArgs,    // Names of arguments to the eval function
                            ConvNumArgs, // Number of arguments
                            ConvRanges )
    {
    } // Range of variables

    double eval( std::vector<double> &args ) override;
};

class DynamicViscosityProp : public Property<double>
{
public:
    DynamicViscosityProp()
        : Property<double>( name_base + "_" + "DynamicViscosity", // Name string
                            source,                               // Reference source
                            ViscParams,                           // Property parameters
                            ViscNumParams,                        // Number of parameters
                            ViscArgs,    // Names of arguments to the eval function
                            ViscNumArgs, // Number of arguments
                            ViscRanges )
    {
    } // Range of variables

    double eval( std::vector<double> &args ) override;
};

class EnthalpyProp : public Property<double>
{
public:
    EnthalpyProp()
        : Property<double>( name_base + "_" + "Enthalpy", // Name string
                            source,                       // Reference source
                            EnthalpyParams,               // Property parameters
                            EnthalpyNumParams,            // Number of parameters
                            EnthalpyArgs,                 // Names of arguments to the eval function
                            EnthalpyNumArgs,              // Number of arguments
                            EnthalpyRanges )
    {
    } // Range of variables

    double eval( std::vector<double> &args ) override;

private:
    double MfpSolve( double, double, double, double );
    double Residual( double, double, double );
    static const double Mfp_atol;          // absolute tolerance for MFP solve
    static const double Mfp_rtol;          // relative tolerance for MFP solve
    static const double Mfp_ftol;          // function tolerance for MFP solve
    static const unsigned int Mfp_maxIter; // maximum number of iterations for MFP solve
};
const double EnthalpyProp::Mfp_atol          = 1.0e-7;  // absolute tolerance for MFP solve
const double EnthalpyProp::Mfp_rtol          = 1.0e-7;  // relative tolerance for MFP solve
const double EnthalpyProp::Mfp_ftol          = 1.0e-14; // function tolerance for MFP solve
const unsigned int EnthalpyProp::Mfp_maxIter = 1000; // maximum number of iterations for MFP solve


//=================== Functions =====================================================

inline double TemperatureProp::eval( std::vector<double> &args )
{
    double H = args[0]; // local enthalpy in J/kg
    double P = args[1]; // local pressure in Pa
    double T;           // temperature in Kelvin

    // convert SI units to units used in correlation
    double PaToPsi     = 1.45037738e-4;
    double JKgToBtuLbm = 4.29922614e-4;
    P                  = P * PaToPsi;
    H                  = H * JKgToBtuLbm;

    const double P_crit = 3208.2; // critical pressure [psi]
    const double H_crit = 906;    // specific enthalpy at critical pressure [Btu/lbm]

    if ( P < 0.1 )
        AMP_ERROR( "Liquid water temperature called with pressure below 0.1 psi." );
    if ( H <= 0 )
        AMP_ERROR( "Liquid water temperature called with enthalpy below 0 Btu/lbm." );

    // extract parameters from parameter array
    std::valarray<double> Param = get_parameters();
    double ct1[2][4], ct2[5][5], ct3[5][5], ct4[5][5];
    int offset = 0;
    for ( int i = 0; i < 2; i++ )
        for ( int j = 0; j < 4; j++ )
            ct1[i][j] = Param[offset + 4 * i + j];

    offset += 4 * 2;
    for ( int i = 0; i < 5; i++ )
        for ( int j = 0; j < 5; j++ )
            ct2[i][j] = Param[offset + 5 * i + j];

    offset += 5 * 5;
    for ( int i = 0; i < 5; i++ )
        for ( int j = 0; j < 5; j++ )
            ct3[i][j] = Param[offset + 5 * i + j];

    offset += 5 * 5;
    for ( int i = 0; i < 5; i++ )
        for ( int j = 0; j < 5; j++ )
            ct4[i][j] = Param[offset + 5 * i + j];

    // calculate temperature
    T = 0;

    if ( P >= P_crit ) {
        if ( H <= H_crit ) {
            // Eq. III.1-6b
            for ( int i = 0; i < 5; i++ )
                for ( int j = 0; j < 5; j++ )
                    T = T + ct2[i][j] * pow( P, (double) i ) * pow( H, (double) j );
        } else {
            // Eq. III.1-6d
            for ( int i = 0; i < 5; i++ )
                for ( int j = 0; j < 5; j++ )
                    T = T + ct4[i][j] * pow( P, (double) i ) * pow( H, (double) j );
        }
    } else {
        // Evaluate saturated liquid/vapor properties
        SaturatedLiquidEnthalpyProp hf_obj;
        std::vector<double> PVec( 1, P / PaToPsi ); // convert back to Pa
        double Hf = hf_obj.eval( PVec );
        Hf        = Hf * JKgToBtuLbm;

        SaturatedVaporEnthalpyProp hg_obj;
        double Hg = hg_obj.eval( PVec );
        Hg        = Hg * JKgToBtuLbm;

        if ( H <= Hf ) {
            // Eq. III.1-6a
            for ( int i = 0; i < 2; i++ )
                for ( int j = 0; j < 4; j++ )
                    T = T + ct1[i][j] * pow( P, (double) i ) * pow( H, (double) j );
        } else if ( H <= Hg ) {
            // Eq. III.1-6a with H=Hf
            for ( int i = 0; i < 2; i++ )
                for ( int j = 0; j < 4; j++ )
                    T = T + ct1[i][j] * pow( P, (double) i ) * pow( Hf, (double) j );
        } else {
            // Eq. III.1-6c
            for ( int i = 0; i < 5; i++ )
                for ( int j = 0; j < 5; j++ )
                    T = T + ct3[i][j] * pow( P, (double) i ) * pow( H, (double) j );
        }
    }
    NULL_USE( H_crit );

    // convert result to SI units
    T = 5. / 9. * ( T - 32. ) + 273.15; // [F] to [K]

    return T;
}

inline double SaturatedLiquidEnthalpyProp::eval( std::vector<double> &args )
{
    double P = args[0]; // local pressure in Pa
    double Hf;          // saturated liquid enthalpy [J/kg]

    // convert SI units to units used in correlation
    double PaToPsi     = 1.45037738e-4;
    double JKgToBtuLbm = 4.29922614e-4;
    P                  = P * PaToPsi;
    double Pmin        = HfSatPminVal * PaToPsi; // [Pa] to [psi]

    const double P_crit = 3208.2; // critical pressure [psi]

    // extract parameters from parameter array
    std::valarray<double> Param = get_parameters();
    double a[9], b[9], c[9];
    for ( int i = 0; i < 9; i++ ) {
        a[i] = Param[i];
        b[i] = Param[9 + i];
        c[i] = Param[18 + i];
    }

    // calculate saturated liquid enthalpy
    Hf = 0;
    if ( P >= Pmin && P < 900 ) // liquid region
    {
        // Eq. III.1-4a
        for ( int i = 0; i < 9; i++ )
            Hf = Hf + a[i] * pow( log( P ), (double) i );
    } else if ( P >= 900 && P < 2450 ) {
        // Eq. III.1-4b
        for ( int i = 0; i < 9; i++ )
            Hf = Hf + b[i] * pow( log( P ), (double) i );
    } else if ( P >= 2450 && P < P_crit ) {
        // Eq. III.1-4c
        for ( int i = 0; i < 9; i++ )
            Hf = Hf + c[i] * pow( pow( ( P_crit - P ), 0.41 ), (double) i );
    } else {
        AMP_ERROR( "Saturated Liquid Enthalpy: Pressure out of range of correlation." );
    }

    // convert result to SI units
    Hf = Hf / JKgToBtuLbm;

    return Hf;
}

inline double SaturatedVaporEnthalpyProp::eval( std::vector<double> &args )
{
    double P = args[0]; // local pressure in Pa
    double Hf;          // saturated liquid enthalpy [J/kg]

    // convert SI units to units used in correlation
    double PaToPsi     = 1.45037738e-4;
    double JKgToBtuLbm = 4.29922614e-4;
    P                  = P * PaToPsi;
    double Pmin        = HgSatPminVal * PaToPsi;

    const double P_crit = 3208.2; // critical pressure [psi]

    // extract parameters from parameter array
    std::valarray<double> Param = get_parameters();
    double a[12], b[9], c[7];
    for ( int i = 0; i < 12; i++ ) {
        a[i] = Param[i];
    }
    for ( int i = 0; i < 9; i++ ) {
        b[i] = Param[12 + i];
    }
    for ( int i = 0; i < 7; i++ ) {
        c[i] = Param[21 + i];
    }

    // calculate saturated vapor enthalpy
    Hf = 0;
    if ( P >= Pmin && P < 1300 ) // liquid region
    {
        // Eq. III.1-5a
        for ( int i = 0; i < 12; i++ )
            Hf = Hf + a[i] * pow( log( P ), (double) i );
    } else if ( P >= 1300 && P < 2600 ) {
        // Eq. III.1-5b
        for ( int i = 0; i < 9; i++ )
            Hf = Hf + b[i] * pow( log( P ), (double) i );
    } else if ( P >= 2600 && P < P_crit ) {
        // Eq. III.1-5c
        for ( int i = 0; i < 7; i++ )
            Hf = Hf + c[i] * pow( pow( ( P_crit - P ), 0.41 ), (double) i );
    } else {
        AMP_ERROR( "Saturated Vapor Enthalpy: Pressure out of range of correlation." );
    }

    // convert result to SI units
    Hf = Hf / JKgToBtuLbm;

    return Hf;
}

inline double SpecificVolumeProp::eval( std::vector<double> &args )
{
    double H = args[0]; // local enthalpy in J/kg
    double P = args[1]; // local pressure in Pa
    double V;           // specific volume in m3/kg

    // convert SI units to units used in correlation
    double PaToPsi     = 1.45037738e-4;
    double JKgToBtuLbm = 4.29922614e-4;
    P                  = P * PaToPsi;
    H                  = H * JKgToBtuLbm;

    const double P_crit = 3208.2; // critical pressure [psi]
    const double H_crit = 906;    // specific enthalpy at critical pressure [Btu/lbm]

    if ( P < 0.1 )
        AMP_ERROR( "Liquid water specific volume called with pressure less than 0.1 psi." );
    if ( H <= 0 )
        AMP_ERROR( "Liquid water specific volume called with enthalpy at or below 0 Btu/lbm." );

    // extract parameters from parameter array
    std::valarray<double> Param = get_parameters();
    double cn0[3][3], cn1[3][5], cn2[4][3], cn3[4][4], cp[3], cx[4], ct[3], cj[4], d;
    int offset = 0;
    for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
            cn0[i][j] = Param[offset + 3 * i + j];

    offset += 3 * 3;
    for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 5; j++ )
            cn1[i][j] = Param[offset + 5 * i + j];

    offset += 3 * 5;
    for ( int i = 0; i < 4; i++ )
        for ( int j = 0; j < 3; j++ )
            cn2[i][j] = Param[offset + 3 * i + j];

    offset += 4 * 3;
    for ( int i = 0; i < 4; i++ )
        for ( int j = 0; j < 4; j++ )
            cn3[i][j] = Param[offset + 4 * i + j];

    offset += 4 * 4;
    for ( int i = 0; i < 3; i++ )
        cp[i] = Param[offset + i];

    offset += 3;
    for ( int i = 0; i < 4; i++ )
        cx[i] = Param[offset + i];

    offset += 4;
    for ( int i = 0; i < 3; i++ )
        ct[i] = Param[offset + i];

    offset += 3;
    for ( int i = 0; i < 4; i++ )
        cj[i] = Param[offset + i];

    offset += 4;
    d = Param[offset];

    // Start evaluation
    V = 0;

    // Always same correlation if H<=250
    if ( H <= 250 ) {
        // Eq. III.1-8a
        for ( int i = 0; i < 3; i++ ) {
            for ( int j = 2; j < 5; j++ ) {
                int jj = j - 2;
                V += cn0[i][jj] * pow( P, (double) i ) * pow( ( 250 - H ), (double) j );
            }
        }
        double ExpSum = 0;
        for ( int i = 0; i < 3; i++ ) {
            for ( int j = 0; j < 5; j++ ) {
                ExpSum += cn1[i][j] * pow( P, (double) i ) * pow( H, (double) j );
            }
        }
        V += exp( ExpSum );
    }
    // If above critical pressure, make sure to avoid evaluating
    // saturated liquid/vapor
    else if ( P >= P_crit ) {
        if ( H <= H_crit ) {
            // Eq. III.1-8b
            double ExpSum = 0;
            for ( int i = 0; i < 3; i++ ) {
                for ( int j = 0; j < 5; j++ ) {
                    ExpSum += cn1[i][j] * pow( P, (double) i ) * pow( H, (double) j );
                }
            }
            V = exp( ExpSum );
        } else if ( H <= 1050 ) {
            // Eq. III.1-8d

            // pre-compute C values
            std::vector<double> C( 4, 0.0 );
            double ExpSum = 0.0;
            for ( int i = 0; i < 3; ++i ) {
                ExpSum += cp[i] * pow( P, (double) i );
            }
            C[0] = exp( ExpSum );
            for ( int i = -1; i < 3; ++i ) {
                C[1] += cx[i + 1] * pow( P, (double) i );
            }
            for ( int i = 0; i < 3; ++i ) {
                C[2] += ct[i] * pow( P, (double) i );
            }
            C[2] *= C[0];
            for ( int i = -1; i < 3; ++i ) {
                C[3] += cj[i + 1] * pow( P, (double) i );
            }

            for ( int i = 0; i < 4; ++i ) {
                double tmpSum = 0.0;
                for ( int j = 0; j < 4; ++j ) {
                    tmpSum += cn3[i][j] * C[j];
                }
                V += pow( H, (double) i ) * d * tmpSum;
            }
        } else {
            // Eq. III.1-8c
            for ( int i = -1; i < 3; ++i ) {
                for ( int j = 0; j < 3; ++j ) {
                    V += cn2[i + 1][j] * pow( P, (double) i ) * pow( H, (double) j );
                }
            }
        }
    } else {
        // Evaluate saturated liquid/vapor properties
        SaturatedLiquidEnthalpyProp hf_obj;
        std::vector<double> PVec( 1, P / PaToPsi ); // convert back to Pa
        double Hf = hf_obj.eval( PVec );
        Hf        = Hf * JKgToBtuLbm;

        SaturatedVaporEnthalpyProp hg_obj;
        double Hg = hg_obj.eval( PVec );
        Hg        = Hg * JKgToBtuLbm;

        // Liquid
        if ( H <= Hf ) {
            // Eq. III.1-8b
            double ExpSum = 0.0;
            for ( int i = 0; i < 3; i++ ) {
                for ( int j = 0; j < 5; j++ ) {
                    ExpSum += cn1[i][j] * pow( P, (double) i ) * pow( H, (double) j );
                }
            }
            V = exp( ExpSum );
        }
        // Two-phase, compute based on Hf and Hg
        else if ( H <= Hg ) {
            // Vapor quality
            double X = ( H - Hf ) / ( Hg - Hf );

            // Compute Vf from Eq. III.1-8b and Vg from III.1-8c
            double Vf = 0.0, Vg = 0.0;
            double ExpSum = 0.0;
            for ( int i = 0; i < 3; i++ ) {
                for ( int j = 0; j < 5; j++ ) {
                    ExpSum += cn1[i][j] * pow( P, (double) i ) * pow( Hf, (double) j );
                }
            }
            Vf = exp( ExpSum );
            for ( int i = -1; i < 3; ++i ) {
                for ( int j = 0; j < 3; ++j ) {
                    Vg += cn2[i + 1][j] * pow( P, (double) i ) * pow( Hg, (double) j );
                }
            }
            V = Vf + X * ( Vg - Vf );
        }
        // Small region above vapor curve and below H=1050
        else if ( H <= 1050 ) {
            // Eq. III.1-8d

            // pre-compute C values
            std::vector<double> C( 4, 0.0 );
            double ExpSum = 0.0;
            for ( int i = 0; i < 3; ++i ) {
                ExpSum += cp[i] * pow( P, (double) i );
            }
            C[0] = exp( ExpSum );
            for ( int i = -1; i < 3; ++i ) {
                C[1] += cx[i + 1] * pow( P, (double) i );
            }
            for ( int i = 0; i < 3; ++i ) {
                C[2] += ct[i] * pow( P, (double) i );
            }
            C[2] *= C[0];
            for ( int i = -1; i < 3; ++i ) {
                C[3] += cj[i + 1] * pow( P, (double) i );
            }

            for ( int i = 0; i < 4; ++i ) {
                double tmpSum = 0.0;
                for ( int j = 0; j < 4; ++j ) {
                    tmpSum += cn3[i][j] * C[j];
                }
                V += pow( H, (double) i ) * d * tmpSum;
            }
        }
        // Vapor
        else {
            // Eq. III.1-8c
            for ( int i = -1; i < 3; ++i ) {
                for ( int j = 0; j < 3; ++j ) {
                    V += cn2[i + 1][j] * pow( P, (double) i ) * pow( H, (double) j );
                }
            }
        }
    }
    NULL_USE( H_crit );

    // convert result to SI units
    V = V * 6.24279605761446e-2; // [ft3/lbm] to [m3/kg]

    return V;
}

inline double ThermalConductivityProp::eval( std::vector<double> &args )
{
    double T   = args[0]; // local temperature in Kelvin
    double rho = args[1]; // local density in kg/m3
    double k;             // thermal conductivity in W/(K-m)

    // check bounds of inputs
    if ( rho <= 0 )
        AMP_ERROR( "Thermal ocnductivity called with density <= 0 kg/m3." );
    if ( T <= 0 )
        AMP_ERROR( "Thermal conductivity called with temperature <= 0 K." );

    // declare parameters
    double Tstar   = 647.3;                                               // [K]
    double a[4]    = { 1.02811e-2, 2.99621e-2, 1.56146e-2, -4.22464e-3 }; // [W/K-m]
    double b0      = -3.97070e-1;                                         // [W/K-m]
    double b1      = 4.00302e-1;                                          // [W/K-m]
    double b2      = 1.06000;                                             // [W/K-m]
    double B1      = -1.71587e-1;
    double B2      = 2.39219;
    double rhostar = 317.7;      // [kg/m3]
    double d1      = 7.01309e-2; // [W/K-m]
    double d2      = 1.18520e-2; // [W/K-m]
    double d3      = 1.69937e-3; // [W/K-m]
    double d4      = -1.02000;   // [W/K-m]
    double C1      = 6.42857e-1;
    double C2      = -4.11717;
    double C3      = -6.17937;
    double C4      = 3.08976e-3;
    double C5      = 8.22994e-2;

    // calculate temperature
    double Tratiosum = 0.0;
    for ( size_t i = 0; i < 4; i++ )
        Tratiosum = Tratiosum + a[i] * pow( ( T / Tstar ), (double) i );
    double k0 = pow( T / Tstar, 0.5 ) * Tratiosum;

    double kbar = b0 + b1 * ( rho / rhostar ) + b2 * exp( B1 * pow( rho / rhostar + B2, 2 ) );

    double dT = fabs( T / Tstar - 1.0 ) + C4;
    double Q  = 2.0 + C5 * pow( dT, -0.6 );
    double R  = Q + 1.0;
    double S;
    if ( T / Tstar >= 1.0 )
        S = pow( dT, -1 );
    else {
        double C6 = 1.00932e1;
        S         = C6 * pow( dT, -0.6 );
    }
    double dk1 = ( d1 * pow( Tstar / T, 10 ) + d2 ) * pow( rho / rhostar, 1.8 ) *
                 exp( C1 * ( 1.0 - pow( rho / rhostar, 2.8 ) ) );
    double dk2 = d3 * S * pow( rho / rhostar, Q ) *
                 exp( Q / R * ( 1.0 - pow( rho / rhostar, (double) R ) ) );
    double dk3 = d4 * exp( C2 * pow( T / Tstar, 1.5 ) + C3 * pow( rhostar / rho, 5.0 ) );
    double dk  = dk1 + dk2 + dk3;

    k = k0 + kbar + dk;

    return k;
}

// Compute the turbulent heat transfer coefficient
inline double ConvectiveHeatProp::eval( std::vector<double> &args )
{
    double T   = args[0]; // local temperature in Kelvin
    double rho = args[1]; // local density in kg/m3
    double D   = args[2]; // hydraulic diameter in m
    double rey = args[3]; // reynolds number
    double prt = args[4]; // prandtl numner
    double k;             // thermal conductivity in W/(K-m)
    double h;             // Convective heat transfer Coefficient in W/(K-m2)

    // check bounds of inputs
    if ( rho <= 0 )
        AMP_ERROR( "Convective Heat called with density <= 0 kg/m3." );
    if ( T <= 0 )
        AMP_ERROR( "Convective Heat called with temperature <= 0 K." );
    if ( D <= 0 )
        AMP_ERROR( "Convective Heat called with hydaulic diameter <= 0 m." );
    if ( rey <= 0 )
        AMP_ERROR( "Convective Heat called with reynolds # <= 0." );

    // get thermal conductivity
    ThermalConductivityProp tCond;
    std::vector<double> largs( 2 );
    largs[0] = T;
    largs[1] = rho;
    k        = tCond.eval( largs );

    // Get the Nusselt number
    double Nu = 0.023 * pow( rey, 0.8 ) * pow( prt, 0.4 ); // Dittus-Boelter correlation
    if ( Nu < 8.0 )
        AMP_WARNING( "Convective Heat should take into account laminar heat transfer" );

    // Compute the heat transfer coefficient for flow
    h = Nu * k / D;

    // If the Nusselt number is small, we should account for thermal conduction
    // This may require using a different distance than the hydraulic diameter D/2
    if ( Nu < 4.0 )
        AMP_WARNING( "Convective Heat should take into account conduction" );

    return h;
}

inline double DynamicViscosityProp::eval( std::vector<double> &args )
{
    double T   = args[0]; // local temperature in Kelvin
    double rho = args[1]; // local density in kg/m3
    double u;             // dynamic viscosity in Pa-s

    // check bounds of inputs
    if ( rho <= 0 )
        AMP_ERROR( "Dynamic viscosity called with density <= 0 kg/m3." );
    if ( T <= 0 )
        AMP_ERROR( "Dynamic viscosity called with temperature <= 0 K." );

    // extract parameters from parameter array
    std::valarray<double> Param = get_parameters();
    double a[4], b[5][6];
    for ( size_t i = 0; i < 4; i++ )
        a[i] = Param[i];
    for ( size_t i = 0; i < 5; i++ )
        for ( size_t j = 0; j < 6; j++ )
            b[i][j] = Param[4 + 6 * i + j];

    double Tstar   = 647.27;  // [K]
    double rhostar = 317.763; // [kg/m3]

    double sum = 0;
    for ( size_t k = 0; k < 4; k++ )
        sum = sum + a[k] * pow( Tstar / T, (double) k );
    double u0 = 1e-6 * pow( T / Tstar, 0.5 ) * pow( sum, -1.0 );

    double expsum = 0;
    for ( size_t i = 0; i < 5; i++ )
        for ( size_t j = 0; j < 6; j++ )
            expsum = expsum + b[i][j] * pow( Tstar / T - 1, (double) j ) *
                                  pow( rho / rhostar - 1, (double) i );
    u = u0 * exp( rho / rhostar * expsum );

    // According to the source of the correlation, it seems that the viscosity
    // is currently in micro Pa-s rather than Pa-s.  If this correlation is
    // to be used, this value should be verified.

    return u;
}

inline double EnthalpyProp::eval( std::vector<double> &args )
{
    double T = args[0]; // local temperature in Kelvin
    double P = args[1]; // local pressure in Pa
    double h;           // specific enthalpy in J/kg

    double hmin = 0.1;         // Can't be zero
    double hmax = TempHmaxVal; // Max value we can evaluate the Temperature
    // We used to do this solve using Newton's method, however, there
    //  were occasional instabilities when evaluating at large temperatures.
    // Switching to the Method of False Position for improved robustness,
    //  it appears to use approximately twice as many function evaluations
    //  as Newton on average.
    h = MfpSolve( hmin, hmax, T, P );
    return h;
}

inline double EnthalpyProp::Residual( double h, double T, double P )
{
    TemperatureProp temperatureProperty;
    std::vector<double> tempArgs;
    tempArgs.push_back( h );
    tempArgs.push_back( P );
    double tempResult = temperatureProperty.eval( tempArgs );

    return ( T - tempResult );
}

// Use Method of False Position (MFP) to solve for enthalpy
inline double EnthalpyProp::MfpSolve( double hmin, double hmax, double T, double P )
{
    // enthalpy values at left (l), right (r), and middle (m) as well as
    //  residual evaluated at these points
    double l, r, fl, fr;
    l  = hmin;
    r  = hmax;
    fl = Residual( hmin, T, P );
    fr = Residual( hmax, T, P );

    AMP_ASSERT( fl * fr < 0.0 ); // Must have initial bounding box
    unsigned int iter = 0;
    while ( ( r - l ) > Mfp_atol && ( r - l ) / r > Mfp_rtol ) {
        iter++;
        double m  = r - ( fr * ( r - l ) / ( fr - fl ) );
        double fm = Residual( m, T, P );
        // If fm and fl have opposite sides, the root must lie in the left
        //  "half" the new range is [l,m]
        if ( fm * fl < 0.0 ) {
            r  = m;
            fr = fm;
        }
        // Otherwise the root lies in the right "half" and the new range is [m,r]
        else {
            l  = m;
            fl = fm;
        }
        if ( std::abs( fm ) < Mfp_ftol )
            break;
        if ( iter >= Mfp_maxIter ) {
            AMP_ERROR( "MFP solve failed to converge for property function evaluation." );
        }
    }

    // Return point with smallest function value
    return std::abs( fl ) < std::abs( fr ) ? l : r;
}
} // namespace WaterLibrary_NS
//=================== Materials =====================================================

WaterLibrary::WaterLibrary()
{
    d_propertyMap = new std::map<std::string, AMP::shared_ptr<Property<double>>>();
    INSERT_PROPERTY_IN_MAP( Temperature, WaterLibrary_NS );
    INSERT_PROPERTY_IN_MAP( SaturatedLiquidEnthalpy, WaterLibrary_NS );
    INSERT_PROPERTY_IN_MAP( SaturatedVaporEnthalpy, WaterLibrary_NS );
    INSERT_PROPERTY_IN_MAP( SpecificVolume, WaterLibrary_NS );
    INSERT_PROPERTY_IN_MAP( ThermalConductivity, WaterLibrary_NS );
    INSERT_PROPERTY_IN_MAP( DynamicViscosity, WaterLibrary_NS );
    INSERT_PROPERTY_IN_MAP( Enthalpy, WaterLibrary_NS );
    INSERT_PROPERTY_IN_MAP( ConvectiveHeat, WaterLibrary_NS );
}
} // namespace Materials
} // namespace AMP
