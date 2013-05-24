/*
 * \file materials/WaterLibrary.h
 * \brief Material library that contains all the properties of water as a coolant for subchannel flow.
 */

#include "materials/WaterLibrary.h"
#include "materials/Property.h"
#include "materials/Material.h"

#include <string>
#include <valarray>

namespace AMP { 
namespace Materials {
namespace WaterLibrary_NS {

//=================== Constants =====================================================

    static const std::string name_base("WaterLibrary");
    static const std::string source("\
M. P. Paulsen, et. al, \
RETRAN-3D, \
Electric Power Research Institute, Technical Document NP-7450, Volume 1, \
September 1998");

      // Temperature as a function of enthalpy and pressure
    static const unsigned int TempNumArgs   = 2;
    static const unsigned int TempNumParams = 2*4 + 5*5;
    static const double       TempParams[TempNumParams] = {.3276275552e2,   .9763617000e0,   .1857226027e-3, -.4682674330e-6,
                                 .3360880214e-2, -.5595281760e-4,  .1618595991e-6, -.1180204381e-9,
                                 .6390801208e3,  -.3055217235e1,   .8713231868e-2, -.6269403683e-5, -.9844700000e-17,
                                -.4302857237e0,   .2673303442e-2, -.5198380474e-5,  .3037825558e-8,  .3309704045e-12,
                                 .1174524584e-3, -.6839200986e-6,  .1168011772e-8, -.4260074181e-12,-.2732087163e-15,
                                -.1473977290e-7,  .8018858166e-10,-.1164901355e-12, .4559400289e-17, .5825511142e-19,
                                 .7104327342e-12,-.3649500626e-14, .4457387575e-17, .1678398723e-20,-.3756804091e-23};
    static const std::string  TempArgs[TempNumArgs] = {"enthalpy", "pressure"};
    static const double       TempHminVal = 0;  // minimum enthalpy [J/kg]
    static const double       TempHmaxVal = 2107356; // maximum enthalpy; critical enthalpy [J/kg] 
    static const double       TempPminVal = 689.4757; // minimum pressure [Pa] 
    static const double       TempPmaxVal = 22119759.4074; // critical pressure; maximum pressure [Pa]
    static const double       TempRanges[2][2]={{TempHminVal, TempHmaxVal}, {TempPminVal, TempPmaxVal}};

      // Saturated liquid enthalpy as a function of pressure
    static const unsigned int HfSatNumArgs   = 1;
    static const unsigned int HfSatNumParams = 3*9;
    static const double       HfSatParams[HfSatNumParams] = {.6970887859e2, .3337529994e2,  .2318240735e1, .1840599513e0,-.5245502284e-2,
                                 .2878007027e-2,.1753652324e-2,-.4334859620e-3,.3325699282e-4,
                                 .8408618802e6, .3637413208e6, -.4634506669e6, .1130306339e6,-.4350217298e3,
                                -.3898988188e4, .6697399434e3, -.4730726377e2, .1265125057e1,
                                 .9060030436e3,-.1426813520e2,  .1522233257e1,-.6973992961e0,.1743091663e0,
                                -.2319717696e-1,.1694019149e-2,-.6454771710e-4,.1003003098e-5};
    static const std::string  HfSatArgs[HfSatNumArgs] = {"pressure"};
    static const double       HfSatPminVal = 689.4757; // minimum pressure [Pa] 
    static const double       HfSatPmaxVal = 22119759.4074; // critical pressure; maximum pressure [Pa]
    static const double       HfSatRanges[1][2]={{HfSatPminVal, HfSatPmaxVal}};


      // Specific volume as a function of enthalpy and pressure
    static const unsigned int VolNumArgs   = 2;
    static const unsigned int VolNumParams = 3*3 + 3*5;
    static const double       VolParams[VolNumParams] = {1.1956901695e-9,  3.7591804374e-11,-2.4473796276e-13,
                                 1.6173258743e-13,-2.1383283863e-14, 9.3054844544e-17,
                                 7.4927085737e-17, 4.2176141427e-18,-1.1512516405e-20,
                            -.4117961750e1,-.3811294543e-3,  .4308265942e-5, -.9160120130e-8,  .8017924673e-11,
                            -.4816067020e-5,.7744786733e-7, -.6988467605e-9,  .1916720525e-11,-.1760288590e-14,
                            -.1820625039e-8,.1440785930e-10,-.2082170753e-13,-.3603625114e-16, .7407124321e-19};
    static const std::string  VolArgs[VolNumArgs] = {"enthalpy", "pressure"};
    static const double       VolHminVal = 0;  // minimum enthalpy [J/kg]
    static const double       VolHmaxVal = 2107356; // maximum enthalpy; critical enthalpy [J/kg] 
    static const double       VolPminVal = 689.4757; // minimum pressure [Pa] 
    static const double       VolPmaxVal = 22119759.4074; // critical pressure; maximum pressure [Pa]
    static const double       VolRanges[2][2]={{VolHminVal, VolHmaxVal}, {VolPminVal, VolPmaxVal}};

      // Thermal conductivity as a function of temperature and density
    static const unsigned int CondNumArgs   = 2;
    static const unsigned int CondNumParams = 0;
    static const double       CondParams[1] = {0};
    static const std::string  CondArgs[CondNumArgs] = {"temperature", "density"};
    static const double       CondTminVal = 200.0;        // minimum temperature [K]
    static const double       CondTmaxVal = 800;        // maximum temperature [K] (arbitrary "very high" temperature)
    static const double       CondRhominVal = 500;        // minimum density [kg/m3] 
    static const double       CondRhomaxVal = 1200.;    // maximum density [kg/m3] (arbitrary "very high" density)
    static const double       CondRanges[2][2]={{CondTminVal, CondTmaxVal}, {CondRhominVal, CondRhomaxVal}};

      // Convective Heat Coefficient as a function of thermal Conductivity,
        // Diameter, Reynolds Number and Prandtl Number 
    static const unsigned int ConvNumArgs   = 5;
    static const unsigned int ConvNumParams = 0;
    static const double       ConvParams[1] = {0};
    static const std::string  ConvArgs[ConvNumArgs] = {"temperature", "density", "diameter", "reynolds", "prandtl"};
    static const double       ConvTminVal = 200.0;        // minimum temperature [K]
    static const double       ConvTmaxVal = 800;        // maximum temperature [K] (arbitrary "very high" temperature)
    static const double       ConvRhominVal = 500;        // minimum density [kg/m3] 
    static const double       ConvRhomaxVal = 1200.;    // maximum density [kg/m3] (arbitrary "very high" density)
    static const double       ConvDminVal   = 0.0;        // minimum diameter [m]
    static const double       ConvDmaxVal   = 1.0;        // maximum diameter [m] (arbitrary "very high" diameter)
    static const double       ConvReyminVal = 1e3;        // minimum reynolds # [] 
    static const double       ConvReymaxVal = 1e6;      // maximum reynolds # []
    static const double       ConvPrtminVal = 0.87;        // minimum Prandtl # [] 
    static const double       ConvPrtmaxVal = 14.;      // maximum Prandtl # [] 
    static const double       ConvRanges[5][2]={ {ConvTminVal, ConvTmaxVal}, {ConvRhominVal, ConvRhomaxVal}, {ConvDminVal ,ConvDmaxVal }, {ConvReyminVal,ConvReymaxVal}, {ConvPrtminVal,ConvPrtmaxVal}};


      // dynamic viscosity as a function of temperature and density
    static const unsigned int ViscNumArgs   = 2;
    static const unsigned int ViscNumParams = 4 + 5*6;
    static const double       ViscParams[ViscNumParams] = {0.0181583,0.0177624,0.0105287,-0.0036744,
                                 0.501938,  0.162888, -0.130356,  0.907919,-0.551119,  0.146543,
                                 0.235622,  0.789393,  0.673665,  1.207552, 0.0670665,-0.0843370,
                                -0.274637, -0.743539, -0.959456, -0.687343,-0.497089,  0.195286,
                                 0.145831,  0.263129,  0.347247,  0.213486, 0.100754, -0.032932,
                                -0.0270448,-0.0253093,-0.0267758,-0.0822904,0.0602253,-0.0202595};
    static const std::string  ViscArgs[ViscNumArgs] = {"temperature", "density"};
    static const double       ViscTminVal = 0.0;        // minimum temperature [K]
    static const double       ViscTmaxVal = 1.0e3;        // maximum temperature [K] (arbitrary "very high" temperature)
    static const double       ViscRhominVal = 500;        // minimum density [kg/m3] 
    static const double       ViscRhomaxVal = 1200.;    // maximum density [kg/m3] (arbitrary "very high" density)
    static const double       ViscRanges[2][2]={{ViscTminVal, ViscTmaxVal}, {ViscRhominVal, ViscRhomaxVal}};

      // enthalpy as a function of temperature and pressure
    static const unsigned int EnthalpyNumArgs   = 2;
    static const unsigned int EnthalpyNumParams = 6;
    static const double       EnthalpyParams[EnthalpyNumParams] = {-256638.942,-203.118982,0.760349801,-3848757.66,
                                    -0.00106377488,0.0000006177396046};
    static const std::string  EnthalpyArgs[EnthalpyNumArgs] = {"temperature", "pressure"};
    static const double       EnthalpyTminVal = 0.0;        // minimum temperature [K]
    static const double       EnthalpyTmaxVal = 1.0e3;        // maximum temperature [K] (arbitrary "very high" temperature)
    static const double       EnthalpyPminVal = 689.4757;        // minimum pressure [Pa] 
    static const double       EnthalpyPmaxVal = 22119759.4074;    // critical pressure; maximum pressure [Pa]
    static const double       EnthalpyRanges[2][2]={{EnthalpyTminVal, EnthalpyTmaxVal}, {EnthalpyPminVal, EnthalpyPmaxVal}};

//=================== Classes =======================================================

    class TemperatureProp : public Property<double>{
    public:
        TemperatureProp() :
            Property<double> (    name_base + "_" + "Temperature",    // Name string
                source,        // Reference source
                TempParams,    // Property parameters
                TempNumParams,  // Number of parameters
                TempArgs,      // Names of arguments to the eval function
                TempNumArgs,    // Number of arguments
                TempRanges ){}    // Range of variables

        virtual double eval( std::vector<double>& args );
    };

    class SaturatedLiquidEnthalpyProp : public Property<double>{
    public:
        SaturatedLiquidEnthalpyProp() :
            Property<double> (    name_base + "_" + "SaturatedLiquidEnthalpy",    // Name string
                source,        // Reference source
                HfSatParams,    // Property parameters
                HfSatNumParams, // Number of parameters
                HfSatArgs,      // Names of arguments to the eval function
                HfSatNumArgs,    // Number of arguments
                HfSatRanges ){}    // Range of variables

        virtual double eval( std::vector<double>& args );
    };

    class SpecificVolumeProp : public Property<double>{
    public:
        SpecificVolumeProp() :
            Property<double> (    name_base + "_" + "SpecificVolume",    // Name string
                source,        // Reference source
                VolParams,    // Property parameters
                VolNumParams,    // Number of parameters
                VolArgs,      // Names of arguments to the eval function
                VolNumArgs,    // Number of arguments
                VolRanges ){}    // Range of variables

        virtual double eval( std::vector<double>& args );
    };

    class ThermalConductivityProp : public Property<double>{
    public:
        ThermalConductivityProp() :
            Property<double> (    name_base + "_" + "ThermalConductivity",    // Name string
                source,        // Reference source
                CondParams,    // Property parameters
                CondNumParams,    // Number of parameters
                CondArgs,      // Names of arguments to the eval function
                CondNumArgs,    // Number of arguments
                CondRanges ){}    // Range of variables

        virtual double eval( std::vector<double>& args );
    };

    class ConvectiveHeatProp : public Property<double>{
    public:
        ConvectiveHeatProp() :
            Property<double> (    name_base + "_" + "ConvectiveHeat",    // Name string
                source,        // Reference source
                ConvParams,    // Property parameters
                ConvNumParams,    // Number of parameters
                ConvArgs,      // Names of arguments to the eval function
                ConvNumArgs,    // Number of arguments
                ConvRanges ){}    // Range of variables

        virtual double eval( std::vector<double>& args );
    };

    class DynamicViscosityProp : public Property<double>{
    public:
        DynamicViscosityProp() :
            Property<double> (    name_base + "_" + "DynamicViscosity",    // Name string
                source,        // Reference source
                ViscParams,    // Property parameters
                ViscNumParams,    // Number of parameters
                ViscArgs,      // Names of arguments to the eval function
                ViscNumArgs,    // Number of arguments
                ViscRanges ){}    // Range of variables

        virtual double eval( std::vector<double>& args );
    };

    class EnthalpyProp : public Property<double>{
    public:
        EnthalpyProp() :
            Property<double> (    name_base + "_" + "Enthalpy",    // Name string
                source,            // Reference source
                EnthalpyParams,        // Property parameters
                EnthalpyNumParams,    // Number of parameters
                EnthalpyArgs,          // Names of arguments to the eval function
                EnthalpyNumArgs,    // Number of arguments
                EnthalpyRanges ){}    // Range of variables

        virtual double eval( std::vector<double>& args );
        double NewtonSolve(double,double,double);
    private:
        double Residual(double,double,double);
        static const double Newton_atol; // absolute tolerance for Newton solve
        static const double Newton_rtol; // relative tolerance for Newton solve
        static const unsigned int Newton_maxIter; // maximum number of iterations for Newton solve
        static const double machinePrecision; // machine precision; used in perturbation for numerical Jacobian
    };
    const double EnthalpyProp::Newton_atol = 1.0e-7; // absolute tolerance for Newton solve
    const double EnthalpyProp::Newton_rtol = 1.0e-7; // relative tolerance for Newton solve
    const unsigned int EnthalpyProp::Newton_maxIter = 1000; // maximum number of iterations for Newton solve
    const double EnthalpyProp::machinePrecision = 1.0e-15; // machine precision; used in perturbation for numerical Jacobian


//=================== Functions =====================================================

    inline double TemperatureProp::eval( std::vector<double>& args ){
            double H            = args[0];  // local enthalpy in J/kg
            double P            = args[1];  // local pressure in Pa
            double T;                       // temperature in Kelvin
        
        // get saturated liquid enthalpy
        SaturatedLiquidEnthalpyProp hf_obj;
                std::vector<double> PVec(1, P);
                double Hf = hf_obj.eval(PVec);

        // convert SI units to units used in correlation
        P = P*1.45037738e-4;    // [Pa] to [psi]
        H = H*4.29922614e-4;    // [J/kg] to [Btu/lbm]
        Hf = Hf*4.29922614e-4;    // [J/kg] to [Btu/lbm]
           
        double P_crit = 3208.2; // critical pressure [psi]
        double H_crit = 906; // specific enthalpy at critical pressure [Btu/lbm]

        // determine if water is in either liquid or critical state
        bool InLiquidRegion = false;        // liquid region
        bool InCriticalRegion = false;        // critical region
        bool InAcceptableRegion = false;    // liquid or critical region
        if (P < P_crit && H <= Hf) InLiquidRegion = true;
        if (P >= P_crit && H <= H_crit) InCriticalRegion = true;
        InAcceptableRegion = InLiquidRegion || InCriticalRegion;
        if (!InAcceptableRegion) AMP_ERROR("Liquid water temperature called, but water not in liquid or critical state.");
        if (P < 0.1) AMP_ERROR("Liquid water temperature called with pressure below 0.1 psi.");
        if (H <= 0) AMP_ERROR("Liquid water temperature called with enthalpy below 0 Btu/lbm.");
    
        // extract parameters from parameter array
        std::valarray<double> Param = get_parameters();
        double a[2][4], b[5][5];
        for (size_t i=0; i<2; i++)
            for (size_t j=0; j<4; j++)
                a[i][j] = Param[4*i+j];
    
        for (size_t i=0; i<5; i++)
            for (size_t j=0; j<5; j++)
                b[i][j] = Param[2*4+5*i+j];
    
        // calculate temperature
        T = 0;
        if (InLiquidRegion) // liquid region
        {
            for (size_t i=0; i<2; i++)
                for (size_t j=0; j<4; j++)
                    T = T + a[i][j]*pow(P,(double)i)*pow(H,(double)j);
        }
        else if (InCriticalRegion) // critical region
        {
            for (size_t i=0; i<5; i++)
                for (size_t j=0; j<5; j++)
                    T = T + b[i][j]*pow(P,(double)i)*pow(H,(double)j);
        }

        // convert result to SI units
        T = 5./9.*(T - 32.) + 273.15; // [F] to [K]

        return T;
    }

    inline double SaturatedLiquidEnthalpyProp::eval( std::vector<double>& args ){
            double P            = args[0];  // local pressure in Pa
            double Hf;                       // saturated liquid enthalpy [J/kg]
        
        // convert SI units to units used in correlation
        P = P*1.45037738e-4;                // [Pa] to [psi]
        double Pmin = HfSatPminVal*1.45037738e-4;    // [Pa] to [psi]
           
        double P_crit = 3208.2; // critical pressure [psi]

        // extract parameters from parameter array
        std::valarray<double> Param = get_parameters();
        double a[9], b[9], c[9];
        for (size_t i=0; i<9; i++)
        {
            a[i] = Param[i];
            b[i] = Param[9+i];
            c[i] = Param[18+i];
        }
    
        // calculate saturated liquid enthalpy
        Hf = 0;
        if (P >= Pmin && P < 900) // liquid region
        {
            for (size_t i=0; i<9; i++)
                Hf = Hf + a[i]*pow(log(P),(double)i);
        }
        else if (P >= 900 && P < 2450)
        {
            for (size_t i=0; i<9; i++)
                Hf = Hf + b[i]*pow(log(P),(double)i);
        }
        else if (P >= 2450 && P < P_crit)
        {
            for (size_t i=0; i<9; i++)
                Hf = Hf + c[i]*pow(pow((P_crit-P),0.41),(double)i);
        }
        else
        {
            AMP_ERROR("Saturated Liquid Enthalpy: Pressure out of range of correlation.");
        }

        // convert result to SI units
        Hf = Hf/4.29922614e-4; // [Btu/lbm] to [J/kg]

        return Hf;
    }

    inline double SpecificVolumeProp::eval( std::vector<double>& args ){
            double H            = args[0];  // local enthalpy in J/kg
            double P            = args[1];  // local pressure in Pa
            double V;                       // specific volume in m3/kg
        
        // get saturated liquid enthalpy
        SaturatedLiquidEnthalpyProp hf_obj;
                std::vector<double> PVec(1, P);
                double Hf = hf_obj.eval(PVec);

        // convert SI units to units used in correlation
        P = P*1.45037738e-4;    // [Pa] to [psi]
        H = H*4.29922614e-4;    // [J/kg] to [Btu/lbm]
        Hf = Hf*4.29922614e-4;    // [J/kg] to [Btu/lbm]
           
        double P_crit = 3208.2; // critical pressure [psi]
        double H_crit = 906; // specific enthalpy at critical pressure [Btu/lbm]

        // determine if water is in either liquid or critical state
        bool InLiquidRegion1 = false;        // liquid region 1
        bool InLiquidRegion2 = false;        // liquid region 2
        bool InCriticalRegion = false;        // critical region
        bool InAcceptableRegion = false;    // liquid 1, liquid 2, or critical region
        if (H <= 250 && H <= Hf) InLiquidRegion1 = true;
        if (H > 250 && H <= Hf && P < P_crit) InLiquidRegion2 = true;
        if (H > 250 && H <= H_crit && P > P_crit) InCriticalRegion = true;
        InAcceptableRegion = InLiquidRegion1 || InLiquidRegion2 || InCriticalRegion;
        if (!InAcceptableRegion) AMP_ERROR("Liquid water specific volume called, but water not in liquid or critical state.");
        if (P < 0.1) AMP_ERROR("Liquid water specific volume called with pressure less than 0.1 psi.");
        if (H <= 0) AMP_ERROR("Liquid water specific volume called with enthalpy at or below 0 Btu/lbm.");
    
        // extract parameters from parameter array
        std::valarray<double> Param = get_parameters();
        double a[3][3], b[3][5];
        for (size_t i=0; i<3; i++)
            for (size_t j=0; j<3; j++)
                a[i][j] = Param[3*i+j];
    
        for (size_t i=0; i<3; i++)
            for (size_t j=0; j<5; j++)
                b[i][j] = Param[3*3+5*i+j];
    
        // calculate specific volume
        V = 0;
        if (InLiquidRegion1) // liquid region 1
        {
            for (size_t i=0; i<3; i++)
            {
                for (size_t j=2; j<5; j++)
                {
                    size_t jj = j-2;
                    V = V + a[i][jj]*pow(P,(double)i)*pow((250-H),(double)j);
                }
            }
            double ExpSum = 0;
            for (size_t i=0; i<3; i++)
            {
                for (size_t j=0; j<5; j++)
              {
                    ExpSum = ExpSum + b[i][j]*pow(P,(double)i)*pow(H,(double)j); 
                }
            }
            V = V + exp(ExpSum);
        }
        else if (InLiquidRegion2 || InCriticalRegion) // liquid region 2 or critical region
        {
            double ExpSum = 0;
            for (size_t i=0; i<3; i++)
            {
                for (size_t j=0; j<5; j++)
              {
                    ExpSum = ExpSum + b[i][j]*pow(P,(double)i)*pow(H,(double)j); 
                }
            }
            V = V + exp(ExpSum);
        }
        // convert result to SI units
        V = V*6.24279605761446e-2; // [ft3/lbm] to [m3/kg]

        return V;
    }

    inline double ThermalConductivityProp::eval( std::vector<double>& args ){
            double T            = args[0];  // local temperature in Kelvin
            double rho          = args[1];  // local density in kg/m3
            double k;                       // thermal conductivity in W/(K-m)
        
        // check bounds of inputs
        if (rho <= 0) AMP_ERROR("Thermal ocnductivity called with density <= 0 kg/m3.");
        if (T <= 0) AMP_ERROR("Thermal conductivity called with temperature <= 0 K.");
    
        // declare parameters
        double Tstar = 647.3; // [K]
        double a[4] = {1.02811e-2, 2.99621e-2, 1.56146e-2, -4.22464e-3}; // [W/K-m]
        double b0 = -3.97070e-1; // [W/K-m]
        double b1 = 4.00302e-1; // [W/K-m]
        double b2 = 1.06000; // [W/K-m]
        double B1 = -1.71587e-1;
        double B2 = 2.39219;
        double rhostar = 317.7; // [kg/m3]
        double d1 = 7.01309e-2; // [W/K-m]
        double d2 = 1.18520e-2; // [W/K-m]
        double d3 = 1.69937e-3; // [W/K-m]
        double d4 = -1.02000; // [W/K-m]
        double C1 = 6.42857e-1;
        double C2 = -4.11717;
        double C3 = -6.17937;
        double C4 = 3.08976e-3;
        double C5 = 8.22994e-2;
    
        // calculate temperature
        double Tratiosum = 0.0;
        for (size_t i=0; i<4; i++)
            Tratiosum = Tratiosum + a[i]*pow((T/Tstar),(double)i);
        double k0 = pow(T/Tstar,0.5)*Tratiosum;

        double kbar = b0 + b1*(rho/rhostar) + b2*exp(B1*pow(rho/rhostar+B2,2));

        double dT = abs(T/Tstar - 1.0) + C4;
        double Q = 2.0 + C5*pow(dT,-0.6);
        double R = Q + 1.0;
        double S;
        if (T/Tstar >= 1.0) S = pow(dT,-1);
        else {
          double C6 = 1.00932e1;
      S = C6*pow(dT,-0.6);
    }
        double dk1 = (d1*pow(Tstar/T,10) + d2) * pow(rho/rhostar,1.8) * exp(C1*(1.0 - pow(rho/rhostar,2.8)));
        double dk2 = d3*S*pow(rho/rhostar,Q)*exp(Q/R*(1.0-pow(rho/rhostar,(double)R)));
        double dk3 = d4*exp(C2*pow(T/Tstar,1.5)+C3*pow(rhostar/rho,5.0));
        double dk = dk1 + dk2 + dk3;

        k = k0 + kbar + dk;

        return k;
    }

    // Compute the turbulent heat transfer coefficient
    inline double ConvectiveHeatProp::eval( std::vector<double>& args ){
            double T            = args[0];  // local temperature in Kelvin
            double rho          = args[1];  // local density in kg/m3
            double D            = args[2];  // hydraulic diameter in m
            double rey          = args[3];  // reynolds number 
            double prt          = args[4];  // prandtl numner 
            double k;                       // thermal conductivity in W/(K-m)
            double h;                       // Convective heat transfer Coefficient in W/(K-m2)
        
        // check bounds of inputs
        if (rho <= 0) AMP_ERROR("Convective Heat called with density <= 0 kg/m3.");
        if (T <= 0) AMP_ERROR("Convective Heat called with temperature <= 0 K.");
        if (D <= 0) AMP_ERROR("Convective Heat called with hydaulic diameter <= 0 m.");
        if (rey <= 0) AMP_ERROR("Convective Heat called with reynolds # <= 0.");
        
        // get thermal conductivity  
        ThermalConductivityProp tCond;
                std::vector<double> largs(2);
                largs[0]= T;
                largs[1]= rho;
                k = tCond.eval(largs);

        // Get the Nusselt number
        double Nu = 0.023 * pow(rey, 0.8) * pow(prt, 0.4);  // Dittus-Boelter correlation
        if ( Nu < 8.0 ) AMP_WARNING("Convective Heat should take into account laminar heat transfer");

        // Compute the heat transfer coefficient for flow
        h  =  Nu*k/D;

        // If the Nusselt number is small, we should account for thermal conduction
        // This may require using a different distance than the hydraulic diameter D/2
        if ( Nu < 4.0 ) AMP_WARNING("Convective Heat should take into account conduction");

        return h;
    }

    inline double DynamicViscosityProp::eval( std::vector<double>& args ){
            double T            = args[0];  // local temperature in Kelvin
            double rho          = args[1];  // local density in kg/m3
            double u;                       // dynamic viscosity in Pa-s
        
        // check bounds of inputs
        if (rho <= 0) AMP_ERROR("Dynamic viscosity called with density <= 0 kg/m3.");
        if (T <= 0) AMP_ERROR("Dynamic viscosity called with temperature <= 0 K.");
    
        // extract parameters from parameter array
        std::valarray<double> Param = get_parameters();
        double a[4], b[5][6];
        for (size_t i=0; i<4; i++)
            a[i] = Param[i];
        for (size_t i=0; i<5; i++)
            for (size_t j=0; j<6; j++)
                b[i][j] = Param[4+6*i+j];

        double Tstar = 647.27;        // [K]
        double rhostar = 317.763;    // [kg/m3]

        double sum = 0;
        for (size_t k=0; k<4; k++)
            sum = sum + a[k]*pow(Tstar/T,(double)k);
        double u0 = 1e-6*pow(T/Tstar,0.5)*pow(sum,-1.0);

        double expsum = 0;
        for (size_t i=0; i<5; i++)
            for (size_t j=0; j<6; j++)
                expsum = expsum + b[i][j]*pow(Tstar/T-1,(double)j)*pow(rho/rhostar-1,(double)i);
        u = u0*exp(rho/rhostar*expsum);

        return u;
    }

    inline double EnthalpyProp::eval( std::vector<double>& args ){
            double T            = args[0];  // local temperature in Kelvin
            double P            = args[1];  // local pressure in Pa
            double h;                       // specific enthalpy in J/kg

        // convert SI units to units used in correlation
        double P_brit = P*1.45037738e-4;    // [Pa] to [psi]
        double P_crit = 3208.2; // critical pressure [psi]
        double h_guess = 0.0; // enthalpy guess for Newton solve
        if (P_brit < P_crit){
            // guess some h < hf and > 0
            // get hf
            SaturatedLiquidEnthalpyProp liquidEnthalpyProperty;
            std::vector<double> liqEnthalpyArgs;
            liqEnthalpyArgs.push_back(P);
            double hf = liquidEnthalpyProperty.eval(liqEnthalpyArgs);
            // pick h_guess such that: 0 < h_guess < hf
            h_guess = hf/2.0;
        } else { // P_brit >= P_crit
            // guess some h <= h_crit
            h_guess = 700.0/4.29922614e-4; // h_crit = 906 Btu/lb -> choose 700 Btu/lb and convert to SI
        }
           
        h = NewtonSolve(h_guess, T, P);
        return h;
    }

    inline double EnthalpyProp::Residual(double h, double T, double P){
        TemperatureProp temperatureProperty;
                std::vector<double> tempArgs;
        tempArgs.push_back(h);
        tempArgs.push_back(P);
                double tempResult = temperatureProperty.eval(tempArgs);

        return (T - tempResult);
    }

    inline double EnthalpyProp::NewtonSolve(double guess, double param1, double param2)
    {
        double x_new = guess;
        double x_old = guess;
        bool converged = false;
        for (unsigned int iter=1; iter<=Newton_maxIter; ++iter){
            x_old = x_new;
                        double b_perturb = sqrt(machinePrecision);
            double perturbation = (1.0+x_new)*b_perturb;
            // numerical Jacobian with forward perturbation
            double J = (Residual(x_old+perturbation,param1,param2) - Residual(x_old,param1,param2))/perturbation;
            double dx = -1.0*Residual(x_old,param1,param2)/J;
            x_new = x_old + dx;
            // check convergence
            double abs_err = std::abs(x_new - x_old); // absolute error
            double rel_err = 0.0; // relative error
            if (x_old != 0.0){ // test to ensure no division by zero
                rel_err = std::abs((x_new - x_old)/x_old);
            }
            if ((abs_err < Newton_atol) && (rel_err < Newton_rtol)){
                converged = true;
                break;
            }
        }
        if (!converged){
            AMP_ERROR("Newton solve failed to converge for property function evaluation.");
        }
        return x_new;
    }
}
//=================== Materials =====================================================

WaterLibrary::WaterLibrary()
{
        d_propertyMap = new std::map<std::string, boost::shared_ptr<Property<double> > >();
        INSERT_PROPERTY_IN_MAP(Temperature,    WaterLibrary_NS);
        INSERT_PROPERTY_IN_MAP(SaturatedLiquidEnthalpy, WaterLibrary_NS);
        INSERT_PROPERTY_IN_MAP(SpecificVolume, WaterLibrary_NS);
        INSERT_PROPERTY_IN_MAP(ThermalConductivity, WaterLibrary_NS);
        INSERT_PROPERTY_IN_MAP(DynamicViscosity, WaterLibrary_NS);
        INSERT_PROPERTY_IN_MAP(Enthalpy, WaterLibrary_NS);
        INSERT_PROPERTY_IN_MAP(ConvectiveHeat, WaterLibrary_NS);
}


} 
}
