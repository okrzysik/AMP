/*
 * ThermalDiffusionCoefficientProp.h
 *
 *  Created on: Sep 14, 2010
 *	  Author: gad
 */

/**
 * These context-dependent classes handle the thermal diffusion coefficient in a generic way.
 * They depend on previous definitions of Fick and Soret coefficients and optionally their
 * derivatives.
 * Parameter variables thermalDiffusionParams and numberThDiffParams must be defined.
 * Argument names list is thermDiffArgs, number is numberThermDiffArgs, ranges are thermDiffRanges.
 * Length of parameter array must be the same as the combined Fick and Soret parameter arrays.
 * The thermal diffusion coefficient calculated herein is the FickCoefficient * SoretCoefficient.
 * Simply include this file in the scope of your material namespace in your material header file.
 * Specializations may be handled within specific material header files in lieu of including this
 * one.
 */

//=================== Classes =======================================================

class ThermalDiffusionCoefficientProp : public AMP::Materials::Property
{
public:
    ThermalDiffusionCoefficientProp()
        : AMP::Materials::Property( std::string( name_base ) +
                                        "_ThermalDiffusionCoefficient", // Name string
                                    source,                             // Reference source
                                    thermalDiffusionParams,             // Property parameters
                                    thermDiffArgs, // Names of arguments to the eval function
                                    thermDiffRanges ),
          d_ExtraParams( 2 )
    {
    }

    virtual void set_parameters( std::vector<double> params )
    {
        std::vector<double> fickParams  = d_FickProp.get_parameters();
        std::vector<double> soretParams = d_SoretProp.get_parameters();
        AMP_ASSERT( params.size() == fickParams.size() + soretParams.size() );
        for ( size_t i = 0; i < fickParams.size(); i++ )
            fickParams[i] = params[i];
        for ( size_t i = 0; i < soretParams.size(); i++ )
            soretParams[i] = params[fickParams.size() + i];
        d_FickProp.set_parameters( std::move( fickParams ) );
        d_SoretProp.set_parameters( std::move( soretParams ) );
    }

    double eval( std::vector<double> &args ) override;

protected:
    virtual void getExtraParameters( std::vector<double> &args )
    {
        d_ExtraParams[0] = d_FickProp.eval( args );
        d_ExtraParams[1] = d_SoretProp.eval( args );
    }

private:
    FickCoefficientProp d_FickProp;
    SoretCoefficientProp d_SoretProp;
    std::vector<double> d_ExtraParams;
};

// define this macro in your host class if you want derivatives in your class and have
// derivatives of Fick and Soret coefficients.
#ifdef THERMAL_DIFFUSION_DERIVATIVE

class DxThermalDiffusionCoefficientProp : public AMP::Materials::Property
{
public:
    DxThermalDiffusionCoefficientProp()
        : AMP::Materials::Property( std::string( name_base ) +
                                        "_DxThermalDiffusionCoefficient", // Name string
                                    source,                               // Reference source
                                    thermalDiffusionParams,               // Property parameters
                                    thermDiffArgs, // Names of arguments to the eval function
                                    thermDiffRanges ),
          d_ExtraParams( 4 )
    {
    } // Number of arguments

    virtual void set_parameters( const double *params, const unsigned int nparams )
    {
        std::vector<double> fickParams  = d_FickProp.get_parameters();
        std::vector<double> soretParams = d_SoretProp.get_parameters();
        AMP_ASSERT( nparams == fickParams.size() + soretParams.size() );
        for ( size_t i = 0; i < fickParams.size(); i++ )
            fickParams[i] = params[i];
        for ( size_t i = 0; i < soretParams.size(); i++ )
            soretParams[i] = params[fickParams.size() + i];
        d_FickProp.set_parameters( std::move( fickParams ) );
        d_SoretProp.set_parameters( std::move( soretParams ) );
    }

    double eval( std::vector<double> &args ) override;

protected:
    virtual void getExtraParameters( std::vector<double> &args )
    {
        d_ExtraParams[0] = d_FickProp.eval( args );
        d_ExtraParams[1] = d_SoretProp.eval( args );
        d_ExtraParams[2] = d_DxFickProp.eval( args );
        d_ExtraParams[3] = d_DxSoretProp.eval( args );
    }

private:
    FickCoefficientProp d_FickProp;
    SoretCoefficientProp d_SoretProp;
    DxFickCoefficientProp d_DxFickProp;
    DxSoretCoefficientProp d_DxSoretProp;
    std::vector<double> d_ExtraParams;
};

class DTThermalDiffusionCoefficientProp : public AMP::Materials::Property
{
public:
    DTThermalDiffusionCoefficientProp()
        : AMP::Materials::Property( std::string( name_base ) +
                                        "_DTThermalDiffusionCoefficient", // Name string
                                    source,                               // Reference source
                                    thermalDiffusionParams,               // Property parameters
                                    thermDiffArgs, // Names of arguments to the eval function
                                    thermDiffRanges ),
          d_ExtraParams( 4 )
    {
    } // Number of arguments

    virtual void set_parameters( std::vector<double> params )
    {
        std::vector<double> fickParams  = d_FickProp.get_parameters();
        std::vector<double> soretParams = d_SoretProp.get_parameters();
        AMP_ASSERT( params.size() == fickParams.size() + soretParams.size() );
        for ( size_t i = 0; i < fickParams.size(); i++ )
            fickParams[i] = params[i];
        for ( size_t i = 0; i < soretParams.size(); i++ )
            soretParams[i] = params[fickParams.size() + i];
        d_FickProp.set_parameters( std::move( fickParams ) );
        d_SoretProp.set_parameters( std::move( soretParams ) );
    }

    double eval( std::vector<double> &args ) override;

protected:
    virtual void getExtraParameters( std::vector<double> &args )
    {
        d_ExtraParams[0] = d_FickProp.eval( args );
        d_ExtraParams[1] = d_SoretProp.eval( args );
        d_ExtraParams[2] = d_DTFickProp.eval( args );
        d_ExtraParams[3] = d_DTSoretProp.eval( args );
    }

private:
    FickCoefficientProp d_FickProp;
    SoretCoefficientProp d_SoretProp;
    DTFickCoefficientProp d_DTFickProp;
    DTSoretCoefficientProp d_DTSoretProp;
    std::vector<double> d_ExtraParams;
};

#endif

//=================== Functions =====================================================

inline double ThermalDiffusionCoefficientProp::eval( std::vector<double> &args )
{
    this->getExtraParameters( args );
    double fick  = d_ExtraParams[0];
    double soret = d_ExtraParams[1];
    return fick * soret;
}

#ifdef THERMAL_DIFFUSION_DERIVATIVE

inline double DxThermalDiffusionCoefficientProp::eval( std::vector<double> &args )
{
    this->getExtraParameters( args );
    double fick    = d_ExtraParams[0];
    double soret   = d_ExtraParams[1];
    double fickDx  = d_ExtraParams[2];
    double soretDx = d_ExtraParams[3];
    return fickDx * soret + fick * soretDx;
}

inline double DTThermalDiffusionCoefficientProp::eval( std::vector<double> &args )
{
    this->getExtraParameters( args );
    double fick    = d_ExtraParams[0];
    double soret   = d_ExtraParams[1];
    double fickDT  = d_ExtraParams[2];
    double soretDT = d_ExtraParams[3];
    return fickDT * soret + fick * soretDT;
}

#endif
