#ifndef included_AMP_ThermalDiffusionCoefficientProp
#define included_AMP_ThermalDiffusionCoefficientProp

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


class ThermalDiffusionCoefficientProp : public AMP::Materials::Property
{
public:
    ThermalDiffusionCoefficientProp( std::shared_ptr<Property> FickProp,
                                     std::shared_ptr<Property> SoretProp )
        : AMP::Materials::Property( std::string( name_base ) +
                                        "_ThermalDiffusionCoefficient", // Name string
                                    Units(),                            // Units
                                    source,                             // Reference source
                                    thermalDiffusionParams,             // Property parameters
                                    thermDiffArgs, // Names of arguments to the eval function
                                    thermDiffRanges ),
          d_FickProp( FickProp ),
          d_SoretProp( SoretProp ),
          d_ExtraParams( 2 )
    {
    }

    virtual void set_parameters( std::vector<double> params )
    {
        std::vector<double> fickParams  = d_FickProp->get_parameters();
        std::vector<double> soretParams = d_SoretProp->get_parameters();
        AMP_ASSERT( params.size() == fickParams.size() + soretParams.size() );
        for ( size_t i = 0; i < fickParams.size(); i++ )
            fickParams[i] = params[i];
        for ( size_t i = 0; i < soretParams.size(); i++ )
            soretParams[i] = params[fickParams.size() + i];
        d_FickProp->set_parameters( std::move( fickParams ) );
        d_SoretProp->set_parameters( std::move( soretParams ) );
    }

    double eval( const std::vector<double> &args ) override
    {
        this->getExtraParameters( args );
        double fick  = d_ExtraParams[0];
        double soret = d_ExtraParams[1];
        return fick * soret;
    }

protected:
    virtual void getExtraParameters( const std::vector<double> &args )
    {
        d_ExtraParams[0] = d_FickProp->evalDirect( args );
        d_ExtraParams[1] = d_SoretProp->evalDirect( args );
    }

private:
    std::shared_ptr<Property> d_FickProp;
    std::shared_ptr<Property> d_SoretProp;
    std::vector<double> d_ExtraParams;
};


class DxThermalDiffusionCoefficientProp : public AMP::Materials::Property
{
public:
    DxThermalDiffusionCoefficientProp( std::shared_ptr<Property> FickProp,
                                       std::shared_ptr<Property> SoretProp,
                                       std::shared_ptr<Property> DxFickProp,
                                       std::shared_ptr<Property> DxSoretProp )
        : AMP::Materials::Property( std::string( name_base ) +
                                        "_DxThermalDiffusionCoefficient", // Name string
                                    Units(),                              // Units
                                    source,                               // Reference source
                                    thermalDiffusionParams,               // Property parameters
                                    thermDiffArgs, // Names of arguments to the eval function
                                    thermDiffRanges ),
          d_FickProp( FickProp ),
          d_SoretProp( SoretProp ),
          d_DxFickProp( DxFickProp ),
          d_DxSoretProp( DxSoretProp ),
          d_ExtraParams( 4 )
    {
    } // Number of arguments

    virtual void set_parameters( const double *params, const unsigned int nparams )
    {
        std::vector<double> fickParams  = d_FickProp->get_parameters();
        std::vector<double> soretParams = d_SoretProp->get_parameters();
        AMP_ASSERT( nparams == fickParams.size() + soretParams.size() );
        for ( size_t i = 0; i < fickParams.size(); i++ )
            fickParams[i] = params[i];
        for ( size_t i = 0; i < soretParams.size(); i++ )
            soretParams[i] = params[fickParams.size() + i];
        d_FickProp->set_parameters( std::move( fickParams ) );
        d_SoretProp->set_parameters( std::move( soretParams ) );
    }

    double eval( const std::vector<double> &args ) override
    {
        this->getExtraParameters( args );
        double fick    = d_ExtraParams[0];
        double soret   = d_ExtraParams[1];
        double fickDx  = d_ExtraParams[2];
        double soretDx = d_ExtraParams[3];
        return fickDx * soret + fick * soretDx;
    }

protected:
    virtual void getExtraParameters( const std::vector<double> &args )
    {
        d_ExtraParams[0] = d_FickProp->evalDirect( args );
        d_ExtraParams[1] = d_SoretProp->evalDirect( args );
        d_ExtraParams[2] = d_DxFickProp->evalDirect( args );
        d_ExtraParams[3] = d_DxSoretProp->evalDirect( args );
    }

private:
    std::shared_ptr<Property> d_FickProp;
    std::shared_ptr<Property> d_SoretProp;
    std::shared_ptr<Property> d_DxFickProp;
    std::shared_ptr<Property> d_DxSoretProp;
    std::vector<double> d_ExtraParams;
};

class DTThermalDiffusionCoefficientProp : public AMP::Materials::Property
{
public:
    DTThermalDiffusionCoefficientProp( std::shared_ptr<Property> FickProp,
                                       std::shared_ptr<Property> SoretProp,
                                       std::shared_ptr<Property> DTFickProp,
                                       std::shared_ptr<Property> DTSoretProp )
        : AMP::Materials::Property( std::string( name_base ) +
                                        "_DTThermalDiffusionCoefficient", // Name string
                                    Units(),                              // Units
                                    source,                               // Reference source
                                    thermalDiffusionParams,               // Property parameters
                                    thermDiffArgs, // Names of arguments to the eval function
                                    thermDiffRanges ),
          d_FickProp( FickProp ),
          d_SoretProp( SoretProp ),
          d_DTFickProp( DTFickProp ),
          d_DTSoretProp( DTSoretProp ),
          d_ExtraParams( 4 )
    {
    } // Number of arguments

    virtual void set_parameters( std::vector<double> params )
    {
        std::vector<double> fickParams  = d_FickProp->get_parameters();
        std::vector<double> soretParams = d_SoretProp->get_parameters();
        AMP_ASSERT( params.size() == fickParams.size() + soretParams.size() );
        for ( size_t i = 0; i < fickParams.size(); i++ )
            fickParams[i] = params[i];
        for ( size_t i = 0; i < soretParams.size(); i++ )
            soretParams[i] = params[fickParams.size() + i];
        d_FickProp->set_parameters( std::move( fickParams ) );
        d_SoretProp->set_parameters( std::move( soretParams ) );
    }

    double eval( const std::vector<double> &args ) override
    {
        this->getExtraParameters( args );
        double fick    = d_ExtraParams[0];
        double soret   = d_ExtraParams[1];
        double fickDT  = d_ExtraParams[2];
        double soretDT = d_ExtraParams[3];
        return fickDT * soret + fick * soretDT;
    }

protected:
    virtual void getExtraParameters( const std::vector<double> &args )
    {
        d_ExtraParams[0] = d_FickProp->evalDirect( args );
        d_ExtraParams[1] = d_SoretProp->evalDirect( args );
        d_ExtraParams[2] = d_DTFickProp->evalDirect( args );
        d_ExtraParams[3] = d_DTSoretProp->evalDirect( args );
    }

private:
    std::shared_ptr<Property> d_FickProp;
    std::shared_ptr<Property> d_SoretProp;
    std::shared_ptr<Property> d_DTFickProp;
    std::shared_ptr<Property> d_DTSoretProp;
    std::vector<double> d_ExtraParams;
};


#endif
