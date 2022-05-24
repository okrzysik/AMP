#ifndef included_AMP_ScalarProperty
#define included_AMP_ScalarProperty

#include "AMP/materials/Property.h"


//! Macro to register a scalar property
#define registerScalarProperty( PROPERTY, VALUE, ... )                        \
    d_propertyMap[PROPERTY] = std::make_shared<ScalarProperty>(               \
        AMP::Utilities::demangle( typeid( *this ).name() ) + "::" + PROPERTY, \
        VALUE,                                                                \
        __VA_ARGS__ );

//! Macro to register a polynomial based property
#define registerPolynomialProperty( PROPERTY, ... )                 \
    d_propertyMap[PROPERTY] = std::make_shared<PolynomialProperty>( \
        AMP::Utilities::demangle( typeid( *this ).name() ) + "::" + PROPERTY, __VA_ARGS__ );


namespace AMP::Materials {


//! Scalar property class (fixed value)
class ScalarProperty final : public Property
{
public:
    ScalarProperty( std::string name,
                    double value,
                    const AMP::Units &unit = AMP::Units(),
                    std::string source     = "" )
        : Property( std::move( name ), unit, std::move( source ) ), d_value( value )
    {
    }
    double eval( const std::vector<double> & ) override { return d_value; }

private:
    double d_value;
};


//! Polynomial based property class
class PolynomialProperty final : public Property
{
public:
    PolynomialProperty() {}
    PolynomialProperty( std::string name,
                        std::string source,
                        std::vector<double> params,
                        std::vector<std::string> args             = {},
                        std::vector<std::array<double, 2>> ranges = {},
                        const AMP::Units &unit                    = AMP::Units() )
        : Property( std::move( name ),
                    unit,
                    std::move( source ),
                    std::move( params ),
                    std::move( args ),
                    std::move( ranges ) )
    {
        if ( d_params.size() > 1 )
            AMP_ASSERT( d_arguments.size() == 1 );
        else
            AMP_ASSERT( d_arguments.empty() );
    }
    double eval( const std::vector<double> &args ) override
    {
        if ( d_params.size() == 1 )
            return d_params[0];
        double y  = 0;
        double x  = args[0];
        double x2 = 1.0;
        for ( size_t i = 0; i < d_params.size(); i++ ) {
            y += d_params[i] * x2;
            x2 *= x;
        }
        return y;
    }
};


} // namespace AMP::Materials

#endif
