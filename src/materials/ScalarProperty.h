#ifndef included_AMP_ScalarProperty
#define included_AMP_ScalarProperty

#include "AMP/materials/Property.h"

#include <cstring>


namespace AMP::Materials {


//! Scalar property class (fixed value)
class ScalarProperty final : public Property
{
public:
    ScalarProperty( std::string name,
                    double value,
                    const AMP::Units &unit = AMP::Units(),
                    std::string source     = "" )
        : Property( std::move( name ), { 1 }, unit, std::move( source ) )
    {
        d_value.resize( 1 );
        d_value( 0 ) = value;
    }
    ScalarProperty( std::string name,
                    AMP::Array<double> value,
                    const AMP::Units &unit = AMP::Units(),
                    std::string source     = "" )
        : Property( std::move( name ), value.size(), unit, std::move( source ) ),
          d_value( std::move( value ) )
    {
    }
    void eval( AMP::Array<double> &result, const AMP::Array<double> & ) const override
    {
        size_t N1 = d_value.length();
        size_t N2 = result.size( d_value.ndim() + 1 );
        for ( size_t i = 0; i < N2; i++ ) {
            memcpy( &result( 0, i ), d_value.data(), N1 * sizeof( double ) );
            result.fill( d_value );
        }
    }

private:
    AMP::Array<double> d_value;
};


//! Polynomial based property class
class PolynomialProperty : public Property
{
public:
    PolynomialProperty() {}
    PolynomialProperty( std::string name,
                        std::string source,
                        const AMP::Units &unit                    = {},
                        std::vector<double> params                = {},
                        std::vector<std::string> args             = {},
                        std::vector<std::array<double, 2>> ranges = {},
                        std::vector<AMP::Units> argUnits          = {} )
        : Property( std::move( name ),
                    { 1 },
                    unit,
                    std::move( source ),
                    std::move( args ),
                    std::move( ranges ),
                    std::move( argUnits ) ),
          d_p( params )
    {
        if ( d_p.size() > 1 )
            AMP_ASSERT( d_arguments.size() == 1 );
        else
            AMP_ASSERT( d_arguments.empty() );
    }
    void eval( AMP::Array<double> &result, const AMP::Array<double> &args ) const override
    {
        if ( d_p.size() == 1 )
            return result.fill( d_p[0] );
        for ( size_t i = 0; i < result.length(); i++ ) {
            double x  = args( i );
            double y  = 0;
            double x2 = 1.0;
            for ( size_t i = 0; i < d_p.size(); i++ ) {
                y += d_p[i] * x2;
                x2 *= x;
            }
            result( i ) = y;
        }
    }

private:
    std::vector<double> d_p;
};


} // namespace AMP::Materials

#endif
