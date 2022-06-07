#ifndef included_AMP_ScalarProperty
#define included_AMP_ScalarProperty

#include "AMP/materials/Property.h"


namespace AMP::Materials {


//! Scalar property class (fixed value)
class ScalarProperty final : public Property
{
public:
    ScalarProperty( std::string name,
                    double value,
                    const AMP::Units &unit = AMP::Units(),
                    std::string source     = "" )
        : Property( std::move( name ), { 1 }, unit, std::move( source ) ), d_value( value )
    {
    }
    void eval( AMP::Array<double> &result, const AMP::Array<double> & ) const override
    {
        result.fill( d_value );
    }

private:
    double d_value;
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


//! Scalar Vector Property
class ScalarVectorProperty : public Property
{
public:
    explicit ScalarVectorProperty( const std::string &name,
                                   double value,
                                   const std::string &source = "" )
        : Property( name, { 1 }, {}, source, {}, {} ), d_value( value )
    {
    }

    void eval( AMP::Array<double> &result, const AMP::Array<double> & ) const override
    {
        result.fill( d_value );
    }

private:
    double d_value;
};


//! Scalar Tensor Property
class ScalarTensorProperty : public Property
{
public:
    explicit ScalarTensorProperty( const std::string &name,
                                   const std::string &source,
                                   const AMP::ArraySize &dims,
                                   const std::vector<double> &params )
        : Property( name, dims, {}, source, {}, {} ), d_params( params )
    {
        AMP_INSIST( d_params.size() == dims[0] * dims[1],
                    "dimensions and number of parameters don't match" );
    }

    void eval( AMP::Array<double> &result, const AMP::Array<double> & ) const override
    {
        for ( size_t i = 0; i < d_dim[0]; i++ ) {
            for ( size_t j = 0; j < d_dim[1]; j++ ) {
                for ( size_t k = 0; k < result.size( 2 ); k++ ) {
                    result( i, j, k ) = d_params[i * d_dim[1] + j];
                }
            }
        }
    }

private:
    std::vector<double> d_params;
};

} // namespace AMP::Materials

#endif
