#ifndef included_AMP_ScalarProperty
#define included_AMP_ScalarProperty

#include "AMP/materials/Property.h"
#include "AMP/materials/TensorProperty.h"
#include "AMP/materials/VectorProperty.h"


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
    double eval( const std::vector<double> & ) const override { return d_value; }

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
    double eval( const std::vector<double> &args ) const override
    {
        if ( d_p.size() == 1 )
            return d_p[0];
        double y  = 0;
        double x  = args[0];
        double x2 = 1.0;
        for ( size_t i = 0; i < d_p.size(); i++ ) {
            y += d_p[i] * x2;
            x2 *= x;
        }
        return y;
    }

private:
    std::vector<double> d_p;
};


//! Scalar Vector Property
class ScalarVectorProperty : public VectorProperty
{
public:
    explicit ScalarVectorProperty( const std::string &name,
                                   double value,
                                   const std::string &source = "" )
        : VectorProperty( name, source, {}, {}, 1 ), d_value( value )
    {
    }

    std::vector<double> evalVector( const std::vector<double> & ) override { return { d_value }; }

private:
    double d_value;
};


//! Scalar Tensor Property
class ScalarTensorProperty : public TensorProperty
{
public:
    explicit ScalarTensorProperty( const std::string &name,
                                   const std::string &source,
                                   const std::vector<size_t> &dims,
                                   const std::vector<double> &params )
        : TensorProperty( name, source, {}, {}, dims ), d_params( params )
    {
        AMP_INSIST( d_params.size() == dims[0] * dims[1],
                    "dimensions and number of parameters don't match" );
    }

    std::vector<std::vector<double>> evalTensor( const std::vector<double> & ) override
    {
        std::vector<std::vector<double>> result( d_dimensions[0],
                                                 std::vector<double>( d_dimensions[1] ) );
        for ( size_t i = 0; i < d_dimensions[0]; i++ )
            for ( size_t j = 0; j < d_dimensions[1]; j++ )
                result[i][j] = d_params[i * d_dimensions[1] + j];
        return result;
    }

private:
    std::vector<double> d_params;
};

} // namespace AMP::Materials

#endif
