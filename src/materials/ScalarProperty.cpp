#include "AMP/materials/ScalarProperty.h"
#include "AMP/utils/MathExpr.h"
#include "AMP/utils/Utilities.h"


namespace AMP::Materials {


/*******************************************************************
 *  StringProperty                                                  *
 *******************************************************************/
StringProperty::StringProperty( std::string name, std::string value, std::string source )
    : Property( std::move( name ), { 1 }, {}, std::move( source ) ), d_value( std::move( value ) )
{
}
void StringProperty::eval( AMP::Array<double> &, const AMP::Array<double> & ) const
{
    AMP_ERROR( "numerically evaluating StringProperty is not supported" );
}


/*******************************************************************
 *  ScalarProperty                                                  *
 *******************************************************************/
ScalarProperty::ScalarProperty( std::string name,
                                double value,
                                const AMP::Units &unit,
                                std::string source )
    : Property( std::move( name ), { 1 }, unit, std::move( source ) )
{
    d_value.resize( 1 );
    d_value( 0 ) = value;
}
ScalarProperty::ScalarProperty( std::string name,
                                AMP::Array<double> value,
                                const AMP::Units &unit,
                                std::string source )
    : Property( std::move( name ), value.size(), unit, std::move( source ) ),
      d_value( std::move( value ) )
{
}
void ScalarProperty::eval( AMP::Array<double> &result, const AMP::Array<double> & ) const
{
    size_t N1 = d_value.length();
    size_t N2 = result.size( d_value.ndim() );
    AMP_ASSERT( N1 * N2 == result.length() );
    for ( size_t i = 0; i < N2; i++ )
        memcpy( &result( 0, i ), d_value.data(), N1 * sizeof( double ) );
}


/*******************************************************************
 *  PolynomialProperty                                              *
 *******************************************************************/
PolynomialProperty::PolynomialProperty( std::string name,
                                        std::string source,
                                        const AMP::Units &unit,
                                        std::vector<double> params,
                                        std::vector<std::string> args,
                                        std::vector<std::array<double, 2>> ranges,
                                        std::vector<AMP::Units> argUnits )
    : Property( std::move( name ),
                { 1 },
                unit,
                std::move( source ),
                std::move( args ),
                std::move( ranges ),
                std::move( argUnits ) ),
      d_p( std::move( params ) )
{
    if ( d_p.size() > 1 )
        AMP_ASSERT( d_arguments.size() == 1 );
    else
        AMP_ASSERT( d_arguments.empty() );
}
void PolynomialProperty::eval( AMP::Array<double> &result, const AMP::Array<double> &args ) const
{
    if ( d_p.size() == 1 )
        return result.fill( d_p[0] );
    for ( size_t i = 0; i < result.length(); i++ ) {
        double x  = args( i );
        double y  = 0;
        double x2 = 1.0;
        for ( size_t j = 0; j < d_p.size(); j++ ) {
            y += d_p[j] * x2;
            x2 *= x;
        }
        result( i ) = y;
    }
}


/*******************************************************************
 *  ScalarPropInterpolatedPropertyerty                              *
 *******************************************************************/
InterpolatedProperty::InterpolatedProperty( std::string name,
                                            const AMP::Units &unit,
                                            const std::string &var_name,
                                            std::vector<double> x,
                                            std::vector<double> y,
                                            const std::array<double, 2> range,
                                            const AMP::Units &argUnit,
                                            double default_value,
                                            std::string source )
    : Property( std::move( name ),
                { 1 },
                unit,
                std::move( source ),
                { var_name },
                { range },
                { argUnit } ),
      d_x( std::move( x ) ),
      d_y( std::move( y ) )
{
    set_defaults( { default_value } );
    AMP_ASSERT( d_x.size() == d_y.size() );
    AMP::Utilities::quicksort( d_x, d_y );
    for ( size_t i = 1; i < d_x.size(); i++ )
        AMP_INSIST( d_x[i] != d_x[i - 1], "Values of interpolant must be unique" );
}
void InterpolatedProperty::eval( AMP::Array<double> &result, const AMP::Array<double> &args ) const
{
    for ( size_t i = 0; i < result.length(); i++ )
        result( i ) = AMP::Utilities::linear( d_x, d_y, args( i ) );
}


/*******************************************************************
 *  EquationProperty                                                *
 *******************************************************************/
std::vector<std::array<double, 2>> getDefaultRanges( std::vector<std::array<double, 2>> ranges,
                                                     const std::vector<std::string> &vars )
{
    if ( ranges.empty() )
        ranges.resize( vars.size(), { { -1e100, 1e100 } } );
    return ranges;
}
EquationProperty::EquationProperty( std::string name,
                                    std::shared_ptr<const MathExpr> eq,
                                    const AMP::Units &unit,
                                    std::vector<std::array<double, 2>> ranges,
                                    std::vector<AMP::Units> argUnits,
                                    std::string source )
    : Property( std::move( name ),
                { 1 },
                unit,
                std::move( source ),
                eq->getVars(),
                getDefaultRanges( std::move( ranges ), eq->getVars() ),
                std::move( argUnits ) ),
      d_eq( eq )
{
    AMP_ASSERT( d_eq );
}
EquationProperty::EquationProperty( std::string name,
                                    const std::string &expression,
                                    const AMP::Units &unit,
                                    std::vector<std::string> args,
                                    std::vector<std::array<double, 2>> ranges,
                                    std::vector<AMP::Units> argUnits,
                                    std::string source )
    : Property( std::move( name ),
                { 1 },
                unit,
                std::move( source ),
                std::move( args ),
                getDefaultRanges( std::move( ranges ), args ),
                std::move( argUnits ) ),
      d_eq( std::make_shared<MathExpr>( expression, args ) )
{
    AMP_ASSERT( d_eq );
}
void EquationProperty::eval( AMP::Array<double> &result, const AMP::Array<double> &args ) const
{
    AMP_ASSERT( d_eq );
    for ( size_t i = 0; i < result.length(); i++ ) {
        if ( args.empty() )
            result( i ) = ( *d_eq )();
        else
            result( i ) = ( *d_eq )( &args( 0, i ) );
    }
}


} // namespace AMP::Materials
