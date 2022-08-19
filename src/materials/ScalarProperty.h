#ifndef included_AMP_ScalarProperty
#define included_AMP_ScalarProperty

#include "AMP/materials/Property.h"
#include "AMP/utils/MathExpr.h"

#include <cstring>


namespace AMP::Materials {


//! Scalar property class (fixed value)
class ScalarProperty final : public Property
{
public:
    ScalarProperty( std::string name,
                    double value,
                    const AMP::Units &unit = AMP::Units(),
                    std::string source     = "" );
    ScalarProperty( std::string name,
                    AMP::Array<double> value,
                    const AMP::Units &unit = AMP::Units(),
                    std::string source     = "" );
    void eval( AMP::Array<double> &result, const AMP::Array<double> & ) const override;
    inline const AMP::Array<double> &getValue() const { return d_value; }

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
                        std::vector<AMP::Units> argUnits          = {} );
    void eval( AMP::Array<double> &result, const AMP::Array<double> &args ) const override;

private:
    std::vector<double> d_p;
};


//! Interpolated property class
class InterpolatedProperty final : public Property
{
public:
    InterpolatedProperty( std::string name,
                          const AMP::Units &unit,
                          const std::string &var_name,
                          std::vector<double> x,
                          std::vector<double> y,
                          const std::array<double, 2> range,
                          const AMP::Units &argUnit,
                          double default_value,
                          std::string source = "" );
    void eval( AMP::Array<double> &result, const AMP::Array<double> & ) const override;

private:
    std::vector<double> d_x;
    std::vector<double> d_y;
};


//! Equation based property class
class EquationProperty : public Property
{
public:
    EquationProperty() {}
    EquationProperty( std::string name,
                      std::shared_ptr<const MathExpr> eq,
                      const AMP::Units &unit                    = {},
                      std::vector<std::array<double, 2>> ranges = {},
                      std::vector<AMP::Units> argUnits          = {},
                      std::string source                        = "" );
    EquationProperty( std::string name,
                      const std::string &expression,
                      const AMP::Units &unit                    = {},
                      std::vector<std::string> args             = {},
                      std::vector<std::array<double, 2>> ranges = {},
                      std::vector<AMP::Units> argUnits          = {},
                      std::string source                        = "" );
    void eval( AMP::Array<double> &result, const AMP::Array<double> &args ) const override;

private:
    std::shared_ptr<const MathExpr> d_eq;
};


} // namespace AMP::Materials

#endif
