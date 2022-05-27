// Define the materials in AMP
#ifndef included_AMP_Materials_CylindricallySymmetric
#define included_AMP_Materials_CylindricallySymmetric

#include "AMP/materials/TensorProperty.h"


namespace AMP::Materials {


/** full cylindrically symmetric tensor diffusion coefficient
 *
 * The parameters are set as follows:
 * params[0] = number of parameters for radial
 * params[1]...params[ params[0] ] = parameters for radial
 * the rest are for the longitudinal
 * AuxiliaryInteger data "derivative" has values 0, 1, 2 for
 * zeroth, r- and z- derivatives, respectively.
 */
class CylindricallySymmetricTensor : public TensorProperty
{
public:
    CylindricallySymmetricTensor( const std::string &name,
                                  std::vector<double> params = { 1, 1, 1 } );

    void set_parameters_and_number( std::vector<double> params ) override;

    std::vector<std::vector<double>> evalTensor( const std::vector<double> &args ) override;

private:
    std::vector<double> d_paramsRadial;
    std::vector<double> d_paramsLongitudinal;
};


} // namespace AMP::Materials

#endif
