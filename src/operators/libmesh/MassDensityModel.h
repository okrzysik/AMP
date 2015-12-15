#ifndef included_AMP_MassDensityModel
#define included_AMP_MassDensityModel

#include <string>
#include <valarray>
#include "operators/ElementPhysicsModel.h"
#include "utils/shared_ptr.h"
#include "materials/Material.h"
#include "materials/Property.h"
#include "utils/ManufacturedSolution.h"
#include "libmesh/point.h"

namespace AMP {
namespace Operator {

typedef ElementPhysicsModelParameters MassDensityModelParameters;

class MassDensityModel: public ElementPhysicsModel
{
public:
    // This was moved from the AMP namespace due to a conflict with
    // MechanicsConstants.h.  --WKC
    enum MassEquation
    {
        Mechanics, Thermal, Chemical, Manufactured, UnknownMassEquation
    };

    enum ManufacturedEquation
    {
        ThermalSrc, FickSrc, SoretSrc, FickSoretSrc, UnknownManufacturedEquation
    };

    explicit MassDensityModel(const AMP::shared_ptr<MassDensityModelParameters>& params);

    virtual ~MassDensityModel()
    {
    }

    virtual void getDensityMechanics(std::vector<double> & result,
            const std::vector<double>& T, const std::vector<double>& U,
            const std::vector<double>& B);

    virtual void getDensityThermal(std::vector<double> & result,
            const std::vector<double>& T, const std::vector<double>& U,
            const std::vector<double>& B);

    virtual void getDensityChemical(std::vector<double> & result,
            const std::vector<double>& T, const std::vector<double>& U,
            const std::vector<double>& B);

    virtual void getDensityManufactured(std::vector<double> & result,
            const std::vector<double>& T, const std::vector<double>& U,
            const std::vector<double>& B, const std::vector<libMesh::Point>& xyz);

    //For LinearOperator's Reset/Init

    virtual void preLinearAssembly()
    {
    }

    virtual void postLinearAssembly()
    {
    }

    virtual void preLinearElementOperation()
    {
    }

    virtual void postLinearElementOperation()
    {
    }

    virtual void preLinearGaussPointOperation()
    {
    }

    virtual void postLinearGaussPointOperation()
    {
    }

    AMP::shared_ptr<ManufacturedSolution> getManufacturedSolution()
            {return d_ManufacturedSolution;}

    MassEquation getEquation(){return d_equation;}

protected:

    AMP::Materials::Material::shared_ptr d_material;

    /**
     * \brief Use a bilogarithmic scaling of material arguments
     *
     * See documentation for DiffusionTransportModel for background.
     * The Chemical mass matrix is 1 of not scaling.
     */
    bool d_UseBilogScaling;

    /**
     * \brief the material argument to which the bilogarithmic transformation applies
     * These must be one of the values returned by Material::get_arg_names().
     */
    std::string d_BilogVariable;

    MassEquation d_equation;

private:

    std::vector<double> d_BilogRange;

    size_t d_BilogIndex;

    AMP::shared_ptr<ManufacturedSolution> d_ManufacturedSolution;

    ManufacturedEquation d_ManufacturedEquation;

    std::string d_ManufacturedVariable;

    bool d_ManufacturedUseTemp;
    bool d_ManufacturedUseConc;

    std::string d_PropertyName;
    std::vector<double> d_Parameters;

};

}
}

#endif

