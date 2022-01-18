
#ifndef included_AMP_ThermalStrainMaterialModel
#define included_AMP_ThermalStrainMaterialModel

#include "MechanicsMaterialModel.h"

#include <memory>

#include <vector>

namespace AMP::Operator {

class ThermalStrainMaterialModel : public MechanicsMaterialModel
{
public:
    explicit ThermalStrainMaterialModel( std::shared_ptr<MechanicsMaterialModelParameters> );

    virtual ~ThermalStrainMaterialModel() {}

    void getConstitutiveMatrix( double *& ) override;

    void getConstitutiveMatrixUpdatedLagrangian( double[6][6], double[3][3] ) override;

    void getStressForUpdatedLagrangian( double currentStress[6] ) override
    {
        for ( int i = 0; i < 6; i++ ) {
            currentStress[i] = d_tmp1Stress[( 6 * d_gaussPtCnt ) + i];
        }
    }

    void getExternalStress( double *& ) override;

    void getInternalStress( const std::vector<std::vector<double>> &, double *& ) override;

    void getInternalStress_UL( const std::vector<std::vector<double>> &,
                               double *&,
                               double[3][3],
                               double[3][3],
                               double ) override;

    void preLinearAssembly() override { d_gaussPtCnt = 0; }

    void postLinearGaussPointOperation() override { d_gaussPtCnt++; }

    void preNonlinearInit( bool, bool ) override;

    void nonlinearInitGaussPointOperation( double ) override;

    void preNonlinearAssembly() override { d_gaussPtCnt = 0; }

    void postNonlinearAssemblyGaussPointOperation() override { d_gaussPtCnt++; }

    void preNonlinearReset() override { d_gaussPtCnt = 0; }

    void postNonlinearResetGaussPointOperation() override { d_gaussPtCnt++; }

    void nonlinearResetGaussPointOperation( const std::vector<std::vector<double>> & ) override;

    void nonlinearResetGaussPointOperation_UL( const std::vector<std::vector<double>> &,
                                               double[3][3],
                                               double[3][3] ) override;

    void globalReset() override;

    void postNonlinearReset() override;

    void preNonlinearJacobian() override { d_gaussPtCnt = 0; }

    void postNonlinearJacobianGaussPointOperation() override { d_gaussPtCnt++; }

    void nonlinearJacobianGaussPointOperation( const std::vector<std::vector<double>> & ) override;

    void nonlinearJacobianGaussPointOperation_UL( const std::vector<std::vector<double>> &,
                                                  double[3][3],
                                                  double[3][3] ) override;

protected:
    void Thermal_Strain_Gauss_Point( const double *stra_np1,
                                     const double Temp,
                                     double *stre_np1,
                                     const std::vector<std::vector<double>> &strain,
                                     double R_n[3][3],
                                     double R_np1[3][3] );

    void computeEvalv( const std::vector<std::vector<double>> &strain );

    void constructConstitutiveMatrix( const double, const double );

    void constructConstitutiveMatrixUpdatedLagrangian( const double, const double );

    double default_TEMPERATURE;

    double default_BURNUP;

    double default_OXYGEN_CONCENTRATION;

    // Thermal expansion coefficient.
    std::vector<double> d_alpha;

    std::vector<double> d_E;

    std::vector<double> d_Nu;

    double default_E;

    double default_Nu;

    double default_alpha;

    unsigned int d_gaussPtCnt;

    double d_constitutiveMatrix[6][6];

    double d_constitutiveMatrix_UL[6][6];

    std::vector<double> d_EquilibriumStress;

    std::vector<double> d_EquilibriumStrain;

    std::vector<double> d_EquilibriumTemperature;

    std::vector<double> d_tmp1Stress;

    std::vector<double> d_tmp1Strain;

    std::vector<double> d_tmp1Temperature;

    std::vector<double> d_detULF;

    bool d_resetReusesRadialReturn;

    bool d_jacobianReusesRadialReturn;

    // If Source then TRUE, if Material_Model then FALSE
    bool d_Is_Source;

private:
};
} // namespace AMP::Operator

#endif
