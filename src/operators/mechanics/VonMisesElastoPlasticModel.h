
#ifndef included_AMP_VonMisesElastoPlasticModel
#define included_AMP_VonMisesElastoPlasticModel

#include "MechanicsMaterialModel.h"

#include "utils/shared_ptr.h"

#include <vector>

namespace AMP {
namespace Operator {


class VonMisesElastoPlasticModel : public MechanicsMaterialModel
{
public:
    explicit VonMisesElastoPlasticModel(
        const AMP::shared_ptr<MechanicsMaterialModelParameters> & );

    virtual ~VonMisesElastoPlasticModel() {}

    void getConstitutiveMatrix(double *&) override;

    void getConstitutiveMatrixUpdatedLagrangian(double[6][6],
                                                double[3][3]) override;

    void getStressForUpdatedLagrangian(double currentStress[6]) override {
      for (int i = 0; i < 6; i++) {
        currentStress[i] = d_tmp1Stress[(6 * d_gaussPtCnt) + i];
      }
    }

    void getInternalStress(const std::vector<std::vector<double>> &,
                           double *&) override;

    void getInternalStress_UL(const std::vector<std::vector<double>> &,
                              double *&, double[3][3], double[3][3],
                              double) override;

    void getEffectiveStress(double *&) override;

    void getEquivalentStrain(double *&) override;

    void preLinearAssembly() override { d_gaussPtCnt = 0; }

    void postLinearGaussPointOperation() override { d_gaussPtCnt++; }

    void preNonlinearInit(bool, bool) override;

    void nonlinearInitGaussPointOperation(double) override;

    void preNonlinearAssembly() override {
      Plastic_Gauss_Point = 0;
      d_gaussPtCnt = 0;
    }

    void postNonlinearAssembly() override;

    void postNonlinearAssemblyGaussPointOperation() override { d_gaussPtCnt++; }

    void preNonlinearReset() override { d_gaussPtCnt = 0; }

    void postNonlinearResetGaussPointOperation() override { d_gaussPtCnt++; }

    void nonlinearResetGaussPointOperation(
        const std::vector<std::vector<double>> &) override;

    void nonlinearResetGaussPointOperation_UL(
        const std::vector<std::vector<double>> &, double[3][3],
        double[3][3]) override;

    void globalReset() override;

    void postNonlinearReset() override;

    void preNonlinearJacobian() override { d_gaussPtCnt = 0; }

    void postNonlinearJacobianGaussPointOperation() override { d_gaussPtCnt++; }

    void nonlinearJacobianGaussPointOperation(
        const std::vector<std::vector<double>> &) override;

    void nonlinearJacobianGaussPointOperation_UL(
        const std::vector<std::vector<double>> &, double[3][3],
        double[3][3]) override;

    unsigned int getLocalPlasticGaussPointCount() { return Plastic_Gauss_Point; }

    unsigned int getLocalGaussPointCount() { return Total_Gauss_Point; }

    double getFractionPlastic()
    {
        double Plastic_Fraction = ( (double) Plastic_Gauss_Point ) / ( (double) Total_Gauss_Point );
        Plastic_Fraction        = Plastic_Fraction * 100.0;
        return Plastic_Fraction;
    }

protected:
    void radialReturn( const double *stra_np1,
                       double *stre_np1,
                       double *ystre_np1,
                       double *eph_bar_plas_np1,
                       const std::vector<std::vector<double>> &strain,
                       double R_n[3][3],
                       double R_np1[3][3] );

    void constructConstitutiveMatrix();

    double default_TEMPERATURE;

    double default_BURNUP;

    double default_OXYGEN_CONCENTRATION;

    std::vector<double> d_E;

    std::vector<double> d_Nu;

    std::vector<double> d_detULF;

    double default_E;

    double default_Nu;

    double d_H;

    double d_Sig0;

    int mat_name;

    unsigned int
        Total_Gauss_Point; /**< Total how many gauss points are there in this simulation. */

    unsigned int Plastic_Gauss_Point; /**< How many gauss points have reached plasticity at the
                                         current stage. */

    unsigned int d_gaussPtCnt;

    double d_constitutiveMatrix[6][6];

    std::vector<double> d_EquilibriumStress;

    std::vector<double> d_EquilibriumStrain;

    std::vector<double> d_EquilibriumYieldStress;

    std::vector<double> d_EquilibriumEffectivePlasticStrain;

    std::vector<double> d_Lambda;

    std::vector<int> d_ElPl;

    std::vector<double> d_tmp1Stress;

    std::vector<double> d_tmp1Strain;

    std::vector<double> d_tmp1YieldStress;

    std::vector<double> d_tmp1EffectivePlasticStrain;

    std::vector<double> d_tmp2Stress;

    std::vector<double> d_tmp2YieldStress;

    std::vector<double> d_tmp2EffectivePlasticStrain;

    bool d_resetReusesRadialReturn;

    bool d_jacobianReusesRadialReturn;

private:
};
}
}

#endif
