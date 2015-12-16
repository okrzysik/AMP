
#ifndef included_AMP_VonMises_IsotropicKinematicHardening
#define included_AMP_VonMises_IsotropicKinematicHardening

#include "MechanicsMaterialModel.h"

#include "utils/shared_ptr.h"

#include <vector>

namespace AMP {
namespace Operator {

class VonMises_IsotropicKinematicHardening : public MechanicsMaterialModel {
public:
    explicit VonMises_IsotropicKinematicHardening(
        const AMP::shared_ptr<MechanicsMaterialModelParameters> & );

    virtual ~VonMises_IsotropicKinematicHardening() {}

    void getConstitutiveMatrix( double *& );

    void getInternalStress( const std::vector<std::vector<double>> &, double *& );

    void preLinearAssembly() { d_gaussPtCnt = 0; }

    void postLinearGaussPointOperation() { d_gaussPtCnt++; }

    void preNonlinearInit( bool, bool );

    void nonlinearInitGaussPointOperation( double );

    void preNonlinearAssembly()
    {
        Plastic_Gauss_Point = 0;

        d_gaussPtCnt = 0;
    }

    void postNonlinearAssemblyGaussPointOperation() { d_gaussPtCnt++; }

    void postNonlinearAssembly();

    void preNonlinearReset() { d_gaussPtCnt = 0; }

    void postNonlinearResetGaussPointOperation() { d_gaussPtCnt++; }

    void nonlinearResetGaussPointOperation( const std::vector<std::vector<double>> & );

    void globalReset();

    void postNonlinearReset();

    void preNonlinearJacobian() { d_gaussPtCnt = 0; }

    void postNonlinearJacobianGaussPointOperation() { d_gaussPtCnt++; }

    void nonlinearJacobianGaussPointOperation( const std::vector<std::vector<double>> & );

protected:
    void radialReturn( const double *stra_np1,
                       double *stre_np1,
                       double *back_stress_np1,
                       double *ystre_np1,
                       double *eph_bar_plas_np1,
                       const std::vector<std::vector<double>> &strain );

    void computeEvalv( const std::vector<std::vector<double>> &strain );

    double default_TEMPERATURE;

    double default_BURNUP;

    double default_OXYGEN_CONCENTRATION;

    std::vector<double> d_E;

    std::vector<double> d_Nu;

    double default_E;

    double default_Nu;

    double d_delta;

    double d_K_0;

    double d_K_inf;

    double d_H;

    double d_beta;

    double d_Ep;

    double d_Kin;

    double d_Sig_0;

    unsigned int d_gaussPtCnt;

    unsigned int
        Total_Gauss_Point; /**< Total how many gauss points are there in this simulation. */

    unsigned int Plastic_Gauss_Point; /**< How many gauss points have reached plasticity at the
                                         current stage. */

    double d_constitutiveMatrix[6][6];

    std::vector<double> d_EquilibriumStress;

    std::vector<double> d_EquilibriumBackStress;

    std::vector<double> d_EquilibriumStrain;

    std::vector<double> d_EquilibriumYieldStress;

    std::vector<double> d_EquilibriumEffectivePlasticStrain;

    std::vector<double> d_Lambda;

    std::vector<int> d_ElPl;

    std::vector<double> d_tmp1Stress;

    std::vector<double> d_tmp1BackStress;

    std::vector<double> d_tmp1Strain;

    std::vector<double> d_tmp1YieldStress;

    std::vector<double> d_tmp1EffectivePlasticStrain;

    std::vector<double> d_tmp2Stress;

    std::vector<double> d_tmp2BackStress;

    std::vector<double> d_tmp2YieldStress;

    std::vector<double> d_tmp2EffectivePlasticStrain;

    bool d_resetReusesRadialReturn;

    bool d_jacobianReusesRadialReturn;

    bool d_CM_Test;

    bool d_TW_Test;

private:
};
}
}

#endif
