
#ifndef included_AMP_IsotropicElasticModel
#define included_AMP_IsotropicElasticModel

#include "MechanicsMaterialModel.h"

#include <memory>

namespace AMP::Operator {

class IsotropicElasticModel : public MechanicsMaterialModel
{
public:
    explicit IsotropicElasticModel( std::shared_ptr<MechanicsMaterialModelParameters> );

    virtual ~IsotropicElasticModel() {}

    void getConstitutiveMatrix( double *& ) override;

    void getConstitutiveMatrixUpdatedLagrangian( double[6][6], double[3][3] ) override;

    void getStressForUpdatedLagrangian( double currentStress[6] ) override
    {
        for ( int i = 0; i < 6; i++ ) {
            // currentStress[i] = d_EquilibriumStress[(6*d_gaussPtCnt) + i];
            currentStress[i] = d_tmp1Stress[( 6 * d_gaussPtCnt ) + i];
        }
    }

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

    double getYoungsModulus()
    {
        if ( d_E.empty() ) {
            return default_E;
        } else {
            return d_E[d_gaussPtCnt];
        }
    }

    double getPoissonsRatio()
    {
        if ( d_Nu.empty() ) {
            return default_Nu;
        } else {
            return d_Nu[d_gaussPtCnt];
        }
    }

    void nonlinearJacobianGaussPointOperation( const std::vector<std::vector<double>> & ) override;

    void nonlinearJacobianGaussPointOperation_UL( const std::vector<std::vector<double>> &,
                                                  double[3][3],
                                                  double[3][3] ) override;

protected:
    void computeEvalv( const std::vector<std::vector<double>> & );

    void calculateStress( const std::vector<std::vector<double>> &, double *& );

    void calculateStress( const std::vector<std::vector<double>> &,
                          double *&,
                          double[3][3],
                          double[3][3] );

    void constructConstitutiveMatrix( const double, const double );

    void constructConstitutiveMatrixUpdatedLagrangian( const double, const double );

    double default_TEMPERATURE;

    double default_BURNUP;

    double default_OXYGEN_CONCENTRATION;

    std::vector<double> d_E;

    std::vector<double> d_Nu;

    std::vector<double> d_detULF;

    double default_E;

    double default_Nu;

    unsigned int d_gaussPtCnt;

    double d_constitutiveMatrix[6][6];

    double d_constitutiveMatrix_UL[6][6];

    std::vector<double> d_EquilibriumStress;

    std::vector<double> d_EquilibriumStrain;

    std::vector<double> d_tmp1Stress;

    std::vector<double> d_tmp1Strain;

    bool d_resetReusesRadialReturn;

    bool d_jacobianReusesRadialReturn;

    // FILE *fout;

private:
};
} // namespace AMP::Operator

#endif
