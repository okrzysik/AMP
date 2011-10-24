
#ifndef included_AMP_ThermalStrainMaterialModel
#define included_AMP_ThermalStrainMaterialModel

#include "MechanicsMaterialModel.h"

#include "boost/shared_ptr.hpp"

#include <vector>

namespace AMP {
namespace Operator {

  class ThermalStrainMaterialModel : public MechanicsMaterialModel 
  {
    public :

      ThermalStrainMaterialModel(const boost::shared_ptr<MechanicsMaterialModelParameters>& );

      ~ThermalStrainMaterialModel() { }

      void getConstitutiveMatrix(double*& );

      void getConstitutiveMatrixUpdatedLagrangian(double[6][6], double[3][3]);

      void getStressForUpdatedLagrangian(double currentStress[6]) {
        for(int i = 0; i < 6; i++) {
          currentStress[i] = d_tmp1Stress[(6*d_gaussPtCnt) + i];
        }
      }

      void getExternalStress(double*& );

      void getInternalStress(const std::vector<std::vector<double> >& , double*& );

      void getInternalStress_UL(const std::vector<std::vector<double> >& , double*&, double[3][3], double[3][3], double);

      void preLinearAssembly() {
        d_gaussPtCnt = 0;
      }

      void postLinearGaussPointOperation() {
        d_gaussPtCnt++;
      }

      void preNonlinearInit(bool, bool);

      void nonlinearInitGaussPointOperation(double); 

      void preNonlinearAssembly() {
        d_gaussPtCnt = 0;
      }

      void postNonlinearAssemblyGaussPointOperation() {
        d_gaussPtCnt++;
      }

      void preNonlinearReset() {
        d_gaussPtCnt = 0;
      }

      void postNonlinearResetGaussPointOperation() {
        d_gaussPtCnt++;
      }

      void nonlinearResetGaussPointOperation(const std::vector<std::vector<double> >&); 

      void nonlinearResetGaussPointOperation_UL(const std::vector<std::vector<double> >&, double[3][3], double[3][3] );

      void globalReset();

      void postNonlinearReset();

      void preNonlinearJacobian() {
        d_gaussPtCnt = 0;
      }

      void postNonlinearJacobianGaussPointOperation() {
        d_gaussPtCnt++;
      }

      void nonlinearJacobianGaussPointOperation(const std::vector<std::vector<double> >&);

      void nonlinearJacobianGaussPointOperation_UL(const std::vector<std::vector<double> >&, double[3][3], double[3][3] );

    protected :

      void Thermal_Strain_Gauss_Point(const double* stra_np1, const double Temp, double* stre_np1, 
          const std::vector<std::vector<double> >& strain, double R_n[3][3], double R_np1[3][3]);

      void computeEvalv(const std::vector<std::vector<double> >& strain);

      void constructConstitutiveMatrix(const double, const double);

      void constructConstitutiveMatrixUpdatedLagrangian(const double, const double);

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

    private :

  };

}
}

#endif

