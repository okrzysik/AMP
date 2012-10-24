
#ifndef included_AMP_ElasticDamageThermalStrainModel
#define included_AMP_ElasticDamageThermalStrainModel

#include "MechanicsMaterialModel.h"

#include "boost/shared_ptr.hpp"

#include <vector>

namespace AMP {
namespace Operator {

  class ElasticDamageThermalStrainModel : public MechanicsMaterialModel 
  {
    public :

      ElasticDamageThermalStrainModel(const boost::shared_ptr<MechanicsMaterialModelParameters>& );

      virtual ~ElasticDamageThermalStrainModel() { }

      void getConstitutiveMatrix(double*& ); 

      void getExternalStress(double*& );

      void getInternalStress(const std::vector<std::vector<double> >& , double*& );

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

      void globalReset();

      void postNonlinearReset();

      void preNonlinearJacobian() {
        d_gaussPtCnt = 0;
      }

      void postNonlinearJacobianGaussPointOperation() {
        d_gaussPtCnt++;
      }

      void nonlinearJacobianGaussPointOperation(const std::vector<std::vector<double> >&);

      std::vector<double> d_EquilibriumDamage;

      std::vector<double> d_tmp1Damage;

      std::vector<double> d_tmp2Damage;

    protected :

      void Thermal_Strain_Gauss_Point(const double* stra_np1, const double Temp, double* stre_np1, 
          const std::vector<std::vector<double> >& strain);

      void computeEvalv(const std::vector<std::vector<double> >& strain);

      void constructConstitutiveMatrix(const double, const double);

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

      double d_DamageThreshold;

      double d_CriticalDamageThreshold;

      unsigned int d_gaussPtCnt;

      double d_constitutiveMatrix[6][6];

      double d_initialConstitutiveMatrix[6][6];

      std::vector<double> d_EquilibriumStress;

      std::vector<double> d_EquilibriumStrain;

      std::vector<double> d_EquilibriumTemperature;

      std::vector<double> d_tmp1Stress;

      std::vector<double> d_tmp1Strain;

      std::vector<double> d_tmp1Temperature;

      std::vector<double> d_tmp2Stress;

      std::vector<double> d_EquilibriumTau;

      std::vector<double> d_tmp1Tau;

      std::vector<double> d_tmp2Tau;

      std::vector<double> d_EquilibriumDamageThreshold;

      std::vector<double> d_tmp1DamageThreshold;

      std::vector<double> d_tmp2DamageThreshold;

      std::vector<double> d_tmp1ThermalStrain_Axial;

      std::vector<double> d_tmp1ThermalStrain_Radial;

      std::vector<double> d_EquilibriumThermalStrain_Axial;

      std::vector<double> d_EquilibriumThermalStrain_Radial;

      std::vector<double> d_InitialDamageVec;

      std::vector<double> d_CriticalDamageVec;

      bool d_resetReusesRadialReturn;

      bool d_jacobianReusesRadialReturn;

      // If Source then TRUE, if Material_Model then FALSE
      bool d_Is_Source;

    private :

  };

}
}

#endif



