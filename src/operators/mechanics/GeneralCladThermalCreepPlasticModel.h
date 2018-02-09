
#ifndef included_AMP_GeneralCladThermalCreepPlasticModel
#define included_AMP_GeneralCladThermalCreepPlasticModel

#include "MechanicsMaterialModel.h"

#include "AMP/utils/shared_ptr.h"

#include <vector>

namespace AMP {
namespace Operator {

/**
  \brief Non-linear elasto-plastic mechanics with only isotropic options for thermal and creep
  strains.

  This class represents the mechanics non-linear elasto-plastic material model. This has a
  capability of
  capturing the thermal expansion of the clad in a isotropic sense. Other than
  thermal expansion, creep and plastic deformation of the clad can be captured using this model as
  well.
  */
class GeneralCladThermalCreepPlasticModel : public MechanicsMaterialModel
{
public:
    /**
      Constructor. This reads the values of some parameters from the database object contained in
      the
      parameter object (params), or sets those parameters to default values if they are not
      specified
      in the database. Some of these parameters might be overwritten by values from the materials
      library.
      The parameters are as follows:

      (1) USE_MATERIALS_LIBRARY {d_useMaterialsLibrary} (FALSE by default) - This must be set to
      TRUE for the
      current material.

      (2) Material {matname} - Name of the material that is being used to calculate the material
      properties for
      the current material model. This must be a clad material.

      (3) Use_Thermal_Strain {d_UseThermalStrain} (TRUE by default) - This parameter is set to TRUE
      if the strain due to
      thermal
      expansion has to be taken into account. Without thermal strain, clad model does not work.

      (4) Use_Creep_Strain {d_UseCreepStrain} (FALSE by default) - This parameter is set to TRUE to
      take into
      consideration
      the irrecoverable strain due to thermal and stress induced creep inside the clad.

      (5) Linear_Strain_Hardening {d_H} [5.0e8 (Pa) by default] - This one signifies the strain
      hardening coefficient
      during the plastic deformation.

      (6) Plastic_Strain_Exponent {d_n} [1.0 by default] - The equivalent plastic strain exponent
      when calculating
      the modified yield stress during plastic deformation.

      (7) Youngs_Modulus {default_E} [2.08e11 (Pa) by default] - The default value of Young's
      Modulus.

      (8) Poissons_Ratio {default_Nu} [0.23 by default] - The default value of Poisson's Ratio.

      (9) THERMAL_EXPANSION_COEFFICIENT {default_alpha} [2.0e-6 (/ K) by default] - Isotropic
      thermal expansion
      coefficient.

      (11) Elastic_Yield_Stress {d_Sig0} [3.45e8 (Pa) by default] - The limit where elastic
      deformation ends and
      plastic deformation starts.

      (12) Default_Oxygen_Concentration {default_OXYGEN_CONCENTRATION} [0.0 (unit-less) by default]
      - This states the
      default
      oxygen concentration if the "oxygen_concentration" variable is inactive.

      (13) Default_Temperature {default_TEMPERATURE} [310.0 (K) by default] - This gives the default
      temperature if the
      "temperature" variable is inactive.

      (14) Default_Burnup {default_BURNUP} [0.0 (GWd/MTU) by default] - This gives the default
      burnup at each gauss
      point
      if the "burnup" variable is inactive.

      (15) Creep_Delta_Time {d_Delta_Time} [1.0 days (d) by default] - This is used in the
      calculation of the
      creep strain inside the clad.
      */
    explicit GeneralCladThermalCreepPlasticModel(
        const AMP::shared_ptr<MechanicsMaterialModelParameters> & );

    /**
      Destructor.
      */
    virtual ~GeneralCladThermalCreepPlasticModel() {}

    /**
      Calculates the linearized elasto-plastic stress-strain constitutive matrix at each gauss
      point which enters the computation of the jacobian for the element.
      */
    void getConstitutiveMatrix( double *& ) override;

    /**
      Given a certain strain, this function calculates the corresponding stress taking into
      consideration the
      thermal swelling and creep effects. Creep is modeled in an implicit fashion.
      */
    void getInternalStress( const std::vector<std::vector<double>> &, double *& ) override;

    /**
      This function is called just before the linear assembly.
      */
    void preLinearAssembly() override { d_gaussPtCnt = 0; }

    /**
      This is called after the apply() at each gauss point is called.
      */
    void postLinearGaussPointOperation() override { d_gaussPtCnt++; }

    /**
      This function is called before the NonlinearInit is invoked. Here the memory
      assigned to all the gauss point vectors are cleared.
      */
    void preNonlinearInit( bool, bool ) override;

    /**
      This function initializes all the gauss point vectors with the default values.
      */
    void nonlinearInitGaussPointOperation( double ) override;

    /**
      This is called before every non-linear assembly. Here the gauss point count and the total
      number of gauss points which have reached plasticity, are initialized to zero.
      */
    void preNonlinearAssembly() override
    {

        Plastic_Gauss_Point = 0;

        d_gaussPtCnt = 0;
    }

    /**
      How many gauss points have reached plasticity, is calculated in this function (done after each
      assembly).
      */
    void postNonlinearAssembly() override;

    /**
      After the apply at each gauss point , this function increments the gauss point count by one.
      */
    void postNonlinearAssemblyGaussPointOperation() override { d_gaussPtCnt++; }

    /**
      This function initializes the gauss point count for the nonlinear reset.
      */
    void preNonlinearReset() override { d_gaussPtCnt = 0; }

    /**
      Incrementing the gauss point count by unity is conducted in this function.
      */
    void postNonlinearResetGaussPointOperation() override { d_gaussPtCnt++; }

    /**
      At each gauss point the "radialReturn" function is called from this function.
      */
    void nonlinearResetGaussPointOperation( const std::vector<std::vector<double>> & ) override;

    /**
      This function updates the old equilibrium values with the new one under certain conditions.
      */
    void globalReset() override;

    /**
      If the updating is not done in the global reset, it is done in this function.
      */
    void postNonlinearReset() override;

    /**
      Initializes the gauss point count before the nonlinear jacobian is called.
      */
    void preNonlinearJacobian() override { d_gaussPtCnt = 0; }

    /**
      The gauss point count is incremented after computing the nonlinear jacobian at every gauss
      point.
      */
    void postNonlinearJacobianGaussPointOperation() override { d_gaussPtCnt++; }

    /**
      This function updates all the material parameters and calls the radialReturn before
      calculating the jacobian.
      */
    void nonlinearJacobianGaussPointOperation( const std::vector<std::vector<double>> & ) override;

protected:
    /**
      This function calls the evalv functions in the materials library and updates the material
      parameters, such as,
      Young's Modulus, Poisson's Ratio and Thermal Expansion Coefficient. These material parameters
      corresponding
      to the current temperature, burnup and oxygen concentration is extracted from the materials
      library.
      */
    void computeEvalv( const std::vector<std::vector<double>> & );

    /**
      Given the total strain, this function calculates the stress which is used to compute the
      internal force.
      The formula is as follows, \f$f_{int} = {int}_{V} B^T \sigma dv\f$.
      */
    void radialReturn( const double *stra_np1,
                       double *stre_np1,
                       double *ystre_np1,
                       double *eph_bar_plas_np1 );

    /**
      This function is used to calculate the implicit creep strain rate. It returns a creep strain
      rate at a given
      temperature and effective stress (which enters the function as input parameters).
      */
    double help_compute_E1( const double Temp_np1, double effective_stress );

    /**
      This function calculates the sign of a double which it takes as an input.
      */
    int det_sign( double a );

    /**
      This is used to evaluate the implicit creep strain rate. This function returns the first
      derivative
      of E1 with respect to the effective stress.
      */
    double compute_dE1_dsig_e( double effective_stress_trial,
                               double effective_stress,
                               double Temp_np1,
                               double G );

    /**
      This function calculates the creep strain increment (\f$\dot{epsilon}^c\f$) at each time step
      using an
      implicit scheme. The formula is a general one: \f$\dot{\epsilon}^c = A\sigma_e^n
      exp(\frac{-Q}{RT})\f$.
      */
    void computeCreepStrain( const double Temp_np1,
                             const double stress_n[6],
                             const double creep_strain_prev,
                             double delta_creep_strain[6],
                             double net_stra_np1[6] );

    double default_TEMPERATURE; /**< The default temperature which is passed to the materials
                                   library when evalv is
                                   called. */

    double default_OXYGEN_CONCENTRATION; /**< The default oxygen_concentration which is passed to
                                            the materials library
                                            when evalv is called. */

    double default_BURNUP; /**< The default burnup which is passed to the materials library when
                              evalv is called. */

    double default_E; /**< The default value of the Young's Modulus. */

    double default_Nu; /**< The default value of the Poisson'd Ratio. */

    double d_H; /**< Default value for the linear strain hardening coefficient. */

    double d_n; /**< Default value for the equivalent plastic strain exponent. */

    double d_Sig0; /**< The initial Yield Stress. */

    double default_alpha; /**< Default value of the thermal expansion coefficient. */

    double d_Delta_Time; /**< The time increment (to be used in clad creep). */

    double d_constitutiveMatrix[6]
                               [6]; /**< Stores the 6x6 constitutive matrix for each gauss point. */

    unsigned int d_gaussPtCnt; /**< A counter that keeps track of the number of the gauss points. */

    unsigned int
        Total_Gauss_Point; /**< Total how many gauss points are there in this simulation. */

    unsigned int Plastic_Gauss_Point; /**< How many gauss points have reached plasticity at the
                                         current stage. */

    std::vector<double>
        d_E; /**< A gauss point vector containing the Young's Moduli at each gauss point. */

    std::vector<double>
        d_Nu; /**< A gauss point vector containing the Poisson's Ratio at each gauss point. */

    std::vector<double> d_alpha; /**< A gauss point vector containing the Thermal Expansion
                                    Coefficient at each gauss point. */

    std::vector<double> d_StrengthCoeff; /**< A gauss point vector containing the linear strain
                                            hardening coefficient at each gauss
                                            point. */

    std::vector<double> d_PlasticExponent; /**< A gauss point vector containing the equivalent
                                              plastic strain exponent
                                              at each gauss point. */

    std::vector<double>
        d_EquilibriumStress; /**< A gauss point vector which stores the equilibrium stress. */

    std::vector<double>
        d_EquilibriumStrain; /**< A gauss point vector which stores the equilibrium strain. */

    std::vector<double>
        d_EquilibriumYieldStress; /**< A gauss point vector which stores the yield stress at the
                                     previous equilibrium configuration. */

    std::vector<double>
        d_EquilibriumEffectivePlasticStrain; /**< A gauss point vector which stores the equivalent
                                                plastic strain at the previous equilibrium
                                                configuration. */

    std::vector<double> d_EquilibriumCreepStrain; /**< A gauss point vector which stores the creep
                                                     strain (scalar) at
                                                     the previous equilibrium configuration. */

    std::vector<double> d_EquilibriumTemperature; /**< A gauss point vector which stores the
                                                     equilibrium temperature. */

    std::vector<double> d_EquilibriumThermalStrain; /**< A gauss point vector containing the
                                                       equilibrium Thermal Strain at each gauss
                                                       point. */

    std::vector<double>
        d_Lambda; /**< A gauss point vector which contains a plasticity parameter at each gauss
                     point. */

    std::vector<int> d_ElPl; /**< A gauss point vector which contains the information whether a
                                particular gauss has
                                reached plastic state or still in the elastic regime. */

    std::vector<double>
        d_tmp1Stress; /**< A gauss point vector which stores the stress at the current step. */

    std::vector<double>
        d_tmp1Strain; /**< A gauss point vector which stores the strain at the current step. */

    std::vector<double> d_tmp1YieldStress; /**< A gauss point vector which stores the yield stress
                                              at the current configuration. */

    std::vector<double>
        d_tmp1EffectivePlasticStrain; /**< A gauss point vector which stores the equivalent plastic
                                         strain at the current configuration. */

    std::vector<double>
        d_tmp1CreepStrain; /**< A gauss point vector which stores the creep strain (scalar) at the
                              current configuration. */

    std::vector<double>
        d_tmp1Temperature; /**< A gauss point vector which stores the temperature at the current
                              step. */

    std::vector<double> d_tmp1ThermalStrain; /**< A gauss point vector containing the current
                                                Thermal Strain at each gauss point. */

    std::vector<double>
        d_tmp2Stress; /**< A gauss point vector which stores the stress at the current step while
                         jacobian is being calculated. */

    std::vector<double> d_tmp2YieldStress; /**< A gauss point vector which stores the yield stress
                                              at the current step
                                              while jacobian is being calculated. */

    std::vector<double> d_tmp2EffectivePlasticStrain; /**< A gauss point vector which stores the
                                                         equivalent plastic strain at the
                                                         current step while jacobian is being
                                                         calculated. */

    bool d_resetReusesRadialReturn;

    bool d_jacobianReusesRadialReturn;

    bool d_Is_Init_Called; /**< Checks whether init has been called already or not. */

    bool d_UseThermalStrain; /**< This variable determines whether to use thermal strain or not. */

    bool d_UseCreepStrain; /**< This variable determines whether to use creep strain or not. */

private:
};
} // namespace Operator
} // namespace AMP

#endif
