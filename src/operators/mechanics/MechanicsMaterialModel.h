#ifndef included_AMP_MechanicsMaterialModel
#define included_AMP_MechanicsMaterialModel

#include <cstring>

#include "AMP/materials/Material.h"
#include "AMP/operators/ElementPhysicsModel.h"
#include "AMP/operators/mechanics/UpdatedLagrangianUtils.h"
#include <memory>

namespace AMP {
namespace Operator {

typedef ElementPhysicsModelParameters MechanicsMaterialModelParameters;

/**
 * An abstract base class for representing the mechanics material models. The
 derived classes can represent both linear and nonlinear material behavior.
 * Derived classes must implement the function named getConstitutiveMatrix,
 which provides the ConstitutiveMatrix that is used
 * in the construction of the Jacobian. Non-linear material models must
 * implement the function getInternalStress. Linear material models that
 use an implicit method must implement the getInternalStress function too.
 */
class MechanicsMaterialModel : public ElementPhysicsModel
{
public:
    /**
     *  Constructor. This reads the value for the key USE_MATERIALS_LIBRARY (false by default)
     from the database object contained in the parameter object, params. This key specifies
     whether or not the AMP::materials interface is used in this model.
     */
    explicit MechanicsMaterialModel( std::shared_ptr<MechanicsMaterialModelParameters> params )
        : ElementPhysicsModel( params )
    {
        d_useMaterialsLibrary = params->d_db->getWithDefault( "USE_MATERIALS_LIBRARY", false );

        d_useUpdatedLagrangian = params->d_db->getWithDefault( "USE_UPDATED_LAGRANGIAN", false );

        d_useJaumannRate = params->d_db->getWithDefault( "USE_JAUMANN_RATE", false );

        d_useContinuumTangent = params->d_db->getWithDefault( "USE_CONTINUUM_TANGENT", false );

        d_checkCladOrPellet = params->d_db->getWithDefault( "IS_IT_CLAD", false );

        if ( d_useMaterialsLibrary == true ) {
            AMP_INSIST( ( params->d_db->keyExists( "Material" ) ), "Key ''Material'' is missing!" );
            std::string matname = params->d_db->getString( "Material" );
            d_material =
                AMP::voodoo::Factory<AMP::Materials::Material>::instance().create( matname );
        }

        d_currentTime  = 0;
        d_previousTime = 0;
    }

    /**
     * Destructor.
     */
    virtual ~MechanicsMaterialModel() {}

    /**
     * Calculates the constitutive matrix for the material model. This matrix
     * is used for the construction of the Jacobian during the solve process.
     */
    virtual void getConstitutiveMatrix( double *& ) {}

    /**
     * Calculates the constitutive matrix for the material model in Updated Lagrangian. This matrix
     * is used for the construction of the Jacobian during the solve process.
     */
    virtual void getConstitutiveMatrixUpdatedLagrangian( double[6][6], double[3][3] ) {}

    /**
     * Returns the 6x1 stress vector at the current gauss point.
     * Used in UpdatedLagrangian calculations.
     */
    virtual void getStressForUpdatedLagrangian( double[6] ) {}

    /**
     * Given a strain, the stress state is calculated in ths function. It is
     * necessary for non-linear material models or linear models with
     * implicit solver.
     */
    virtual void getInternalStress( const std::vector<std::vector<double>> &, double *& ) {}

    virtual void getInternalStress_UL(
        const std::vector<std::vector<double>> &, double *&, double[3][3], double[3][3], double )
    {
    }

    /**
     * Used to print the effective stress at any point of the simulation.
     */
    virtual void getEffectiveStress( double *& ) {}

    /**
     * Used to print the equivalent plastic or creep
     * or thermal strain at any point of the simulation.
     */
    virtual void getEquivalentStrain( double *& ) {}

    /**
     * Used for linear material models if the problem is being solved in an
     * explicit method. This function has not been implemented, because all
     * the linear material model problems are being solved in implicit way.
     */
    virtual void getExternalStress( double *& ) {}

    /*
       A bunch of virtual functions have been implemented to allow the
       programmer to initialize, reset or update any of the variables within
       the material model from outside the gauss point calculation.
       Most of these functions are not implemented. The are added whenever it
       is required.
       */

    // For LinearOperator's Reset/Init

    virtual void preLinearAssembly() {}

    virtual void postLinearAssembly() {}

    virtual void preLinearElementOperation() {}

    virtual void postLinearElementOperation() {}

    virtual void preLinearGaussPointOperation() {}

    virtual void postLinearGaussPointOperation() {}

    // For NonlinearOperator's Init

    virtual void preNonlinearInit( bool, bool ) {}

    virtual void postNonlinearInit() {}

    virtual void preNonlinearInitElementOperation() {}

    virtual void postNonlinearInitElementOperation() {}

    virtual void preNonlinearInitGaussPointOperation() {}

    virtual void postNonlinearInitGaussPointOperation() {}

    /**
     * Initializes all the variables with zero, except the temperature
     * variable which has some non-zero value initially (something like
     * room temperature). The input argument is the initial_Temperature.
     */
    virtual void nonlinearInitGaussPointOperation( double ) {}

    // For NonlinearOperator's Apply

    virtual void preNonlinearAssembly() {}

    virtual void postNonlinearAssembly() {}

    virtual void preNonlinearAssemblyElementOperation() {}

    virtual void postNonlinearAssemblyElementOperation() {}

    virtual void preNonlinearAssemblyGaussPointOperation() {}

    virtual void postNonlinearAssemblyGaussPointOperation() {}

    // For NonlinearOperator's Reset

    virtual void globalReset() {}

    virtual void preNonlinearReset() {}

    virtual void postNonlinearReset() {}

    virtual void preNonlinearResetElementOperation() {}

    virtual void postNonlinearResetElementOperation() {}

    virtual void preNonlinearResetGaussPointOperation() {}

    virtual void postNonlinearResetGaussPointOperation() {}

    /**
     * In the implicit solution technique, once the solver converges,
     * the previous equilibrium values are replaced by the current converged
     * values in this function. The input is a vector of all the variables at
     * that particular gauss point.
     */
    virtual void nonlinearResetGaussPointOperation( const std::vector<std::vector<double>> & ) {}

    virtual void nonlinearResetGaussPointOperation_UL( const std::vector<std::vector<double>> &,
                                                       double[3][3],
                                                       double[3][3] )
    {
    }

    // For NonlinearOperator's GetJacobianParameters

    virtual void preNonlinearJacobian() {}

    virtual void postNonlinearJacobian() {}

    virtual void preNonlinearJacobianElementOperation() {}

    virtual void postNonlinearJacobianElementOperation() {}

    virtual void preNonlinearJacobianGaussPointOperation() {}

    virtual void postNonlinearJacobianGaussPointOperation() {}

    virtual void nonlinearJacobianGaussPointOperation( const std::vector<std::vector<double>> & ) {}

    virtual void nonlinearJacobianGaussPointOperation_UL( const std::vector<std::vector<double>> &,
                                                          double[3][3],
                                                          double[3][3] )
    {
    }

    void updateTime( double currTime )
    {
        d_previousTime = d_currentTime;
        d_currentTime  = currTime;
    }

    std::shared_ptr<AMP::Materials::Material> getMaterial() { return d_material; }

protected:
    double d_currentTime; /**< The time at present. */

    double d_previousTime; /**< Time at the previous step. */

    bool d_useMaterialsLibrary; /**< A flag that is true if the AMP::materials
                                  library is used in this model and false otherwise.  */

    bool d_useUpdatedLagrangian; /**< Flag to check whether to use updated lagrangian or not. */

    bool d_useJaumannRate; /**< Flag to check whether to use Jaumann rate in updated lagrangian or
                              not. */

    bool d_useContinuumTangent; /**< Flag to check whether to use Continuum tangent is elasto
                                   plasticity or not. */

    bool d_checkCladOrPellet;

    std::shared_ptr<AMP::Materials::Material>
        d_material; /**< Shared pointer to the materials object. */

private:
};
} // namespace Operator
} // namespace AMP

#endif
