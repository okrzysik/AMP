
#ifndef included_AMP_MechanicsNonlinearElement
#define included_AMP_MechanicsNonlinearElement

#include <vector>

#include "utils/shared_ptr.h"

/* AMP files */
#include "MechanicsConstants.h"
#include "MechanicsElement.h"

namespace AMP {
namespace Operator {


/**
  A class for representing the element level computation performed within a
  nonlinear finite element operator for modelling solid mechanics.
 */
class MechanicsNonlinearElement : public MechanicsElement
{
public:
    /**
      This is primarily for use with the updateMaterialModel() function.
      */
    enum MaterialUpdateType { RESET, JACOBIAN };

    //! Constructor.
    explicit MechanicsNonlinearElement( const AMP::shared_ptr<ElementOperationParameters> &params )
        : MechanicsElement( params ), d_elementOutputVector( nullptr )
    {
        d_JxW = &( d_fe->get_JxW() );

        d_dphi = &( d_fe->get_dphi() );

        d_phi = &( d_fe->get_phi() );

        d_xyz = &( d_fe->get_xyz() );

        AMP_INSIST( ( d_useJaumannRate == false ),
                    "Jaumann rate with small strain does not make any sense." );
    }

    //! Destructor.
    virtual ~MechanicsNonlinearElement() {}

    /**
      This function is used by MechanicsNonlinearFEOperator to pass the address
      of the element Input and Output vector to this class.
      @param [in] elementInputVectors Element input vector
      @param [in] elementOutputVector Element residual vector
     */
    void setElementVectors( const std::vector<std::vector<double>> &elementInputVectors,
                            std::vector<double> &elementOutputVector )
    {
        d_elementInputVectors = elementInputVectors;
        d_elementOutputVector = &( elementOutputVector );
    }

    /**
      Element residual vector computation.
     */
    void apply() override {
      if (d_useReducedIntegration) {
        apply_Reduced();
      } else {
        apply_Normal();
      }
    }

    /**
      This function is used to initialize the material model at all Gauss points within the current
      element.
      Some material models require the reference values for some fields such as temperature and so
      these values are also passed to the material model through this function.
      Currently, the only field required by the mechanics material models during initialization is
      the reference
      temperature.
      @param [in] initTempVector Reference temperature at the nodes of the current element. This can
      be
      empty if the material model does not require a reference temperature.
     */
    void initMaterialModel( const std::vector<double> &initTempVector );

    /**
      This function is used to update the equilibrium values of stress, strain, temperature, burnup
      etc that are
      stored in the material model class. This is typically required at the end of each loading
      step. This function is
      typically
      called from within the MechanicsNonlinearFEOperator's reset function, which is evaluated at
      the end of each
      loading step.
      This function also gets called from within the MechanicsNonlinearFEOperator's
      getJacobianParameters
      function if the jacobian is evaluated at a state (displacement, temperature etc) different
      from
      that of the last call to MechanicsNonlinearFEOperator's apply function.
      @param [in] elementInputVectors Field (Displacement, Temperature, Burnup etc) values at the
      nodes of the current element.
     */
    template <MaterialUpdateType type>
    void updateMaterialModel( const std::vector<std::vector<double>> &elementInputVectors );

    /**
      Writes the stess and strain values at the Gauss points within the current element to the file
      The 6 components of stress and strain at each Gauss point are arranged in the order:
      xx, yy, zz, yz, xz and  xy.
      @param [in] fp File pointer
      @param [in] inputVec Input vector (Displacement, Temperature, Burnup etc) at the nodes of the
      current element.
     */
    void printStressAndStrain( FILE *fp, const std::vector<std::vector<double>> &inputVec );

    /**
      Computes the stress and strain values at the Gauss points within the current element
      The 6 components of stress and strain at each Gauss point are arranged in the order:
      xx, yy, zz, yz, xz and  xy.
      @param [in] inputVec Input vector (Displacement, Temperature, Burnup etc) at the nodes of the
      current element.
      @param [out] stressVec Stresses at the Gauss points of the current element.
      @param [out] strainVec Strains at the Gauss points of the current element.
     */
    void computeStressAndStrain( const std::vector<std::vector<double>> &inputVec,
                                 std::vector<double> &stressVec,
                                 std::vector<double> &strainVec );

protected:
    template <MaterialUpdateType type>
    void materialModelPreNonlinearElementOperation();

    template <MaterialUpdateType type>
    void materialModelPreNonlinearGaussPointOperation();

    template <MaterialUpdateType type>
    void materialModelNonlinearGaussPointOperation( const std::vector<std::vector<double>> & );

    template <MaterialUpdateType type>
    void materialModelPostNonlinearGaussPointOperation();

    template <MaterialUpdateType type>
    void materialModelPostNonlinearElementOperation();

    /**
      Element residual vector computation using normal integration scheme.
     */
    void apply_Normal();

    /**
      Element residual vector computation using reduced integration scheme.
     */
    void apply_Reduced();

    const std::vector<Real> *d_JxW; /**< Product of the determinant of Jacobian and the quadrature
                                    weight at the Gauss points in the current element. */

    const std::vector<std::vector<RealGradient>>
        *d_dphi; /**< Spatial Derivatives of the shape functions at
                  the Gauss points in the current element. */

    const std::vector<std::vector<Real>> *d_phi; /**< Shape functions at
                                                        the Gauss points in the current element. */

    const std::vector<Point> *d_xyz; /**< Locations of the Gauss points in the current element. */

    std::vector<std::vector<double>>
        d_elementInputVectors; /**< Element input vectors
                                      (Displacement, temperature, burnup
                                      etc). */

    std::vector<double> *d_elementOutputVector; /**< Element residual vector */

private:
};


template <>
inline void MechanicsNonlinearElement::materialModelPreNonlinearElementOperation<
    MechanicsNonlinearElement::RESET>()
{
    d_materialModel->preNonlinearResetElementOperation();
}

template <>
inline void MechanicsNonlinearElement::materialModelPreNonlinearElementOperation<
    MechanicsNonlinearElement::JACOBIAN>()
{
    d_materialModel->preNonlinearJacobianElementOperation();
}

template <>
inline void MechanicsNonlinearElement::materialModelPreNonlinearGaussPointOperation<
    MechanicsNonlinearElement::RESET>()
{
    d_materialModel->preNonlinearResetGaussPointOperation();
}

template <>
inline void MechanicsNonlinearElement::materialModelPreNonlinearGaussPointOperation<
    MechanicsNonlinearElement::JACOBIAN>()
{
    d_materialModel->preNonlinearJacobianGaussPointOperation();
}

template <>
inline void MechanicsNonlinearElement::materialModelPostNonlinearElementOperation<
    MechanicsNonlinearElement::RESET>()
{
    d_materialModel->postNonlinearResetElementOperation();
}

template <>
inline void MechanicsNonlinearElement::materialModelPostNonlinearElementOperation<
    MechanicsNonlinearElement::JACOBIAN>()
{
    d_materialModel->postNonlinearJacobianElementOperation();
}

template <>
inline void MechanicsNonlinearElement::materialModelPostNonlinearGaussPointOperation<
    MechanicsNonlinearElement::RESET>()
{
    d_materialModel->postNonlinearResetGaussPointOperation();
}

template <>
inline void MechanicsNonlinearElement::materialModelPostNonlinearGaussPointOperation<
    MechanicsNonlinearElement::JACOBIAN>()
{
    d_materialModel->postNonlinearJacobianGaussPointOperation();
}

template <>
inline void MechanicsNonlinearElement::materialModelNonlinearGaussPointOperation<
    MechanicsNonlinearElement::RESET>( const std::vector<std::vector<double>> &strain )
{
    d_materialModel->nonlinearResetGaussPointOperation( strain );
}

template <>
inline void MechanicsNonlinearElement::materialModelNonlinearGaussPointOperation<
    MechanicsNonlinearElement::JACOBIAN>( const std::vector<std::vector<double>> &strain )
{
    d_materialModel->nonlinearJacobianGaussPointOperation( strain );
}

template <MechanicsNonlinearElement::MaterialUpdateType type>
void MechanicsNonlinearElement::updateMaterialModel(
    const std::vector<std::vector<double>> &elementInputVectors )
{
    const std::vector<std::vector<RealGradient>> &dphi = ( *d_dphi );

    const std::vector<std::vector<Real>> &phi = ( *d_phi );

    d_fe->reinit( d_elem );

    materialModelPreNonlinearElementOperation<type>();

    const unsigned int num_nodes = d_elem->n_nodes();

    for ( unsigned int qp = 0; qp < d_qrule->n_points(); qp++ ) {
        materialModelPreNonlinearGaussPointOperation<type>();

        /* Compute Strain From Given Displacement */

        double dudx = 0;
        double dudy = 0;
        double dudz = 0;
        double dvdx = 0;
        double dvdy = 0;
        double dvdz = 0;
        double dwdx = 0;
        double dwdy = 0;
        double dwdz = 0;

        for ( unsigned int k = 0; k < num_nodes; k++ ) {
            dudx +=
                ( elementInputVectors[Mechanics::DISPLACEMENT][( 3 * k ) + 0] * dphi[k][qp]( 0 ) );
            dudy +=
                ( elementInputVectors[Mechanics::DISPLACEMENT][( 3 * k ) + 0] * dphi[k][qp]( 1 ) );
            dudz +=
                ( elementInputVectors[Mechanics::DISPLACEMENT][( 3 * k ) + 0] * dphi[k][qp]( 2 ) );

            dvdx +=
                ( elementInputVectors[Mechanics::DISPLACEMENT][( 3 * k ) + 1] * dphi[k][qp]( 0 ) );
            dvdy +=
                ( elementInputVectors[Mechanics::DISPLACEMENT][( 3 * k ) + 1] * dphi[k][qp]( 1 ) );
            dvdz +=
                ( elementInputVectors[Mechanics::DISPLACEMENT][( 3 * k ) + 1] * dphi[k][qp]( 2 ) );

            dwdx +=
                ( elementInputVectors[Mechanics::DISPLACEMENT][( 3 * k ) + 2] * dphi[k][qp]( 0 ) );
            dwdy +=
                ( elementInputVectors[Mechanics::DISPLACEMENT][( 3 * k ) + 2] * dphi[k][qp]( 1 ) );
            dwdz +=
                ( elementInputVectors[Mechanics::DISPLACEMENT][( 3 * k ) + 2] * dphi[k][qp]( 2 ) );
        } // end for k

        std::vector<std::vector<double>> fieldsAtGaussPt( Mechanics::TOTAL_NUMBER_OF_VARIABLES );

        // Strain
        fieldsAtGaussPt[Mechanics::DISPLACEMENT].push_back( dudx );
        fieldsAtGaussPt[Mechanics::DISPLACEMENT].push_back( dvdy );
        fieldsAtGaussPt[Mechanics::DISPLACEMENT].push_back( dwdz );
        fieldsAtGaussPt[Mechanics::DISPLACEMENT].push_back( 0.5 * ( dvdz + dwdy ) );
        fieldsAtGaussPt[Mechanics::DISPLACEMENT].push_back( 0.5 * ( dudz + dwdx ) );
        fieldsAtGaussPt[Mechanics::DISPLACEMENT].push_back( 0.5 * ( dudy + dvdx ) );

        if ( !( elementInputVectors[Mechanics::TEMPERATURE].empty() ) ) {
            double valAtGaussPt = 0;
            for ( unsigned int k = 0; k < num_nodes; k++ ) {
                valAtGaussPt += ( elementInputVectors[Mechanics::TEMPERATURE][k] * phi[k][qp] );
            } // end for k
            fieldsAtGaussPt[Mechanics::TEMPERATURE].push_back( valAtGaussPt );
        }

        if ( !( elementInputVectors[Mechanics::BURNUP].empty() ) ) {
            double valAtGaussPt = 0;
            for ( unsigned int k = 0; k < num_nodes; k++ ) {
                valAtGaussPt += ( elementInputVectors[Mechanics::BURNUP][k] * phi[k][qp] );
            } // end for k
            fieldsAtGaussPt[Mechanics::BURNUP].push_back( valAtGaussPt );
        }

        if ( !( elementInputVectors[Mechanics::OXYGEN_CONCENTRATION].empty() ) ) {
            double valAtGaussPt = 0;
            for ( unsigned int k = 0; k < num_nodes; k++ ) {
                valAtGaussPt +=
                    ( elementInputVectors[Mechanics::OXYGEN_CONCENTRATION][k] * phi[k][qp] );
            } // end for k
            fieldsAtGaussPt[Mechanics::OXYGEN_CONCENTRATION].push_back( valAtGaussPt );
        }

        if ( !( elementInputVectors[Mechanics::LHGR].empty() ) ) {
            double valAtGaussPt = 0;
            for ( unsigned int k = 0; k < num_nodes; k++ ) {
                valAtGaussPt += ( elementInputVectors[Mechanics::LHGR][k] * phi[k][qp] );
            } // end for k
            fieldsAtGaussPt[Mechanics::LHGR].push_back( valAtGaussPt );
        }

        materialModelNonlinearGaussPointOperation<type>( fieldsAtGaussPt );

        materialModelPostNonlinearGaussPointOperation<type>();
    } // end for qp

    materialModelPostNonlinearElementOperation<type>();
}
}
}

#endif
