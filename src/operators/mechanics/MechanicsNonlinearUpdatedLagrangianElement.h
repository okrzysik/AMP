
#ifndef included_AMP_MechanicsNonlinearUpdatedLagrangianElement
#define included_AMP_MechanicsNonlinearUpdatedLagrangianElement

#include <vector>

#include "utils/shared_ptr.h"

/* AMP files */
#include "MechanicsConstants.h"
#include "MechanicsElement.h"
#include "UpdatedLagrangianUtils.h"

namespace AMP {
namespace Operator {


/**
A class for representing the element level computation performed within a
nonlinear updated lagrangian finite element operator for modelling solid mechanics.
*/
class MechanicsNonlinearUpdatedLagrangianElement : public MechanicsElement {
public:
    /**
      This is primarily for use with the updateMaterialModel() function.
     */
    enum MaterialUpdateType { RESET, JACOBIAN };

    //! Constructor.
    explicit MechanicsNonlinearUpdatedLagrangianElement(
        const AMP::shared_ptr<ElementOperationParameters> &params )
        : MechanicsElement( params ), d_elementOutputVector( NULL )
    {
        d_JxW = &( d_fe->get_JxW() );

        d_dphi = &( d_fe->get_dphi() );

        d_phi = &( d_fe->get_phi() );

        d_xyz = &( d_fe->get_xyz() );

        d_onePointShearIntegration = false;

        d_gaussPtCnt = 0;
    }

    //! Destructor.
    virtual ~MechanicsNonlinearUpdatedLagrangianElement() {}

    /**
      This function is used by MechanicsNonlinearFEOperator to pass the address
      of the element Input and Output vector to this class.
      @param [in] elementInputVectors Element input vector
      @param [in] elementInputVectors_pre Element input vector for previous time
      @param [in] elementOutputVector Element residual vector
     */
    void setElementVectors( const std::vector<std::vector<double>> &elementInputVectors,
                            const std::vector<std::vector<double>> &elementInputVectors_pre,
                            std::vector<double> &elementOutputVector )
    {
        d_elementInputVectors     = elementInputVectors;
        d_elementInputVectors_pre = elementInputVectors_pre;
        d_elementOutputVector     = &( elementOutputVector );
    }

    /**
      Element residual vector computation.
     */
    void apply()
    {
        if ( d_useReducedIntegration ) {
            apply_Reduced();
        }
        else {
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
      @param [in] elementInputVectors_pre Field (Displacement, Temperature, Burnup etc) values at
      the
      nodes of the current element at the previous time-step.
     */
    template <MaterialUpdateType type>
    void updateMaterialModel( const std::vector<std::vector<double>> &elementInputVectors,
                              const std::vector<std::vector<double>> &elementInputVectors_pre );

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

    /**
      Computes the deformation gradient at n-th, (n+1)-th and (n+1/2)-th time step.
     */
    void computeDeformationGradient( const std::vector<std::vector<RealGradient>> &dphi,
                                     const std::vector<Point> &xyz,
                                     unsigned int num_nodes,
                                     unsigned int qp,
                                     double F[3][3] );

    /**
      Initializes the reference x, y and z coordinates.
     */
    void initializeReferenceXYZ( std::vector<double> &elementRefXYZ );

    /**
      Assign the reference x, y and z coordinates for the current element.
     */
    void assignReferenceXYZ( std::vector<double> elementRefXYZ )
    {
        d_elementRefXYZ = elementRefXYZ;
    }

    void preNonlinearElementInit();

    void zeroOutGaussPointCount() { d_gaussPtCnt = 0; }

    void resetElementInfo()
    {
        if ( ( d_useJaumannRate == false ) && ( d_useFlanaganTaylorElem == true ) ) {
            d_leftStretchV_n = d_leftStretchV_np1;
            d_rotationR_n    = d_rotationR_np1;
        }
    }

protected:
    template <MaterialUpdateType type>
    void materialModelPreNonlinearElementOperation();

    template <MaterialUpdateType type>
    void materialModelPreNonlinearGaussPointOperation();

    template <MaterialUpdateType type>
    void materialModelNonlinearGaussPointOperation( const std::vector<std::vector<double>> &,
                                                    double[3][3],
                                                    double[3][3] );

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

    std::vector<std::vector<double>> d_elementInputVectors; /**< Element input vectors
                                                                  (Displacement, temperature, burnup
                                                                  etc). */

    std::vector<std::vector<double>>
        d_elementInputVectors_pre; /**< Element input vectors at the previous time step
                                       (Displacement, temperature, burnup etc). */

    std::vector<double> *d_elementOutputVector; /**< Element residual vector */

    bool d_onePointShearIntegration;

    std::vector<double> d_elementRefXYZ;

    std::vector<double> d_leftStretchV_n;

    std::vector<double> d_leftStretchV_np1;

    std::vector<double> d_rotationR_n;

    std::vector<double> d_rotationR_np1;

    int d_gaussPtCnt;

private:
};


template <>
inline void MechanicsNonlinearUpdatedLagrangianElement::materialModelPreNonlinearElementOperation<
    MechanicsNonlinearUpdatedLagrangianElement::RESET>()
{
    d_materialModel->preNonlinearResetElementOperation();
}

template <>
inline void MechanicsNonlinearUpdatedLagrangianElement::materialModelPreNonlinearElementOperation<
    MechanicsNonlinearUpdatedLagrangianElement::JACOBIAN>()
{
    d_materialModel->preNonlinearJacobianElementOperation();
}

template <>
inline void
MechanicsNonlinearUpdatedLagrangianElement::materialModelPreNonlinearGaussPointOperation<
    MechanicsNonlinearUpdatedLagrangianElement::RESET>()
{
    d_materialModel->preNonlinearResetGaussPointOperation();
}

template <>
inline void
MechanicsNonlinearUpdatedLagrangianElement::materialModelPreNonlinearGaussPointOperation<
    MechanicsNonlinearUpdatedLagrangianElement::JACOBIAN>()
{
    d_materialModel->preNonlinearJacobianGaussPointOperation();
}

template <>
inline void MechanicsNonlinearUpdatedLagrangianElement::materialModelPostNonlinearElementOperation<
    MechanicsNonlinearUpdatedLagrangianElement::RESET>()
{
    d_materialModel->postNonlinearResetElementOperation();
}

template <>
inline void MechanicsNonlinearUpdatedLagrangianElement::materialModelPostNonlinearElementOperation<
    MechanicsNonlinearUpdatedLagrangianElement::JACOBIAN>()
{
    d_materialModel->postNonlinearJacobianElementOperation();
}

template <>
inline void
MechanicsNonlinearUpdatedLagrangianElement::materialModelPostNonlinearGaussPointOperation<
    MechanicsNonlinearUpdatedLagrangianElement::RESET>()
{
    d_materialModel->postNonlinearResetGaussPointOperation();
}

template <>
inline void
MechanicsNonlinearUpdatedLagrangianElement::materialModelPostNonlinearGaussPointOperation<
    MechanicsNonlinearUpdatedLagrangianElement::JACOBIAN>()
{
    d_materialModel->postNonlinearJacobianGaussPointOperation();
}

template <>
inline void MechanicsNonlinearUpdatedLagrangianElement::materialModelNonlinearGaussPointOperation<
    MechanicsNonlinearUpdatedLagrangianElement::RESET>(
    const std::vector<std::vector<double>> &strain, double R_n[3][3], double R_np1[3][3] )
{
    d_materialModel->nonlinearResetGaussPointOperation_UL( strain, R_n, R_np1 );
}

template <>
inline void MechanicsNonlinearUpdatedLagrangianElement::materialModelNonlinearGaussPointOperation<
    MechanicsNonlinearUpdatedLagrangianElement::JACOBIAN>(
    const std::vector<std::vector<double>> &strain, double R_n[3][3], double R_np1[3][3] )
{
    d_materialModel->nonlinearJacobianGaussPointOperation_UL( strain, R_n, R_np1 );
}


template <MechanicsNonlinearUpdatedLagrangianElement::MaterialUpdateType type>
void MechanicsNonlinearUpdatedLagrangianElement::updateMaterialModel(
    const std::vector<std::vector<double>> &elementInputVectors,
    const std::vector<std::vector<double>> &elementInputVectors_pre )
{
    const std::vector<std::vector<RealGradient>> &dphi = ( *d_dphi );

    const std::vector<std::vector<Real>> &phi = ( *d_phi );

    std::vector<Point> xyz, xyz_n, xyz_np1, xyz_np1o2;

    d_fe->reinit( d_elem );

    materialModelPreNonlinearElementOperation<type>();

    const unsigned int num_nodes = d_elem->n_nodes();

    xyz.resize( num_nodes );
    xyz_n.resize( num_nodes );
    xyz_np1.resize( num_nodes );
    xyz_np1o2.resize( num_nodes );

    Point p1;
    for ( unsigned int ijk = 0; ijk < num_nodes; ijk++ ) {
        p1       = d_elem->point( ijk );
        xyz[ijk] = p1;
    }

    double currX[8], currY[8], currZ[8], dNdx[8], dNdy[8], dNdz[8], detJ[1], delta_u[8], delta_v[8],
        delta_w[8];
    double x_np1o2[8], y_np1o2[8], z_np1o2[8];
    double rsq3              = ( 1.0 / sqrt( 3.0 ) );
    const double currXi[8]   = { -rsq3, rsq3, -rsq3, rsq3, -rsq3, rsq3, -rsq3, rsq3 };
    const double currEta[8]  = { -rsq3, -rsq3, rsq3, rsq3, -rsq3, -rsq3, rsq3, rsq3 };
    const double currZeta[8] = { -rsq3, -rsq3, -rsq3, -rsq3, rsq3, rsq3, rsq3, rsq3 };
    double Bl_np1_bar[6][24], Bl_center[6][24];

    for ( unsigned int ijk = 0; ijk < num_nodes; ijk++ ) {
        xyz_n[ijk]( 0 ) =
            xyz[ijk]( 0 ) + elementInputVectors_pre[Mechanics::DISPLACEMENT][( 3 * ijk ) + 0];
        xyz_n[ijk]( 1 ) =
            xyz[ijk]( 1 ) + elementInputVectors_pre[Mechanics::DISPLACEMENT][( 3 * ijk ) + 1];
        xyz_n[ijk]( 2 ) =
            xyz[ijk]( 2 ) + elementInputVectors_pre[Mechanics::DISPLACEMENT][( 3 * ijk ) + 2];

        currX[ijk] = xyz_np1[ijk]( 0 ) =
            xyz[ijk]( 0 ) + elementInputVectors[Mechanics::DISPLACEMENT][( 3 * ijk ) + 0];
        currY[ijk] = xyz_np1[ijk]( 1 ) =
            xyz[ijk]( 1 ) + elementInputVectors[Mechanics::DISPLACEMENT][( 3 * ijk ) + 1];
        currZ[ijk] = xyz_np1[ijk]( 2 ) =
            xyz[ijk]( 2 ) + elementInputVectors[Mechanics::DISPLACEMENT][( 3 * ijk ) + 2];

        delta_u[ijk] = elementInputVectors[Mechanics::DISPLACEMENT][( 3 * ijk ) + 0] -
                       elementInputVectors_pre[Mechanics::DISPLACEMENT][( 3 * ijk ) + 0];
        delta_v[ijk] = elementInputVectors[Mechanics::DISPLACEMENT][( 3 * ijk ) + 1] -
                       elementInputVectors_pre[Mechanics::DISPLACEMENT][( 3 * ijk ) + 1];
        delta_w[ijk] = elementInputVectors[Mechanics::DISPLACEMENT][( 3 * ijk ) + 2] -
                       elementInputVectors_pre[Mechanics::DISPLACEMENT][( 3 * ijk ) + 2];

        x_np1o2[ijk] = xyz_np1o2[ijk]( 0 ) = xyz_n[ijk]( 0 ) + ( delta_u[ijk] / 2.0 );
        y_np1o2[ijk] = xyz_np1o2[ijk]( 1 ) = xyz_n[ijk]( 1 ) + ( delta_v[ijk] / 2.0 );
        z_np1o2[ijk] = xyz_np1o2[ijk]( 2 ) = xyz_n[ijk]( 2 ) + ( delta_w[ijk] / 2.0 );
    }

    if ( d_useReducedIntegration && d_useJaumannRate ) {
        double sum_detJ = 0.0;
        for ( unsigned int i = 0; i < 6; i++ ) {
            for ( unsigned int j = 0; j < ( 3 * num_nodes ); j++ ) {
                Bl_np1_bar[i][j] = 0.0;
                Bl_center[i][j]  = 0.0;
            }
        }

        for ( unsigned int qp = 0; qp < d_qrule->n_points(); qp++ ) {
            constructShapeFunctionDerivatives( dNdx,
                                               dNdy,
                                               dNdz,
                                               currX,
                                               currY,
                                               currZ,
                                               currXi[qp],
                                               currEta[qp],
                                               currZeta[qp],
                                               detJ );
            sum_detJ += detJ[0];

            for ( unsigned int i = 0; i < 8; i++ ) {
                Bl_np1_bar[0][( 3 * i ) + 0] += ( dNdx[i] * detJ[0] );
                Bl_np1_bar[1][( 3 * i ) + 0] += ( dNdx[i] * detJ[0] );
                Bl_np1_bar[2][( 3 * i ) + 0] += ( dNdx[i] * detJ[0] );

                Bl_np1_bar[0][( 3 * i ) + 1] += ( dNdy[i] * detJ[0] );
                Bl_np1_bar[1][( 3 * i ) + 1] += ( dNdy[i] * detJ[0] );
                Bl_np1_bar[2][( 3 * i ) + 1] += ( dNdy[i] * detJ[0] );

                Bl_np1_bar[0][( 3 * i ) + 2] += ( dNdz[i] * detJ[0] );
                Bl_np1_bar[1][( 3 * i ) + 2] += ( dNdz[i] * detJ[0] );
                Bl_np1_bar[2][( 3 * i ) + 2] += ( dNdz[i] * detJ[0] );
            }
        }

        // std::cout<<"sum_detJ="<<sum_detJ<<std::endl;
        double one3TimesSumDetJ = 1.0 / ( 3.0 * sum_detJ );
        for ( unsigned int i = 0; i < 6; i++ ) {
            for ( unsigned int j = 0; j < ( 3 * num_nodes ); j++ ) {
                Bl_np1_bar[i][j] = Bl_np1_bar[i][j] * one3TimesSumDetJ;
                // std::cout<<"Bl_np1_bar["<<i<<"]["<<j<<"]="<<Bl_np1_bar[i][j]<<std::endl;
            }
        }

        constructShapeFunctionDerivatives(
            dNdx, dNdy, dNdz, currX, currY, currZ, 0.0, 0.0, 0.0, detJ );
        for ( unsigned int qp = 0; qp < d_qrule->n_points(); qp++ ) {
            Bl_center[0][( 3 * qp ) + 0] += ( dNdx[qp] );
            Bl_center[1][( 3 * qp ) + 1] += ( dNdy[qp] );
            Bl_center[2][( 3 * qp ) + 2] += ( dNdz[qp] );
            Bl_center[3][( 3 * qp ) + 1] += ( dNdz[qp] );
            Bl_center[3][( 3 * qp ) + 2] += ( dNdy[qp] );
            Bl_center[4][( 3 * qp ) + 0] += ( dNdz[qp] );
            Bl_center[4][( 3 * qp ) + 2] += ( dNdx[qp] );
            Bl_center[5][( 3 * qp ) + 0] += ( dNdy[qp] );
            Bl_center[5][( 3 * qp ) + 1] += ( dNdx[qp] );
        }
    }

    for ( unsigned int qp = 0; qp < d_qrule->n_points(); qp++ ) {
        materialModelPreNonlinearGaussPointOperation<type>();

        double Identity[3][3];

        // Constructing an identity matrix.
        for ( unsigned int i = 0; i < 3; i++ ) {
            for ( unsigned int j = 0; j < 3; j++ ) {
                Identity[i][j] = 0.0;
            }
            Identity[i][i] = 1.0;
        }

        double Bl_np1[6][24], spin_np1[3][3], el_np1[6];
        double R_n[3][3], R_np1[3][3];
        double e_np1o2_tilda_rotated[3][3];

        if ( d_useJaumannRate == false ) {
            double R_np1o2[3][3];
            if ( d_useFlanaganTaylorElem == false ) {
                double U_n[3][3], U_np1[3][3], U_np1o2[3][3];
                double F_n[3][3], F_np1[3][3], F_np1o2[3][3];
                // The deformation gradients are computed in the next three lines.
                computeDeformationGradient( dphi, xyz_n, num_nodes, qp, F_n );
                computeDeformationGradient( dphi, xyz_np1, num_nodes, qp, F_np1 );
                computeDeformationGradient( dphi, xyz_np1o2, num_nodes, qp, F_np1o2 );

                // Polar decomposition (F=RU) of the deformation gradient is conducted here.
                polarDecompositionFeqRU_Simo( F_n, R_n, U_n );
                polarDecompositionFeqRU_Simo( F_np1, R_np1, U_np1 );
                polarDecompositionFeqRU_Simo( F_np1o2, R_np1o2, U_np1o2 );
            }

            // Gradient of the incremental displacement with respect to the np1o2 configuration.
            double dN_dxnp1o2[8], dN_dynp1o2[8], dN_dznp1o2[8], detJ_np1o2[1], d_np1o2[3][3],
                d_np1o2_temp[3][3];
            constructShapeFunctionDerivatives( dN_dxnp1o2,
                                               dN_dynp1o2,
                                               dN_dznp1o2,
                                               x_np1o2,
                                               y_np1o2,
                                               z_np1o2,
                                               currXi[qp],
                                               currEta[qp],
                                               currZeta[qp],
                                               detJ_np1o2 );

            // Calculate the rate of deformation tensor with respect to the np1o2 configuration.
            computeGradient( dN_dxnp1o2,
                             dN_dynp1o2,
                             dN_dznp1o2,
                             delta_u,
                             delta_v,
                             delta_w,
                             num_nodes,
                             d_np1o2_temp );
            for ( int i = 0; i < 3; i++ ) {
                for ( int j = 0; j < 3; j++ ) {
                    d_np1o2[i][j] = 0.5 * ( d_np1o2_temp[i][j] + d_np1o2_temp[j][i] );
                }
            }

            // Rotate the rate of deformation to the unrotated configuration.
            pullbackCorotational( R_np1o2, d_np1o2, e_np1o2_tilda_rotated );
        }

        if ( d_useJaumannRate == true ) {
            // Calculate the derivatives of the shape functions at the current coordinate.
            constructShapeFunctionDerivatives( dNdx,
                                               dNdy,
                                               dNdz,
                                               currX,
                                               currY,
                                               currZ,
                                               currXi[qp],
                                               currEta[qp],
                                               currZeta[qp],
                                               detJ );

            double Bl_dil[6][24];
            for ( int i = 0; i < 6; i++ ) {
                el_np1[i] = 0.0;
                for ( int j = 0; j < 24; j++ ) {
                    Bl_np1[i][j] = 0.0;
                    if ( d_useReducedIntegration ) {
                        Bl_dil[i][j] = 0.0;
                    }
                }
            }

            for ( unsigned int i = 0; i < num_nodes; i++ ) {
                Bl_np1[0][( 3 * i ) + 0] = dNdx[i];
                Bl_np1[1][( 3 * i ) + 1] = dNdy[i];
                Bl_np1[2][( 3 * i ) + 2] = dNdz[i];
                Bl_np1[3][( 3 * i ) + 1] = dNdz[i];
                Bl_np1[3][( 3 * i ) + 2] = dNdy[i];
                Bl_np1[4][( 3 * i ) + 0] = dNdz[i];
                Bl_np1[4][( 3 * i ) + 2] = dNdx[i];
                Bl_np1[5][( 3 * i ) + 0] = dNdy[i];
                Bl_np1[5][( 3 * i ) + 1] = dNdx[i];

                if ( d_useReducedIntegration ) {
                    double one3 = 1.0 / 3.0;
                    for ( int j = 0; j < 3; j++ ) {
                        Bl_dil[j][( 3 * i ) + 0] = dNdx[i] * one3;
                        Bl_dil[j][( 3 * i ) + 1] = dNdy[i] * one3;
                        Bl_dil[j][( 3 * i ) + 2] = dNdz[i] * one3;
                    }
                }
            }

            // Calculate the spin tensor for the jaumann rate.
            double d_np1[3][3];
            computeGradient( dNdx, dNdy, dNdz, delta_u, delta_v, delta_w, num_nodes, d_np1 );
            for ( int i = 0; i < 3; i++ ) {
                for ( int j = 0; j < 3; j++ ) {
                    spin_np1[i][j] = 0.5 * ( d_np1[i][j] - d_np1[j][i] );
                }
            }

            if ( d_useReducedIntegration ) {
                for ( int i = 0; i < 6; i++ ) {
                    for ( unsigned int j = 0; j < ( 3 * num_nodes ); j++ ) {
                        Bl_np1[i][j] = Bl_np1[i][j] - Bl_dil[i][j] + Bl_np1_bar[i][j];
                        // std::cout<<"Bl_np1["<<i<<"]["<<j<<"]="<<Bl_np1[i][j]<<std::endl;
                    }
                }
            }

            if ( d_onePointShearIntegration ) {
                for ( int i = 3; i < 6; i++ ) {
                    for ( unsigned int j = 0; j < ( 3 * num_nodes ); j++ ) {
                        Bl_np1[i][j] = Bl_center[i][j];
                    }
                }
            }

            for ( int i = 0; i < 6; i++ ) {
                for ( int j = 0; j < 8; j++ ) {
                    el_np1[i] += ( Bl_np1[i][( 3 * j ) + 0] * delta_u[j] );
                    el_np1[i] += ( Bl_np1[i][( 3 * j ) + 1] * delta_v[j] );
                    el_np1[i] += ( Bl_np1[i][( 3 * j ) + 2] * delta_w[j] );
                }
            }
        }

        /* Compute Strain From Given Displacement */

        std::vector<std::vector<double>> fieldsAtGaussPt( Mechanics::TOTAL_NUMBER_OF_VARIABLES );

        // Strain
        if ( d_useJaumannRate == true ) {
            fieldsAtGaussPt[Mechanics::DISPLACEMENT].push_back( el_np1[0] );
            fieldsAtGaussPt[Mechanics::DISPLACEMENT].push_back( el_np1[1] );
            fieldsAtGaussPt[Mechanics::DISPLACEMENT].push_back( el_np1[2] );
            fieldsAtGaussPt[Mechanics::DISPLACEMENT].push_back( el_np1[3] );
            fieldsAtGaussPt[Mechanics::DISPLACEMENT].push_back( el_np1[4] );
            fieldsAtGaussPt[Mechanics::DISPLACEMENT].push_back( el_np1[5] );
        }

        if ( d_useJaumannRate == false ) {
            fieldsAtGaussPt[Mechanics::DISPLACEMENT].push_back( e_np1o2_tilda_rotated[0][0] );
            fieldsAtGaussPt[Mechanics::DISPLACEMENT].push_back( e_np1o2_tilda_rotated[1][1] );
            fieldsAtGaussPt[Mechanics::DISPLACEMENT].push_back( e_np1o2_tilda_rotated[2][2] );
            fieldsAtGaussPt[Mechanics::DISPLACEMENT].push_back( 2.0 * e_np1o2_tilda_rotated[1][2] );
            fieldsAtGaussPt[Mechanics::DISPLACEMENT].push_back( 2.0 * e_np1o2_tilda_rotated[0][2] );
            fieldsAtGaussPt[Mechanics::DISPLACEMENT].push_back( 2.0 * e_np1o2_tilda_rotated[0][1] );
        }

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

        if ( d_useJaumannRate == true ) {
            materialModelNonlinearGaussPointOperation<type>( fieldsAtGaussPt, spin_np1, Identity );
        }
        else {
            materialModelNonlinearGaussPointOperation<type>( fieldsAtGaussPt, R_n, R_np1 );
        }

        materialModelPostNonlinearGaussPointOperation<type>();
    } // end for qp

    materialModelPostNonlinearElementOperation<type>();
}
}
}

#endif
