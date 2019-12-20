
#ifndef included_AMP_MechanicsLinearElement
#define included_AMP_MechanicsLinearElement

#include <vector>

#include <memory>

/* AMP files */
#include "MechanicsElement.h"

namespace AMP {
namespace Operator {

/**
  A class for representing the element level computation performed within a
  linear finite element operator for modelling solid mechanics.
  The linear operator could either be a linear elasticity operator or it could be
  the jacobian of a nonlinear elasticity/elasto-plasticity operator.
 */
class MechanicsLinearElement : public MechanicsElement
{
public:
    //! Constructor.
    explicit MechanicsLinearElement( const std::shared_ptr<ElementOperationParameters> &params )
        : MechanicsElement( params ), d_elementStiffnessMatrix( nullptr )
    {
        d_JxW  = &( d_fe->get_JxW() );
        d_dphi = &( d_fe->get_dphi() );
        d_xyz  = &( d_fe->get_xyz() );
        AMP_INSIST( ( d_useJaumannRate == false ),
                    "Jaumann rate with small strain does not make any sense." );
    }

    //! Destructor.
    virtual ~MechanicsLinearElement() {}

    /**
      This function is used by MechanicsLinearFEOperator to pass the address
      of the element stiffness matrix to this class.
      @param [in] elementStiffnessMatrix Element stiffness matrix
     */
    void setElementStiffnessMatrix( std::vector<std::vector<double>> &elementStiffnessMatrix )
    {
        d_elementStiffnessMatrix = &( elementStiffnessMatrix );
    }

    /**
      Element stiffness matrix computation.
      */
    void apply() override
    {
        if ( d_useReducedIntegration ) {
            apply_Reduced();
        } else {
            apply_Normal();
        }
    }

    /**
      Writes the stess and strain values at the Gauss points within the current element to the file.
      The 6 components of stress and strain at each Gauss point are arranged in the order:
      xx, yy, zz, yz, xz and  xy.
      @param [in] fp File pointer
      @param [in] dispVec Displacements at the nodes of the current element.
     */
    void printStressAndStrain( FILE *fp, const std::vector<double> &dispVec );

    /**
      Computes the stress and strain values at the Gauss points within the current element
      The 6 components of stress and strain at each Gauss point are arranged in the order:
      xx, yy, zz, yz, xz and  xy.
      @param [in] dispVec Displacements at the nodes of the current element.
      @param [out] stressVec Stresses at the Gauss points of the current element.
      @param [out] strainVec Strains at the Gauss points of the current element.
     */
    void computeStressAndStrain( const std::vector<double> &dispVec,
                                 std::vector<double> &stressVec,
                                 std::vector<double> &strainVec );

protected:
    /**
      Element stiffness matrix computation using normal integration scheme.
     */
    void apply_Normal();

    /**
      Element stiffness matrix computation using reduced integration scheme.
     */
    void apply_Reduced();

    const std::vector<libMesh::Real> *d_JxW; /**< Product of the determinant of Jacobian and the quadrature
                                    weight at the Gauss points in the current element. */

    const std::vector<std::vector<libMesh::RealGradient>>
        *d_dphi; /**< Spatial Derivatives of the shape functions at
                  the Gauss points in the current element. */

    const std::vector<libMesh::Point> *d_xyz; /**< Locations of the Gauss points in the current element. */

    std::vector<std::vector<double>> *d_elementStiffnessMatrix; /**< Element stiffness matrix. */

private:
};
} // namespace Operator
} // namespace AMP

#endif
