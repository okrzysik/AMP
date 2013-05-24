
#ifndef included_AMP_MechanicsLinearUpdatedLagrangianElement
#define included_AMP_MechanicsLinearUpdatedLagrangianElement

#include <vector>

#include "boost/shared_ptr.hpp"

/* AMP files */
#include "MechanicsElement.h"
#include "UpdatedLagrangianUtils.h"
#include "MechanicsConstants.h"

namespace AMP {
namespace Operator {


/**
  A class for representing the element level computation performed within a 
  linear finite element operator for modelling solid mechanics.
  The linear operator could either be a linear elasticity operator or it could be
  the jacobian of a nonlinear elasticity/elasto-plasticity operator.
 */
class MechanicsLinearUpdatedLagrangianElement : public MechanicsElement 
{
public :

    //! Constructor.
    MechanicsLinearUpdatedLagrangianElement(const boost::shared_ptr<ElementOperationParameters>& params) :
        MechanicsElement(params),
        d_elementStiffnessMatrix(NULL)
    { 
        d_JxW = &(d_fe->get_JxW());
        d_dphi = &(d_fe->get_dphi());
        d_xyz = &(d_fe->get_xyz());
        d_onePointShearIntegration = false;
    }

    //! Destructor.
    virtual ~MechanicsLinearUpdatedLagrangianElement() {  }

    /**
      This function is used by MechanicsLinearFEOperator to pass the address 
      of the element stiffness matrix to this class. 
      @param [in] elementStiffnessMatrix Element stiffness matrix
      */
    void setElementStiffnessMatrix( std::vector<std::vector<double> > & elementStiffnessMatrix )
    {
        d_elementStiffnessMatrix = &(elementStiffnessMatrix);
    }

    /**
      Element stiffness matrix computation.
     */
    void apply() 
    {
        if(d_useReducedIntegration) {
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
    void printStressAndStrain(FILE* fp, const std::vector<double> & dispVec); 

    /**
      Computes the stress and strain values at the Gauss points within the current element
      The 6 components of stress and strain at each Gauss point are arranged in the order:
      xx, yy, zz, yz, xz and  xy.
      @param [in] dispVec Displacements at the nodes of the current element.
      @param [out] stressVec Stresses at the Gauss points of the current element. 
      @param [out] strainVec Strains at the Gauss points of the current element.
     */
    void computeStressAndStrain(const std::vector<double> & dispVec, 
        std::vector<double> & stressVec, std::vector<double> & strainVec);

    /**
      This function is used by MechanicsLinearFEOperator to pass the address 
      of the element Input vector to this class. 
      @param [in] elementInputVectors Element input vector
     */
    void setElementVectors( const std::vector<std::vector<double> > & elementInputVectors ) 
    {
        d_elementInputVectors = elementInputVectors;
    }

    /**
      Computes the deformation gradient at (n+1)-th time step.
    void computeDeformationGradientLin(const std::vector<std::vector<RealGradient> > & dphi, 
        const std::vector<Point> & xyz, unsigned int num_nodes, unsigned int qp, double F[3][3]);
    */

    /**
      Initializes the reference x, y and z coordinates.
     */
    void initializeReferenceXYZ(std::vector<double> & elementRefXYZ);

    /**
      Assign the reference x, y and z coordinates for
      the current element.
     */
    void assignReferenceXYZ(std::vector<double> elementRefXYZ)
    {
        d_elementRefXYZ = elementRefXYZ;
    }

protected :

    /**
      Element stiffness matrix computation using normal integration scheme.
     */
    void apply_Normal();

    /**
      Element stiffness matrix computation using reduced integration scheme.
     */
    void apply_Reduced();

    const std::vector<Real> *d_JxW; /**< Product of the determinant of Jacobian and the quadrature 
                                    weight at the Gauss points in the current element. */

    const std::vector<std::vector<RealGradient> > *d_dphi; /**< Spatial Derivatives of the shape functions at
                                                           the Gauss points in the current element. */

    const std::vector<Point> *d_xyz; /**< Locations of the Gauss points in the current element. */

    std::vector<std::vector<double> > *d_elementStiffnessMatrix; /**< Element stiffness matrix. */

    std::vector<std::vector<double> > d_elementInputVectors; /**< Element input vectors (Displacement, temperature, burnup etc). */

    bool d_onePointShearIntegration;

    std::vector<double> d_elementRefXYZ;

private :

};


}
}

#endif


