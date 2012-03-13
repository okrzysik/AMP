
#ifndef included_AMP_PressureBoundaryVectorCorrection
#define included_AMP_PressureBoundaryVectorCorrection

#include "BoundaryOperator.h"
#include "PressureBoundaryVectorCorrectionParameters.h"

/* Libmesh files */
#include "fe_type.h"
#include "fe_base.h"
#include "elem.h"
#include "quadrature.h"

#include "enum_order.h"
#include "enum_fe_family.h"
#include "enum_quadrature_type.h"
#include "auto_ptr.h"
#include "string_to_enum.h"

#include <string>
#include <vector>

namespace AMP {
namespace Operator {

  /**
    A class to impose Pressure Boundary Conditions for both Linear and Nonlinear mechanics operator.
    For both the Linear/Nonlinear operator to impose these conditions involves adding the corrections
    to the RHS vector at the appropriate locations. When you do not impose these Pressure condition for 
    the weak formulation, a natural condition is assumed.
    */
  class PressureBoundaryVectorCorrection : public BoundaryOperator
  {
    public :

      /**
        Constructor. This function reads all the parameters required for surface elements.
        */

    PressureBoundaryVectorCorrection(const boost::shared_ptr<PressureBoundaryVectorCorrectionParameters> & params);
      /**
        Set the variable for the vector that will be used with this operator.
        */
      void setVariable(const AMP::LinearAlgebra::Variable::shared_ptr & var) {
        d_variable = var;
      }

      /**
        Destructor
        */
      virtual ~PressureBoundaryVectorCorrection() { }

      /**
        Sets pressure values into the appropriate locations of the output vector (r). 
        It also scales the total pressure by an amount "a". It is used for multiple 
        loading steps.
        */
      virtual void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
              AMP::LinearAlgebra::Vector::shared_ptr  &r, const double a = -1.0, const double b = 1.0);

      /**
        This function reads parameters related to boundary Ids
        */
      virtual void reset(const boost::shared_ptr<OperatorParameters>& params);

      /**
        This function computes the surface integral for either constant pressure values 
        across the boundary.
        */
      void addRHScorrection(AMP::LinearAlgebra::Vector::shared_ptr rhs) ;

      /**
        This function returns a parameter object that can be used to reset the corresponding
        PressureBoundaryVectorCorrection operator.
        */
      boost::shared_ptr<OperatorParameters> getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& );

      boost::shared_ptr<OperatorParameters> getParameters() {
        return d_params;
      }

      void setVariablePressure(const AMP::LinearAlgebra::Vector::shared_ptr &pressure){
        d_variablePressure = pressure->subsetVectorForVariable(d_variable);
      }

//      std::vector<std::vector<std::vector<Point> > > getNormals();

    protected :

      std::vector<short int> d_boundaryIds;

      std::vector<double> d_pressureValues;

      AMP::LinearAlgebra::Vector::shared_ptr d_rhsCorrectionAdd;

      //This must be a simple variable not a dual or multivariable
      AMP::LinearAlgebra::Variable::shared_ptr d_variable;

      bool d_isConstantPressure;

      AMP::LinearAlgebra::Vector::shared_ptr d_Frozen;

      AMP::LinearAlgebra::Vector::shared_ptr d_variablePressure;

      int d_numIds;

      const std::vector<Real> *d_JxW;

      const std::vector<std::vector<Real> > *d_phi;

      const std::vector<std::vector<RealGradient> > *d_dphi;

      const std::vector<Point> *d_normal;

      boost::shared_ptr < ::FEType > d_feType;

      boost::shared_ptr < ::FEBase > d_fe;

      boost::shared_ptr < ::FEBase > d_fe_3d;

      boost::shared_ptr < ::QBase > d_qrule;

      const ::Elem *d_elem;

      const ::Elem *e_elem;

      boost::shared_ptr<PressureBoundaryVectorCorrectionParameters> d_params;

      void createCurrentLibMeshElement();

      void destroyCurrentLibMeshElement();

      void createCurrentLibMeshSide();

      void destroyCurrentLibMeshSide();

      void getDofIndicesForCurrentSide();

      std::vector<AMP::Mesh::MeshElement> d_currNodes;
      std::vector<AMP::Mesh::MeshElement> d_currFaces;

      ::Elem* d_currElemPtr;
      ::Elem* d_currSidePtr;
      
      std::vector<size_t> d_dofIndices; 

      AMP::Discretization::DOFManager::shared_ptr d_dofManager; 


    private :

  };

}
}

#endif

