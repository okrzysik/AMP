
#ifndef included_AMP_NeumannVectorCorrection
#define included_AMP_NeumannVectorCorrection

#include "BoundaryOperator.h"
#include "NeumannVectorCorrectionParameters.h"
#include "RobinPhysicsModel.h"

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

namespace AMP {
namespace Operator {

  /**
    A class to impose Neumann (Flux) Boundary Conditions for both Linear and Nonlinear operator.
    For both the Linear/Nonlinear operator to impose these conditions involves adding the corrections
    to the RHS vector at the appropriate locations. When you do not impose these Neumann condition for
    the weak formulation, a natural condition is assumed.
    This class is also a base class for the Robin Boundary Operator.
    */
  class NeumannVectorCorrection : public BoundaryOperator
  {
    public :

      //! Constructor. This function reads all the parameters required for surface elements.
      NeumannVectorCorrection(const boost::shared_ptr<NeumannVectorCorrectionParameters> & params);

      /**
        Set the variable for the vector that will used with this operator.
        */
      void setVariable(const AMP::LinearAlgebra::Variable::shared_ptr & var) {
        d_variable = var;
      }

      /**
        Destructor
        */
      virtual ~NeumannVectorCorrection() { }

      virtual void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
              AMP::LinearAlgebra::Vector::shared_ptr  &r, const double a = -1.0, const double b = 1.0);

      /**
        This function computes the surface integral for either constant or varrying flux values
        across the boundary.
        */
      void computeRHScorrection(AMP::LinearAlgebra::Vector::shared_ptr rhsCorrection);

      /**
        This function reads parameters related to boundary Ids
        */
      virtual void reset(const boost::shared_ptr<OperatorParameters>& params);

      /**
        Adds a vector to the RHS vector.
        */
      void addRHScorrection(AMP::LinearAlgebra::Vector::shared_ptr rhs); 

      /**
        This function returns a parameter object that can be used to reset the corresponding
        NeumannVectorCorrection operator.
        */
      boost::shared_ptr<OperatorParameters> getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& );

      boost::shared_ptr<OperatorParameters> getParameters() {
        return d_params;
      }

      void setVariableFlux(const AMP::LinearAlgebra::Vector::shared_ptr &flux){
        d_variableFlux = flux->subsetVectorForVariable ( d_variable );
      }

      void setFrozenVector ( AMP::LinearAlgebra::Vector::shared_ptr f );

      boost::shared_ptr<RobinPhysicsModel>  getRobinPhysicsModel () { return d_robinPhysicsModel; }

      std::vector<short int> getBoundaryIds() const { return d_boundaryIds; }

      std::vector<std::vector<unsigned int> > getDofIds() const { return d_dofIds; }

    protected :

      std::vector<short int> d_boundaryIds;

      std::vector<std::vector<double> > d_neumannValues;

      std::vector<std::vector<unsigned int> > d_dofIds;

      AMP::LinearAlgebra::Vector::shared_ptr d_rhsCorrectionAdd;

      //This must be a simple variable not a dual or multivariable
      AMP::LinearAlgebra::Variable::shared_ptr d_variable;

      bool d_isConstantFlux;

      AMP::LinearAlgebra::Vector::shared_ptr d_Frozen;

      AMP::LinearAlgebra::Vector::shared_ptr d_variableFlux;

      std::vector<bool> d_IsCoupledBoundary;

      int d_numIds;

      std::vector<short int> d_numDofIds;

      const std::vector<Real> *d_JxW;

      const std::vector<std::vector<Real> > *d_phi;

      const std::vector<std::vector<RealGradient> > *d_dphi;

      const std::vector<Point> *d_normal;

      boost::shared_ptr < ::FEType > d_feType;

      boost::shared_ptr < ::FEBase > d_fe;

      boost::shared_ptr < ::QBase > d_qrule;

      boost::shared_ptr<NeumannVectorCorrectionParameters> d_params;

      boost::shared_ptr<RobinPhysicsModel> d_robinPhysicsModel;

      std::vector<double> d_beta;

      std::vector<double> d_gamma;

      void createCurrentLibMeshElement();

      void destroyCurrentLibMeshElement();

      void getDofIndicesForCurrentElement();

      std::vector<AMP::Mesh::MeshElement> d_currNodes;

      ::Elem* d_currElemPtr;
      
      std::vector<size_t> d_dofIndices; 

      AMP::Discretization::DOFManager::shared_ptr d_dofManager; 

    private :

  };

}
}

#endif

