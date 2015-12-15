
#ifndef included_AMP_NavierStokesLSWFLinearFEOperator
#define included_AMP_NavierStokesLSWFLinearFEOperator

/* AMP files */
#include "operators/libmesh/LinearFEOperator.h"
#include "operators/flow/NavierStokesConstants.h"
#include "operators/flow/NavierStokesLinearFEOperatorParameters.h"
#include "operators/flow/NavierStokesLSWFLinearElement.h"

#include <vector>

namespace AMP {
  namespace Operator {

    class NavierStokesLSWFLinearFEOperator : public LinearFEOperator 
    {
      public :

        explicit NavierStokesLSWFLinearFEOperator(const AMP::shared_ptr<NavierStokesLinearFEOperatorParameters>& params);

        virtual ~NavierStokesLSWFLinearFEOperator() { }

        void preAssembly(const AMP::shared_ptr<OperatorParameters>& params);

        void postAssembly();

        void preElementOperation(const AMP::Mesh::MeshElement &);

        void postElementOperation();

        AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() {
          return d_inpVariables;
        }

        AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
          return d_outVariables;
        }

      protected :

        void getDofIndicesForCurrentElement(int varId, std::vector<std::vector<size_t> > & dofIds);

        std::vector<std::vector<size_t> > d_type0DofIndices; 
        std::vector<std::vector<size_t> > d_type1DofIndices; 

//        AMP::shared_ptr<AMP::Discretization::DOFManager> d_dofMap[NavierStokes::TOTAL_NUMBER_OF_VARIABLES];
        AMP::shared_ptr<AMP::Discretization::DOFManager> d_dofMap;

        std::vector<std::vector<double> > d_elementStiffnessMatrix; 

        AMP::shared_ptr< NavierStokesLSWFLinearElement > d_flowLSWFLinElem; 

        AMP::shared_ptr<FlowTransportModel> d_transportModel; 

        AMP::LinearAlgebra::Vector::shared_ptr d_inVec;

        AMP::LinearAlgebra::Vector::shared_ptr mySubsetVector(AMP::LinearAlgebra::Vector::shared_ptr vec, 
            AMP::LinearAlgebra::Variable::shared_ptr var);
      private :

        AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_inpVariables; /**< Input variables. */

        AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariables; /**< Output variables. */

    };

  }
}

#endif
