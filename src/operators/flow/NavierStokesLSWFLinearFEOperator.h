
#ifndef included_AMP_NavierStokesLSWFLinearFEOperator
#define included_AMP_NavierStokesLSWFLinearFEOperator

/* AMP files */
#include "operators/LinearFEOperator.h"
#include "operators/flow/NavierStokesConstants.h"
#include "operators/flow/NavierStokesLinearFEOperatorParameters.h"
#include "operators/flow/NavierStokesLSWFLinearElement.h"

#include <vector>

namespace AMP {
  namespace Operator {

    class NavierStokesLSWFLinearFEOperator : public LinearFEOperator 
    {
      public :

        NavierStokesLSWFLinearFEOperator(const boost::shared_ptr<NavierStokesLinearFEOperatorParameters>& params);

        ~NavierStokesLSWFLinearFEOperator() { }

        void preAssembly(const boost::shared_ptr<OperatorParameters>& params);

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

//        boost::shared_ptr<AMP::Discretization::DOFManager> d_dofMap[NavierStokes::TOTAL_NUMBER_OF_VARIABLES];
        boost::shared_ptr<AMP::Discretization::DOFManager> d_dofMap;

        std::vector<std::vector<double> > d_elementStiffnessMatrix; 

        boost::shared_ptr< NavierStokesLSWFLinearElement > d_flowLSWFLinElem; 

        boost::shared_ptr<FlowTransportModel> d_transportModel; 

        AMP::LinearAlgebra::Vector::shared_ptr d_inVec;

        AMP::LinearAlgebra::Vector::shared_ptr mySubsetVector(AMP::LinearAlgebra::Vector::shared_ptr vec, 
            AMP::LinearAlgebra::Variable::shared_ptr var);
      private :

        boost::shared_ptr<AMP::LinearAlgebra::Variable> d_inpVariables; /**< Input variables. */

        boost::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariables; /**< Output variables. */

    };

  }
}

#endif
