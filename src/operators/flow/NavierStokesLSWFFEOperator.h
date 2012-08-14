
#ifndef included_AMP_NavierStokesLSWFFEOperator
#define included_AMP_NavierStokesLSWFFEOperator

/* AMP files */
#include "operators/NonlinearFEOperator.h"
#include "operators/flow/NavierStokesConstants.h"
#include "operators/flow/NavierStokesLSWFFEOperatorParameters.h"
#include "operators/flow/NavierStokesLSWFElement.h"
#include "vectors/MultiVariable.h"
#include "vectors/Variable.h"
#include "discretization/DOF_Manager.h"
#include "ampmesh/MeshElement.h"

#include <vector>

namespace AMP {
  namespace Operator {

    class NavierStokesLSWFFEOperator : public NonlinearFEOperator 
    {
      public :

        NavierStokesLSWFFEOperator(const boost::shared_ptr<NavierStokesLSWFFEOperatorParameters>& params);

        ~NavierStokesLSWFFEOperator() { }

        void preAssembly(AMP::LinearAlgebra::Vector::const_shared_ptr u, boost::shared_ptr<AMP::LinearAlgebra::Vector>  &r);

        void postAssembly();

        void preElementOperation(const AMP::Mesh::MeshElement &);

        void postElementOperation();

        void reset(const boost::shared_ptr<OperatorParameters>& );

        boost::shared_ptr<OperatorParameters> 
          getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& );

/*        
        void setVector(unsigned int id, AMP::LinearAlgebra::Vector::shared_ptr &frozenVec) {
          AMP::LinearAlgebra::Variable::shared_ptr var = d_inpVariables->getVariable(id);
          d_inVec[id] = mySubsetVector(frozenVec, var);
          (d_inVec[id])->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
        }
*/
        AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() {
          return d_inpVariables; 
        }

        AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
          return d_outVariables;
        }

      protected :

        AMP::LinearAlgebra::Vector::shared_ptr mySubsetVector(AMP::LinearAlgebra::Vector::shared_ptr vec, 
            AMP::LinearAlgebra::Variable::shared_ptr var);

        AMP::LinearAlgebra::Vector::const_shared_ptr mySubsetVector(AMP::LinearAlgebra::Vector::const_shared_ptr vec, 
            AMP::LinearAlgebra::Variable::shared_ptr var);

        void getDofIndicesForCurrentElement(int varId, std::vector<std::vector<size_t> > & dofIds);

        std::vector<double> d_elementOutputVector; 

        boost::shared_ptr<NavierStokesLSWFElement> d_nsLSWFElem; 

        boost::shared_ptr<FlowTransportModel> d_transportModel; 

//        std::vector<AMP::LinearAlgebra::Vector::shared_ptr> d_inVec; 
        AMP::LinearAlgebra::Vector::const_shared_ptr d_inVec; 

        AMP::LinearAlgebra::Vector::shared_ptr d_outVec; 

        std::vector<bool> d_isActive; 

        std::vector<bool> d_isFrozen; 

        boost::shared_ptr<AMP::LinearAlgebra::Variable> d_inpVariables;
        boost::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariables; 

//        boost::shared_ptr<AMP::Discretization::DOFManager> d_dofMap[NavierStokes::TOTAL_NUMBER_OF_VARIABLES];
        boost::shared_ptr<AMP::Discretization::DOFManager> d_dofMap;

        std::vector<AMP::Mesh::MeshElement> d_currNodes; 

        std::vector<std::vector<size_t> > d_type0DofIndices;

        std::vector<std::vector<size_t> > d_type1DofIndices;
    };

  }
}

#endif



