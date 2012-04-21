
#ifndef included_AMP_NavierStokesGalWFFEOperator
#define included_AMP_NavierStokesGalWFFEOperator

/* AMP files */
#include "vectors/MultiVariable.h"
#include "vectors/Vector.h"
#include "operators/NonlinearFEOperator.h"

#include "operators/NonlinearFEOperator.h"
#include "operators/flow/NavierStokesConstants.h"
#include "operators/flow/NavierStokesGalWFFEOperatorParameters.h"
#include "operators/flow/NavierStokesGalWFElement.h"

#include <vector>

namespace AMP {
  namespace Operator {

    class NavierStokesGalWFFEOperator : public NonlinearFEOperator 
    {
      public :

        NavierStokesGalWFFEOperator(const boost::shared_ptr<NavierStokesGalWFFEOperatorParameters>& params);

        ~NavierStokesGalWFFEOperator() { }

        void preAssembly(const boost::shared_ptr<AMP::LinearAlgebra::Vector>  &u, boost::shared_ptr<AMP::LinearAlgebra::Vector>  &r);

        void postAssembly();

        void preElementOperation(const AMP::Mesh::MeshElement &);

        void postElementOperation();

        void reset(const boost::shared_ptr<OperatorParameters>& );

        boost::shared_ptr<OperatorParameters> 
          getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& );

        void init();

        void setVector(unsigned int id, AMP::LinearAlgebra::Vector::shared_ptr &frozenVec) {
          d_inVec[id] = frozenVec->subsetVectorForVariable(d_inpVariables->getVariable(id));
          (d_inVec[id])->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
        }

        static AMP::LinearAlgebra::Variable::shared_ptr createInputVariable(const std::string & name, int varId = -1);

        static AMP::LinearAlgebra::Variable::shared_ptr createOutputVariable(const std::string & name, int varId = -1) {
          (void) varId;      
          AMP::LinearAlgebra::Variable::shared_ptr outVar(new AMP::LinearAlgebra::Variable(name) );
          return outVar;
        }

        AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() {
          return d_inpVariables; 
        }

        AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
          return d_outVariables;
        }

        unsigned int numberOfDOFMaps();

        AMP::LinearAlgebra::Variable::shared_ptr getVariableForDOFMap(unsigned int id);

      protected :

        std::vector<unsigned int> d_type0DofIndices[3]; 
        std::vector<unsigned int> d_type1DofIndices; 

        unsigned int d_numNodesForCurrentElement; 

        std::vector<double> d_elementOutputVector; 

        boost::shared_ptr<NavierStokesGalWFElement> d_flowGalWFElem; 

        boost::shared_ptr<FlowTransportModel> d_transportModel; /**< Flow Transport model. */

        std::vector<AMP::LinearAlgebra::Vector::shared_ptr> d_inVec; /**< Input vector. */

        AMP::LinearAlgebra::Vector::shared_ptr d_referenceTemperature; 

        AMP::LinearAlgebra::Vector::shared_ptr d_outVec; 

        std::vector<bool> d_isActive; 

        std::vector<bool> d_isFrozen;

        bool d_coupledFormulation;

      private :

        bool d_isInitialized; /**< A flag that is true if init() has been called and false otherwsie. */

        boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> d_inpVariables; /**< Input variables. */

        boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> d_outVariables; /**< Output variable. */

    };

  }
}

#endif



