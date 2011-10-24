
#ifndef included_AMP_NavierStokesGalWFFEOperator
#define included_AMP_NavierStokesGalWFFEOperator

/* AMP files */
#include "NonlinearFEOperator.h"
#include "NavierStokesConstants.h"
#include "NavierStokesGalWFFEOperatorParameters.h"
#include "NavierStokesGalWFElement.h"
#include "ampmesh/MeshVariable.h"
#include "vectors/MultiVariable.h"

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

      void preElementOperation(const AMP::Mesh::MeshManager::Adapter::Element &, const std::vector<AMP::Mesh::DOFMap::shared_ptr> &);

      void postElementOperation();

      void reset(const boost::shared_ptr<OperatorParameters>& );

      boost::shared_ptr<OperatorParameters> 
        getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& );

      void init();

      void setVector(unsigned int id, AMP::LinearAlgebra::Vector::shared_ptr &frozenVec) {
        d_inVec[id] = frozenVec->subsetVectorForVariable(d_inpVariables->getVariable(id));
        (d_inVec[id])->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
      }

      void setInputVariableName(const std::string & name, int varId = -1) {
        if(varId == -1) {
          d_inpVariables->setName(name);
        } else {
          (d_inpVariables->getVariable(varId))->setName(name);
        }
      }

      void setOutputVariableName(const std::string & name, int varId = -1) {
        (void) varId;      
        if(varId == -1) {
          d_outVariables->setName(name);
        } else {
          (d_outVariables->getVariable(varId))->setName(name);
        }
      }

      static AMP::LinearAlgebra::Variable::shared_ptr createInputVariable(const std::string & name, int varId = -1);

      static AMP::LinearAlgebra::Variable::shared_ptr createOutputVariable(const std::string & name, int varId = -1) {
        (void) varId;      
        AMP::LinearAlgebra::Variable::shared_ptr outVar(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 4>(name) );
        return outVar;
      }

      AMP::LinearAlgebra::Variable::shared_ptr getInputVariable(int varId = -1) {
        if(varId == -1) {
          return d_inpVariables; 
        } else {
          return d_inpVariables->getVariable(varId);
        }
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

      boost::shared_ptr<AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 4> > d_outVariables; /**< Output variable. */

  };

}
}

#endif



