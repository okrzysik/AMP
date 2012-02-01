
#ifndef included_AMP_ConsMomentumGalWFFEOperator
#define included_AMP_ConsMomentumGalWFFEOperator

/* AMP files */
#include "NonlinearFEOperator.h"
#include "NavierStokesConstants.h"
#include "ConsMomentumGalWFFEOperatorParameters.h"
#include "ConsMomentumGalWFElement.h"
#include "vectors/MultiVariable.h"

#include <vector>


namespace AMP {
namespace Operator {


class ConsMomentumGalWFFEOperator : public NonlinearFEOperator 
{
public :

        ConsMomentumGalWFFEOperator(const boost::shared_ptr<ConsMomentumGalWFFEOperatorParameters>& params);

        ~ConsMomentumGalWFFEOperator() { }

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

        boost::shared_ptr<ConsMomentumGalWFElement> d_flowGalWFElem; 

        boost::shared_ptr<FlowTransportModel> d_transportModel; /**< Flow Transport model. */

        std::vector<AMP::LinearAlgebra::Vector::shared_ptr> d_inVec; /**< Input vector. */

        AMP::LinearAlgebra::Vector::shared_ptr d_referenceTemperature; 

        AMP::LinearAlgebra::Vector::shared_ptr d_outVec; 

        std::vector<bool> d_isActive; 

        std::vector<bool> d_isFrozen;

private :

        bool d_isInitialized; /**< A flag that is true if init() has been called and false otherwsie. */

        boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> d_inpVariables; /**< Input variables. */

        boost::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariables; /**< Output variable. */

};


}
}


#endif

