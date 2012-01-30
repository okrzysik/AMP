
#ifndef included_AMP_ConsMassGalWFLinearFEOperator
#define included_AMP_ConsMassGalWFLinearFEOperator

/* AMP files */
#include "LinearFEOperator.h"
#include "NavierStokesConstants.h"
#include "ConsMassGalWFLinearFEOperatorParameters.h"
#include "ConsMassGalWFLinearElement.h"
#include "vectors/MultiVariable.h"

#include <vector>


namespace AMP {
namespace Operator {


class ConsMassGalWFLinearFEOperator : public LinearFEOperator 
{
public :

        ConsMassGalWFLinearFEOperator(const boost::shared_ptr<ConsMassGalWFLinearFEOperatorParameters>& params);

        ~ConsMassGalWFLinearFEOperator() { }

        void preAssembly(const boost::shared_ptr<OperatorParameters>& params);

        void postAssembly();

        void preElementOperation(const AMP::Mesh::MeshElement &);

        void postElementOperation();

        AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() {
          return d_inpVariable;
        }

        AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
          return d_outVariable;
        }

        unsigned int numberOfDOFMaps() {
          return 2;
        }

        AMP::LinearAlgebra::Variable::shared_ptr getVariableForDOFMap(unsigned int id) {
          if(id == 0){
            return (d_outVariable);
          }else {
            return (d_inpVariable);
          }
        }

protected :

        std::vector<unsigned int> d_dofIndices0[3]; 
        std::vector<unsigned int> d_dofIndices1; 

        std::vector<std::vector<double> > d_elementStiffnessMatrix; 

        boost::shared_ptr< ConsMassGalWFLinearElement > d_flowGalWFLinElem; 

        boost::shared_ptr<FlowTransportModel> d_transportModel; 

        std::vector<AMP::LinearAlgebra::Vector::shared_ptr> d_inVec;

        unsigned int d_numNodesForCurrentElement; /**< Number of nodes in the current element. */

private :

        boost::shared_ptr<AMP::LinearAlgebra::Variable> d_inpVariable; /**< Input variables. */

        boost::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariable; /**< Output variable. */

};


}
}


#endif


