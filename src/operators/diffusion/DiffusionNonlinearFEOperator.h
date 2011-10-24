#ifndef included_AMP_DiffusionNonlinearFEOperator
#define included_AMP_DiffusionNonlinearFEOperator

#include "ampmesh/MeshVariable.h"
#include "vectors/MultiVariable.h"
#include "vectors/Vector.h"
#include "operators/NonlinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperatorParameters.h"
#include "operators/diffusion/DiffusionNonlinearElement.h"
#include "operators/diffusion/DiffusionConstants.h"

#include <vector>

namespace AMP {
namespace Operator {

class DiffusionNonlinearFEOperator: public NonlinearFEOperator {
public:

    typedef boost::shared_ptr<DiffusionNonlinearFEOperator> shared_ptr;

    DiffusionNonlinearFEOperator(const boost::shared_ptr<DiffusionNonlinearFEOperatorParameters>& params);

    ~DiffusionNonlinearFEOperator() {
    }

    void preAssembly(const boost::shared_ptr<AMP::LinearAlgebra::Vector> &u, boost::shared_ptr<AMP::LinearAlgebra::Vector> &r);

    void postAssembly();

    void preElementOperation(const AMP::Mesh::MeshManager::Adapter::Element &,
            const std::vector<AMP::Mesh::DOFMap::shared_ptr> &);

    void postElementOperation();

    void reset(const boost::shared_ptr<OperatorParameters>&);

    boost::shared_ptr<OperatorParameters>
    getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>&);

    void setInputVariableName(const std::string & name, int varId = -1) {
        if (varId == -1) {
            d_inpVariables->setName(name);
            d_MeshAdapter->appendMeshNameToVariable ( d_inpVariables );
        } else {
            (d_inpVariables->getVariable(varId))->setName(name);
            d_MeshAdapter->appendMeshNameToVariable ( d_inpVariables->getVariable(varId) );
        }
    }

    void setOutputVariableName(const std::string & name, int varId = -1) {
      (void) varId;
        d_outVariable->setName(name);
        d_MeshAdapter->appendMeshNameToVariable ( d_outVariable );
    }

    AMP::LinearAlgebra::Variable::shared_ptr createInputVariable(const std::string & name,
            int varId = -1) {
        if (varId == -1) {
            return d_inpVariables->cloneVariable(name);
        } else {
            return (d_inpVariables->getVariable(varId))->cloneVariable(name);
        }
    }

    AMP::LinearAlgebra::Variable::shared_ptr createOutputVariable(const std::string & name,
            int varId = -1) { (void) varId;
        return d_outVariable->cloneVariable(name);
    }

    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable(int varId = -1) {
            if(varId == -1) {
                    return d_inpVariables;
            } else {
                    return d_inpVariables->getVariable(varId);
            }
    }

    AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
            return d_outVariable;
    }

    unsigned int numberOfDOFMaps() {
            return 1;
    }

    AMP::LinearAlgebra::Variable::shared_ptr getVariableForDOFMap(unsigned int id) {
      (void) id;
            return d_inpVariables->getVariable(d_PrincipalVariable);
    }

    unsigned int getPrincipalVariableId(){return d_PrincipalVariable;}

    std::vector<unsigned int> getNonPrincipalVariableIds(){
        std::vector<unsigned int> ids;
        for (size_t i=0; i<Diffusion::NUMBER_VARIABLES; i++) {
            if (i != d_PrincipalVariable and d_isActive[i]) ids.push_back(i);
        }
        return ids;
    }

    boost::shared_ptr<DiffusionTransportModel> getTransportModel(){
        return d_transportModel;
    }

    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> getFrozen(){return d_Frozen;}

    /**
       This function is used to set frozen vectors in this operator. This is used when some of the 
       variables are solved for in an uncoupled manner.
       @param [in] id Variable Identifier - One of AMP::Diffusion::TEMPERATURE/BURNUP/OXYGEN_CONCENTRATION
       @param [in] frozenVec Frozen vector
       @see DiffusionConstants.h
        */
    void setVector(unsigned int id, AMP::LinearAlgebra::Vector::shared_ptr &frozenVec)
    {
      d_Frozen[id] = frozenVec->subsetVectorForVariable(d_inpVariables->getVariable(id));
      (d_Frozen[id])->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    }
    
    /**
     * checks input to apply operator for satisfaction of range conditions
     */
    bool isValidInput(AMP::LinearAlgebra::Vector::shared_ptr &u);

protected:

    void init(const boost::shared_ptr<DiffusionNonlinearFEOperatorParameters>& params);

    std::vector<unsigned int> d_DofIndices;

    unsigned int d_numNodesForCurrentElement;

    std::vector<double> d_elementOutputVector;

    boost::shared_ptr<DiffusionNonlinearElement> d_diffNonlinElem;

    boost::shared_ptr<DiffusionTransportModel> d_transportModel;

    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> d_inVec;

    AMP::LinearAlgebra::Vector::shared_ptr d_outVec;

    boost::shared_ptr<std::vector<double> > d_TransportGauss;
    AMP::LinearAlgebra::Vector::shared_ptr d_TransportNodal;

    std::vector<bool> d_isActive;

    std::vector<bool> d_isFrozen;

private:

    boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> d_inpVariables;

    boost::shared_ptr<AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 1> > d_outVariable;

    unsigned int d_PrincipalVariable;

    unsigned int d_numberActive;

    unsigned int d_numberFrozen;

    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> d_Frozen;

    void resetFrozen(const boost::shared_ptr<DiffusionNonlinearFEOperatorParameters> &params);
};

}
}

#endif
