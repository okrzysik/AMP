#ifndef included_AMP_DiffusionNonlinearFEOperator
#define included_AMP_DiffusionNonlinearFEOperator

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

    ~DiffusionNonlinearFEOperator() {}

    void preAssembly(const boost::shared_ptr<AMP::LinearAlgebra::Vector> &u, boost::shared_ptr<AMP::LinearAlgebra::Vector> &r);

    void postAssembly();

    void preElementOperation(const AMP::Mesh::MeshElement &, const std::vector<AMP::Discretization::DOFManager::shared_ptr> &);

    void postElementOperation();

    void reset(const boost::shared_ptr<OperatorParameters>&);

    boost::shared_ptr<OperatorParameters>
    getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>&);

    void setInputVariableName(const std::string & name, int varId = -1);

    void setOutputVariableName(const std::string & name, int varId = -1);

    AMP::LinearAlgebra::Variable::shared_ptr createInputVariable(const std::string & name, int varId = -1);

    AMP::LinearAlgebra::Variable::shared_ptr createOutputVariable(const std::string & name, int varId = -1);

    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable(int varId = -1);

    AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable();

    unsigned int numberOfDOFMaps();

    AMP::LinearAlgebra::Variable::shared_ptr getVariableForDOFMap(unsigned int id);

    unsigned int getPrincipalVariableId();

    std::vector<unsigned int> getNonPrincipalVariableIds();

    boost::shared_ptr<DiffusionTransportModel> getTransportModel();

    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> getFrozen();

    /**
       This function is used to set frozen vectors in this operator. This is used when some of the 
       variables are solved for in an uncoupled manner.
       @param [in] id Variable Identifier - One of AMP::Diffusion::TEMPERATURE/BURNUP/OXYGEN_CONCENTRATION
       @param [in] frozenVec Frozen vector
       @see DiffusionConstants.h
        */
    void setVector(unsigned int id, AMP::LinearAlgebra::Vector::shared_ptr &frozenVec);
    
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

    //boost::shared_ptr<AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 1> > d_outVariable;

    unsigned int d_PrincipalVariable;

    unsigned int d_numberActive;

    unsigned int d_numberFrozen;

    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> d_Frozen;

    void resetFrozen(const boost::shared_ptr<DiffusionNonlinearFEOperatorParameters> &params);
};

}
}

#endif
