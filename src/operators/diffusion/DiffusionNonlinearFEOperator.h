
#ifndef included_AMP_DiffusionNonlinearFEOperator
#define included_AMP_DiffusionNonlinearFEOperator

#include "operators/diffusion/DiffusionConstants.h"
#include "operators/diffusion/DiffusionNonlinearElement.h"
#include "operators/diffusion/DiffusionNonlinearFEOperatorParameters.h"
#include "operators/libmesh/NonlinearFEOperator.h"
#include "vectors/MultiVariable.h"
#include "vectors/Vector.h"

#include <vector>

namespace AMP {
namespace Operator {

class DiffusionNonlinearFEOperator : public NonlinearFEOperator
{
public:
    typedef AMP::shared_ptr<DiffusionNonlinearFEOperator> shared_ptr;

    explicit DiffusionNonlinearFEOperator(
        const AMP::shared_ptr<DiffusionNonlinearFEOperatorParameters> &params );

    virtual ~DiffusionNonlinearFEOperator() {}

    void reset( const AMP::shared_ptr<OperatorParameters> & );

    void setInputVariableName( const std::string &name, int varId = -1 );

    void setOutputVariableName( const std::string &name, int varId = -1 );

    AMP::LinearAlgebra::Variable::shared_ptr createInputVariable( const std::string &name,
                                                                  int varId = -1 );

    AMP::LinearAlgebra::Variable::shared_ptr createOutputVariable( const std::string &name,
                                                                   int varId = -1 );

    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable();

    AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable();

    unsigned int numberOfDOFMaps();

    AMP::LinearAlgebra::Variable::shared_ptr getVariableForDOFMap( unsigned int id );

    unsigned int getPrincipalVariableId();

    std::vector<unsigned int> getNonPrincipalVariableIds();

    AMP::shared_ptr<DiffusionTransportModel> getTransportModel();

    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> getFrozen();

    /**
      This function is used to set frozen vectors in this operator. This is used when some of the
      variables are solved for in an uncoupled manner.
      @param [in] id Variable Identifier - One of
      AMP::Diffusion::TEMPERATURE/BURNUP/OXYGEN_CONCENTRATION
      @param [in] frozenVec Frozen vector
      @see DiffusionConstants.h
      */
    void setVector( unsigned int id, AMP::LinearAlgebra::Vector::shared_ptr frozenVec );

    /**
     * checks input to apply operator for satisfaction of range conditions
     */
    bool isValidInput( AMP::LinearAlgebra::Vector::shared_ptr &u );

protected:
    AMP::shared_ptr<OperatorParameters>
    getJacobianParameters( AMP::LinearAlgebra::Vector::const_shared_ptr u ) override;

    void preAssembly( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                      AMP::LinearAlgebra::Vector::shared_ptr r );

    void postAssembly();

    void preElementOperation( const AMP::Mesh::MeshElement & );

    void postElementOperation();

    void init( const AMP::shared_ptr<DiffusionNonlinearFEOperatorParameters> &params );

    std::vector<double> d_elementOutputVector;

    AMP::shared_ptr<DiffusionNonlinearElement> d_diffNonlinElem;

    AMP::shared_ptr<DiffusionTransportModel> d_transportModel;

    std::vector<AMP::LinearAlgebra::Vector::const_shared_ptr> d_inVec;

    std::vector<AMP::Mesh::MeshElement> d_currNodes;

    AMP::LinearAlgebra::Vector::shared_ptr d_outVec;

    AMP::shared_ptr<std::vector<double>> d_TransportGauss;
    AMP::LinearAlgebra::Vector::shared_ptr d_TransportNodal;

    std::vector<bool> d_isActive;

    std::vector<bool> d_isFrozen;

private:
    AMP::shared_ptr<AMP::LinearAlgebra::MultiVariable> d_inpVariables;

    AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariable;

    unsigned int d_PrincipalVariable;

    unsigned int d_numberActive;

    unsigned int d_numberFrozen;

    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> d_Frozen;

    void resetFrozen( const AMP::shared_ptr<DiffusionNonlinearFEOperatorParameters> params );
};
}
}

#endif
