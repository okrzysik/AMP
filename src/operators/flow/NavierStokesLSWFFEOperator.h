
#ifndef included_AMP_NavierStokesLSWFFEOperator
#define included_AMP_NavierStokesLSWFFEOperator

/* AMP files */
#include "AMP/ampmesh/MeshElement.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/operators/flow/NavierStokesConstants.h"
#include "AMP/operators/flow/NavierStokesLSWFElement.h"
#include "AMP/operators/flow/NavierStokesLSWFFEOperatorParameters.h"
#include "AMP/operators/libmesh/NonlinearFEOperator.h"
#include "AMP/vectors/MultiVariable.h"
#include "AMP/vectors/Variable.h"

#include <vector>

namespace AMP {
namespace Operator {

class NavierStokesLSWFFEOperator : public NonlinearFEOperator
{
public:
    explicit NavierStokesLSWFFEOperator(
        const std::shared_ptr<NavierStokesLSWFFEOperatorParameters> &params );

    virtual ~NavierStokesLSWFFEOperator() {}

    void preAssembly( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                      std::shared_ptr<AMP::LinearAlgebra::Vector> r ) override;

    void postAssembly() override;

    void preElementOperation( const AMP::Mesh::MeshElement & ) override;

    void postElementOperation() override;

    void reset( const std::shared_ptr<OperatorParameters> & ) override;

    /*
            void setVector(unsigned int id, AMP::LinearAlgebra::Vector::shared_ptr frozenVec) {
              AMP::LinearAlgebra::Variable::shared_ptr var = d_inpVariables->getVariable(id);
              d_inVec[id] = mySubsetVector(frozenVec, var);
              (d_inVec[id])->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET
       );
            }
    */
    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() override { return d_inpVariables; }

    AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() override { return d_outVariables; }

protected:
    std::shared_ptr<OperatorParameters>
    getJacobianParameters( AMP::LinearAlgebra::Vector::const_shared_ptr u ) override;

    AMP::LinearAlgebra::Vector::shared_ptr
    mySubsetVector( AMP::LinearAlgebra::Vector::shared_ptr vec,
                    AMP::LinearAlgebra::Variable::shared_ptr var );

    AMP::LinearAlgebra::Vector::const_shared_ptr
    mySubsetVector( AMP::LinearAlgebra::Vector::const_shared_ptr vec,
                    AMP::LinearAlgebra::Variable::shared_ptr var );

    void getDofIndicesForCurrentElement( int varId, std::vector<std::vector<size_t>> &dofIds );

    std::vector<double> d_elementOutputVector;

    std::shared_ptr<NavierStokesLSWFElement> d_nsLSWFElem;

    std::shared_ptr<FlowTransportModel> d_transportModel;

    //        std::vector<AMP::LinearAlgebra::Vector::shared_ptr> d_inVec;
    AMP::LinearAlgebra::Vector::const_shared_ptr d_inVec;

    AMP::LinearAlgebra::Vector::shared_ptr d_outVec;

    std::vector<bool> d_isActive;

    std::vector<bool> d_isFrozen;

    std::shared_ptr<AMP::LinearAlgebra::Variable> d_inpVariables;
    std::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariables;

    //        std::shared_ptr<AMP::Discretization::DOFManager>
    //        d_dofMap[NavierStokes::TOTAL_NUMBER_OF_VARIABLES];
    std::shared_ptr<AMP::Discretization::DOFManager> d_dofMap;

    std::vector<AMP::Mesh::MeshElement> d_currNodes;

    std::vector<std::vector<size_t>> d_type0DofIndices;

    std::vector<std::vector<size_t>> d_type1DofIndices;
};
} // namespace Operator
} // namespace AMP

#endif
