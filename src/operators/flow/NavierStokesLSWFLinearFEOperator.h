
#ifndef included_AMP_NavierStokesLSWFLinearFEOperator
#define included_AMP_NavierStokesLSWFLinearFEOperator

/* AMP files */
#include "AMP/operators/flow/NavierStokesConstants.h"
#include "AMP/operators/flow/NavierStokesLSWFLinearElement.h"
#include "AMP/operators/flow/NavierStokesLinearFEOperatorParameters.h"
#include "AMP/operators/libmesh/LinearFEOperator.h"

#include <vector>

namespace AMP::Operator {

class NavierStokesLSWFLinearFEOperator : public LinearFEOperator
{
public:
    explicit NavierStokesLSWFLinearFEOperator(
        std::shared_ptr<const NavierStokesLinearFEOperatorParameters> params );

    virtual ~NavierStokesLSWFLinearFEOperator() {}

    std::string type() const override { return "NavierStokesLSWFLinearFEOperator"; }

    void preAssembly( std::shared_ptr<const OperatorParameters> params ) override;

    void postAssembly() override;

    void preElementOperation( const AMP::Mesh::MeshElement & ) override;

    void postElementOperation() override;

protected:
    void getDofIndicesForCurrentElement( int varId, std::vector<std::vector<size_t>> &dofIds );

    std::vector<std::vector<size_t>> d_type0DofIndices;
    std::vector<std::vector<size_t>> d_type1DofIndices;

    //        std::shared_ptr<AMP::Discretization::DOFManager>
    //        d_dofMap[NavierStokes::TOTAL_NUMBER_OF_VARIABLES];
    std::shared_ptr<AMP::Discretization::DOFManager> d_dofMap;

    std::vector<std::vector<double>> d_elementStiffnessMatrix;

    std::shared_ptr<NavierStokesLSWFLinearElement> d_flowLSWFLinElem;

    std::shared_ptr<FlowTransportModel> d_transportModel;

    AMP::LinearAlgebra::Vector::shared_ptr d_inVec;

    AMP::LinearAlgebra::Vector::shared_ptr
    mySubsetVector( AMP::LinearAlgebra::Vector::shared_ptr vec,
                    std::shared_ptr<AMP::LinearAlgebra::Variable> var );

private:
};
} // namespace AMP::Operator

#endif
