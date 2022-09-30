
#ifndef included_AMP_NavierStokesGalWFLinearFEOperator
#define included_AMP_NavierStokesGalWFLinearFEOperator

// AMP files
#include "AMP/operators/flow/NavierStokesConstants.h"
#include "AMP/operators/flow/NavierStokesGalWFLinearElement.h"
#include "AMP/operators/flow/NavierStokesLinearFEOperatorParameters.h"
#include "AMP/operators/libmesh/LinearFEOperator.h"
#include "AMP/vectors/MultiVariable.h"

#include <vector>

namespace AMP::Operator {

class NavierStokesGalWFLinearFEOperator : public LinearFEOperator
{
public:
    explicit NavierStokesGalWFLinearFEOperator(
        std::shared_ptr<const NavierStokesLinearFEOperatorParameters> params );

    virtual ~NavierStokesGalWFLinearFEOperator() {}

    std::string type() const override { return "NavierStokesGalWFLinearFEOperator"; }

    void preAssembly( std::shared_ptr<const OperatorParameters> params ) override;

    void postAssembly() override;

    void preElementOperation( const AMP::Mesh::MeshElement & ) override;

    void postElementOperation() override;

    unsigned int numberOfDOFMaps() { return 1; }


protected:
    void gettype0DofIndicesForCurrentElement( int varId, std::vector<std::vector<size_t>> &dofIds );
    void gettype1DofIndicesForCurrentElement( int varId, std::vector<size_t> &dofIds );

    void createHex27LibMeshElement();
    void destroyHex27LibMeshElement();

    std::vector<std::vector<size_t>> d_type0DofIndices; /**< Primary DOF indices */
    std::vector<size_t> d_type1DofIndices;

    std::shared_ptr<AMP::Discretization::DOFManager>
        d_dofMap[NavierStokes::TOTAL_NUMBER_OF_VARIABLES];

    std::vector<std::vector<double>> d_elementStiffnessMatrix;

    std::shared_ptr<NavierStokesGalWFLinearElement> d_flowGalWFLinElem;

    std::shared_ptr<FlowTransportModel> d_transportModel;

    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> d_inVec;

    unsigned int d_numNodesForCurrentElement; /**< Number of nodes in the current element. */

private:
};
} // namespace AMP::Operator

#endif
