
#ifndef included_AMP_NavierStokesGalWFLinearFEOperator
#define included_AMP_NavierStokesGalWFLinearFEOperator

// AMP files
#include "operators/flow/NavierStokesConstants.h"
#include "operators/flow/NavierStokesGalWFLinearElement.h"
#include "operators/flow/NavierStokesLinearFEOperatorParameters.h"
#include "operators/libmesh/LinearFEOperator.h"
#include "vectors/MultiVariable.h"

#include <vector>

namespace AMP {
namespace Operator {

class NavierStokesGalWFLinearFEOperator : public LinearFEOperator
{
public:
    explicit NavierStokesGalWFLinearFEOperator(
        const AMP::shared_ptr<NavierStokesLinearFEOperatorParameters> &params );

    virtual ~NavierStokesGalWFLinearFEOperator() {}

    void preAssembly( const AMP::shared_ptr<OperatorParameters> &params );

    void postAssembly();

    void preElementOperation( const AMP::Mesh::MeshElement & );

    void postElementOperation();

    unsigned int numberOfDOFMaps() { return 1; }

    AMP::LinearAlgebra::Variable::shared_ptr getVariableForDOFMap( unsigned int )
    {
        return d_inputVariable;
    }

protected:
    void gettype0DofIndicesForCurrentElement( int varId, std::vector<std::vector<size_t>> &dofIds );
    void gettype1DofIndicesForCurrentElement( int varId, std::vector<size_t> &dofIds );

    void createHex27LibMeshElement();
    void destroyHex27LibMeshElement();

    std::vector<std::vector<size_t>> d_type0DofIndices; /**< Primary DOF indices */
    std::vector<size_t> d_type1DofIndices;

    AMP::shared_ptr<AMP::Discretization::DOFManager>
        d_dofMap[NavierStokes::TOTAL_NUMBER_OF_VARIABLES];

    std::vector<std::vector<double>> d_elementStiffnessMatrix;

    AMP::shared_ptr<NavierStokesGalWFLinearElement> d_flowGalWFLinElem;

    AMP::shared_ptr<FlowTransportModel> d_transportModel;

    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> d_inVec;

    unsigned int d_numNodesForCurrentElement; /**< Number of nodes in the current element. */

private:
};
} // namespace Operator
} // namespace AMP

#endif
