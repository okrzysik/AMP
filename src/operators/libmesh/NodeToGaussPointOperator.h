#ifndef included_AMP_NodeToGaussPointOperator
#define included_AMP_NodeToGaussPointOperator

#include "discretization/createLibmeshElements.h"
#include "operators/Operator.h"


namespace AMP {
namespace Operator {

class NodeToGaussPointOperator : public Operator
{
public:
    explicit NodeToGaussPointOperator( const AMP::shared_ptr<OperatorParameters> &params );

    virtual ~NodeToGaussPointOperator() {}

    virtual void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                        AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    virtual AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable()
    {
        return d_GaussPtVariable;
    }

    virtual AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() { return d_NodalVariable; }

protected:
    bool d_UseSurfaceElements;
    int d_dim;
    AMP::LinearAlgebra::Variable::shared_ptr d_NodalVariable;
    AMP::LinearAlgebra::Variable::shared_ptr d_GaussPtVariable;
    AMP::Mesh::MeshIterator d_iterator;
    Discretization::createLibmeshElements d_libmeshElements;
    std::vector<std::vector<AMP::Mesh::MeshElementID>> d_nodes;
    std::vector<unsigned short int> d_N_quad;
    std::vector<std::vector<double>> d_phi;
};
}
}

#endif
