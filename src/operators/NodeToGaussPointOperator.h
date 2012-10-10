#ifndef included_AMP_NodeToGaussPointOperator
#define included_AMP_NodeToGaussPointOperator

#include "operators/Operator.h"
#include "discretization/createLibmeshElements.h"


namespace AMP {
namespace Operator {

class NodeToGaussPointOperator : public Operator
{
public :

    NodeToGaussPointOperator (const boost::shared_ptr<OperatorParameters> & params);

    virtual ~NodeToGaussPointOperator() { }

    void apply(AMP::LinearAlgebra::Vector::const_shared_ptr f, AMP::LinearAlgebra::Vector::const_shared_ptr u,
        AMP::LinearAlgebra::Vector::shared_ptr r, const double a = -1.0, const double b = 1.0);

    virtual AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() { return d_GaussPtVariable; }

    virtual AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() {return d_NodalVariable; }

protected :
    bool d_UseSurfaceElements;
    AMP::LinearAlgebra::Variable::shared_ptr    d_NodalVariable;
    AMP::LinearAlgebra::Variable::shared_ptr    d_GaussPtVariable;
    AMP::Mesh::MeshIterator                     d_iterator;
    Discretization::createLibmeshElements       d_libmeshElements;

};


}
}

#endif

