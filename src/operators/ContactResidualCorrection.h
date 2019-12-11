
#ifndef included_AMP_ContactResidualCorrection
#define included_AMP_ContactResidualCorrection

#include "AMP/ampmesh/MeshID.h"
#include "AMP/operators/Operator.h"
#include "AMP/vectors/MultiVariable.h"

namespace AMP {
namespace Operator {

typedef OperatorParameters ContactResidualCorrectionParameters;

class ContactResidualCorrection : public Operator
{
public:
    explicit ContactResidualCorrection(
        const std::shared_ptr<ContactResidualCorrectionParameters> &params )
        : Operator( params )
    {
    }

    virtual ~ContactResidualCorrection() {}

    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    void reset( const std::shared_ptr<OperatorParameters> & ) override {}

    void setMasterVariable( const AMP::LinearAlgebra::Variable::shared_ptr &var )
    {
        d_masterVariable = var;
    }

    void setSlaveVariable( const AMP::LinearAlgebra::Variable::shared_ptr &var )
    {
        d_slaveVariable = var;
    }

    void setMasterMesh( const AMP::Mesh::Mesh::shared_ptr &mesh ) { d_Mesh = mesh; }

    void setSlaveMesh( const AMP::Mesh::Mesh::shared_ptr &mesh ) { d_slaveMesh = mesh; }

    void setMasterNodes( const std::vector<AMP::Mesh::MeshElementID> &vec ) { d_masterNodes = vec; }

    void setSlaveNodes( const std::vector<AMP::Mesh::MeshElementID> &vec ) { d_slaveNodes = vec; }

    void setDofs( const std::vector<std::vector<unsigned int>> &vec ) { d_dofs = vec; }

    AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() override
    {
        std::shared_ptr<AMP::LinearAlgebra::MultiVariable> retVariable(
            new AMP::LinearAlgebra::MultiVariable( "ContactVariable" ) );
        retVariable->add( d_masterVariable );
        retVariable->add( d_slaveVariable );
        return retVariable;
    }

    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() override
    {
        std::shared_ptr<AMP::LinearAlgebra::MultiVariable> retVariable(
            new AMP::LinearAlgebra::MultiVariable( "ContactVariable" ) );
        retVariable->add( d_masterVariable );
        retVariable->add( d_slaveVariable );
        return retVariable;
    }

    AMP::LinearAlgebra::Variable::shared_ptr getMasterVariable() { return d_masterVariable; }

    AMP::LinearAlgebra::Variable::shared_ptr getSlaveVariable() { return d_slaveVariable; }

    std::vector<AMP::Mesh::MeshElementID> getMasterNodes() { return d_masterNodes; }

    std::vector<AMP::Mesh::MeshElementID> getSlaveNodes() { return d_slaveNodes; }

    std::vector<std::vector<unsigned int>> getDofs() { return d_dofs; }

private:
    AMP::LinearAlgebra::Variable::shared_ptr d_masterVariable;
    AMP::LinearAlgebra::Variable::shared_ptr d_slaveVariable;

    std::vector<AMP::Mesh::MeshElementID> d_masterNodes;
    std::vector<AMP::Mesh::MeshElementID> d_slaveNodes;

    std::vector<std::vector<unsigned int>> d_dofs;

    AMP::Mesh::Mesh::shared_ptr d_slaveMesh;
};
} // namespace Operator
} // namespace AMP

#endif
