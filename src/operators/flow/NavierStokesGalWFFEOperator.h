
#ifndef included_AMP_NavierStokesGalWFFEOperator
#define included_AMP_NavierStokesGalWFFEOperator

// AMP files
#include "ampmesh/MeshElement.h"
#include "discretization/DOF_Manager.h"
#include "vectors/MultiVariable.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"

#include "operators/flow/NavierStokesConstants.h"
#include "operators/flow/NavierStokesGalWFElement.h"
#include "operators/flow/NavierStokesGalWFFEOperatorParameters.h"
#include "operators/libmesh/NonlinearFEOperator.h"
#include "operators/libmesh/NonlinearFEOperator.h"

#include <vector>

namespace AMP {
namespace Operator {

class NavierStokesGalWFFEOperator : public NonlinearFEOperator
{
public:
    explicit NavierStokesGalWFFEOperator(
        const AMP::shared_ptr<NavierStokesGalWFFEOperatorParameters> &params );

    virtual ~NavierStokesGalWFFEOperator() {}

    void preAssembly( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                      AMP::shared_ptr<AMP::LinearAlgebra::Vector>
                          r ) override;

    void postAssembly() override;

    void preElementOperation( const AMP::Mesh::MeshElement & ) override;

    void postElementOperation() override;

    void reset( const AMP::shared_ptr<OperatorParameters> & ) override;

    AMP::shared_ptr<OperatorParameters>
    getJacobianParameters( AMP::LinearAlgebra::Vector::const_shared_ptr u_in ) override;

    void init();

    void setVector( unsigned int id, AMP::LinearAlgebra::Vector::shared_ptr frozenVec );

    static AMP::LinearAlgebra::Variable::shared_ptr createInputVariable( const std::string &name,
                                                                         int varId = -1 );

    static AMP::LinearAlgebra::Variable::shared_ptr createOutputVariable( const std::string &name,
                                                                          int varId = -1 )
    {
        (void) varId;
        AMP::LinearAlgebra::Variable::shared_ptr outVar( new AMP::LinearAlgebra::Variable( name ) );
        return outVar;
    }

    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() override { return d_inpVariables; }

    AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() override { return d_outVariables; }

    unsigned int numberOfDOFMaps();

    AMP::LinearAlgebra::Variable::shared_ptr getVariableForDOFMap( unsigned int id );

protected:
    void gettype0DofIndicesForCurrentElement( int varId, std::vector<std::vector<size_t>> &dofIds );
    void gettype1DofIndicesForCurrentElement( int varId, std::vector<size_t> &dofIds );

    std::vector<std::vector<size_t>> d_type0DofIndices; /**< Primary DOF indices */
    std::vector<size_t> d_type1DofIndices;

    unsigned int d_numNodesForCurrentElement;

    std::vector<double> d_elementOutputVector;

    AMP::shared_ptr<NavierStokesGalWFElement> d_flowGalWFElem;

    AMP::shared_ptr<FlowTransportModel> d_transportModel; /**< Flow Transport model. */

    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> d_inVec; /**< Input vector. */

    AMP::LinearAlgebra::Vector::shared_ptr d_referenceTemperature;

    AMP::LinearAlgebra::Vector::shared_ptr d_outVec;

    std::vector<bool> d_isActive;

    std::vector<bool> d_isFrozen;

    std::vector<AMP::Mesh::MeshElement> d_currNodes;

    bool d_coupledFormulation;

    AMP::shared_ptr<AMP::Discretization::DOFManager>
        d_dofMap[NavierStokes::TOTAL_NUMBER_OF_VARIABLES];

private:
    bool d_isInitialized; /**< A flag that is true if init() has been called and false otherwsie. */

    AMP::shared_ptr<AMP::LinearAlgebra::MultiVariable> d_inpVariables; /**< Input variables. */

    AMP::shared_ptr<AMP::LinearAlgebra::MultiVariable> d_outVariables; /**< Output variable. */
};
}
}

#endif
