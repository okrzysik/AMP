
#ifndef included_AMP_VolumeIntegralOperator
#define included_AMP_VolumeIntegralOperator

#include "matrices/Matrix.h"
#include "operators/libmesh/NonlinearFEOperator.h"
#include "operators/libmesh/SourceNonlinearElement.h"
#include "operators/libmesh/VolumeIntegralOperatorParameters.h"
#include "vectors/MultiVariable.h"

#include <vector>

namespace AMP {
namespace Operator {


/**
  This class is used to integrate any arbitrary functional over the entire domain. The
  evaluation of the functional at quadrature points and volume integration are its key
  basic operations.
  */
class VolumeIntegralOperator : public NonlinearFEOperator
{
public:
    /**
      Constructor. This reads the values for the following keys from the database object contained
      in
      the parameter object, params:
      1) Primary(Active)InputVariables - List of active input variables names. The supported
      variable types are:
      NodalScalar, Nodalvector, IntegrationPointScalar, IntegrationPointVector.
      2) AuxillaryInputVariables - List of auxillary input variables names. These are frozen
      variables and are
      temporarily not used in any formulation.
      3) OutputVariable - Name of the output variable
      */
    explicit VolumeIntegralOperator(
        const AMP::shared_ptr<VolumeIntegralOperatorParameters> &params );

    /**
      Destructor.
      */
    virtual ~VolumeIntegralOperator() {}

    /**
      This is used to update the operator between successive solves with the operator.
      */
    void reset( const AMP::shared_ptr<OperatorParameters> & ) override;

    void setAuxVariable( const AMP::shared_ptr<AMP::LinearAlgebra::MultiVariable> &var )
    {
        d_auxVariables = var;
    }

    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() override { return d_inpVariables; }

    AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() override { return d_outVariable; }

    AMP::shared_ptr<SourcePhysicsModel> getSourcePhysicsModel() { return d_sourcePhysicsModel; }

    AMP::shared_ptr<SourceNonlinearElement> getSourceElement() { return d_srcNonlinElem; }

    /*
       AMP::shared_ptr<AMP::LinearAlgebra::Matrix> getLinearizedVolumeIntegralOperator(
       const AMP::shared_ptr<OperatorParameters>& params);
     */

protected:
    /**
      This is used to compute the information required to reset the corresponding Linear (Jacobian)
      operator
      */
    AMP::shared_ptr<OperatorParameters>
    getJacobianParameters( AMP::LinearAlgebra::Vector::const_shared_ptr u ) override;


    /**
      This function is called at the beginning of the FE assembly
      */
    void preAssembly( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                      AMP::LinearAlgebra::Vector::shared_ptr r ) override;

    /**
      This function is called at the end of the FE assembly
      */
    void postAssembly() override;

    /**
      This function is called at the beginning of the element computation
      */
    void preElementOperation( const AMP::Mesh::MeshElement & ) override;

    /**
      This function is called at the end of the element computation
      */
    void postElementOperation() override;

    void getNodeDofIndicesForCurrentElement();

    void init( const AMP::shared_ptr<VolumeIntegralOperatorParameters> &params );

    std::string d_isInputType;

    std::vector<double> d_elementOutputVector; /**< Output vector for the Element Operation. */

    AMP::shared_ptr<SourceNonlinearElement> d_srcNonlinElem; /**< Element Operation. */

    AMP::shared_ptr<SourcePhysicsModel> d_sourcePhysicsModel; /**< Source Physics Model. */

    std::vector<AMP::LinearAlgebra::Vector::const_shared_ptr>
        d_inVec; /**< Input vector for active variables. */

    std::vector<AMP::LinearAlgebra::Vector::const_shared_ptr>
        d_auxVec; /**< Input vector for auxillary variables. */

    AMP::LinearAlgebra::Vector::shared_ptr d_multiAuxPtr;

    AMP::LinearAlgebra::Vector::shared_ptr d_outVec; /**< Output vector. */

    std::vector<size_t> d_dofIndices; /**<Node dof Ids */

    AMP::Discretization::DOFManager::shared_ptr d_elementDofMap;
    AMP::Discretization::DOFManager::shared_ptr d_nodeDofMap;

    std::vector<AMP::Mesh::MeshElement> d_currNodes;

private:
    AMP::shared_ptr<AMP::LinearAlgebra::MultiVariable>
        d_inpVariables; /**< Input Active variables. */

    AMP::shared_ptr<AMP::LinearAlgebra::MultiVariable>
        d_auxVariables; /**< Input Auxillary variables. */

    AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariable; /**< Output variables. */

    // AMP::shared_ptr<AMP::LinearAlgebra::Matrix> d_pDiagonalMatrix;
    // AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_pDiagonalVector;
    // AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_pNullVector;
    // bool d_bMatrixAndVectorsCloned;
};


} // namespace Operator
} // namespace AMP

#endif
