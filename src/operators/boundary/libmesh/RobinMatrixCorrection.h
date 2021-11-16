
#ifndef included_AMP_RobinMatrixCorrection
#define included_AMP_RobinMatrixCorrection

#include "AMP/operators/boundary/BoundaryOperator.h"
#include "AMP/operators/boundary/libmesh/NeumannVectorCorrection.h"
#include "AMP/operators/boundary/libmesh/NeumannVectorCorrectionParameters.h"
#include "AMP/operators/boundary/libmesh/RobinMatrixCorrectionParameters.h"

#include "AMP/discretization/createLibmeshElements.h"

// Libmesh files
#include "libmesh/quadrature.h"

#include <string>

namespace AMP {
namespace Operator {

/**
  A class to impose Robin Boundary conditions for a linear operator. Robin Condition
  is also known as mixed condition of Dirichlet and Neumann Flux conditions. This can
  be written as \f$\alpha k*\frac{\partial u}{\partial n} + \beta h*u = \gamma*c \f$.
  Imposing this condition would involve:
  1) Imposing a Neumann Flux condition on the RHS Vector
  2) Make appropriate matrix corrections on the boundary nodes.
  */
class RobinMatrixCorrection : public BoundaryOperator
{
public:
    /**
       Constructor. This function reads all the parameters required for surface elements.
       This also constructs new NeumannVectorCorrection parameters and calls it reset.
    */
    explicit RobinMatrixCorrection( std::shared_ptr<const RobinMatrixCorrectionParameters> params );

    /**
       Set the variable for the vector that will used with this operator.
    */
    void setVariable( const std::shared_ptr<AMP::LinearAlgebra::Variable> &var )
    {
        d_variable = var;
    }

    /**
       Destructor
    */
    virtual ~RobinMatrixCorrection() {}

    //! Return the name of the operator
    std::string type() const override { return "RobinMatrixCorrection"; }

    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr,
                AMP::LinearAlgebra::Vector::shared_ptr ) override
    {
        // Do Nothing
    }

    /**
       This function reads parameters related to boundary Ids. Since this class allow
       for variable flux values, the parameters stores a vector of values. This vector
       is passed to the NeumannVectorCorrection. This function also does a matrix
       correction on the boundary.
    */
    void reset( std::shared_ptr<const OperatorParameters> params ) override;

    /**
       Adds a Neumann Correction Vector to the RHS vector.
    */
    void addRHScorrection( AMP::LinearAlgebra::Vector::shared_ptr rhs ) override
    {
        d_NeumannCorrection->addRHScorrection( rhs );
    }


protected:
    Discretization::createLibmeshElements libmeshElements;

    std::vector<short int> d_boundaryIds;

    std::vector<std::vector<unsigned int>> d_dofIds;

    std::vector<std::vector<double>> d_robinValues;

    std::shared_ptr<AMP::LinearAlgebra::Variable> d_variable;

    double d_hef; // Convective Coefficient

    double d_alpha; // pre-factor solid flux

    double d_beta;

    double d_gamma;

    const std::vector<libMesh::Real> *d_JxW;

    const std::vector<std::vector<libMesh::Real>> *d_phi;

    std::shared_ptr<libMesh::FEType> d_feType;

    std::shared_ptr<libMesh::FEBase> d_fe;

    std::shared_ptr<libMesh::QBase> d_qrule;

    std::string d_qruleOrderName;

    libMeshEnums::Order d_feTypeOrder;
    libMeshEnums::FEFamily d_feFamily;
    libMeshEnums::QuadratureType d_qruleType;
    libMeshEnums::Order d_qruleOrder;

    AMP::LinearAlgebra::Vector::shared_ptr d_Frozen;

    std::shared_ptr<RobinPhysicsModel> d_robinPhysicsModel;

    std::shared_ptr<AMP::Discretization::DOFManager> d_dofManager;

private:
    std::shared_ptr<NeumannVectorCorrection> d_NeumannCorrection;
    std::shared_ptr<NeumannVectorCorrectionParameters> d_NeumannParams;
};
} // namespace Operator
} // namespace AMP

#endif
