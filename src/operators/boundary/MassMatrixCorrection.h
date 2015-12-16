#ifndef included_AMP_MassMatrixCorrection
#define included_AMP_MassMatrixCorrection

#include "BoundaryOperator.h"
#include "DirichletMatrixCorrectionParameters.h"

#include "vectors/Variable.h"

namespace AMP {
namespace Operator {


typedef DirichletMatrixCorrectionParameters MassMatrixCorrectionParameters;


/**
  A class used to impose HOMOGENEOUS Dirichlet boundary conditions for a mass operator. This is a
  minor modification of
  the class
  DirichletMatrixCorrection. A feature to handle nonhomogeneous boundary condition needs to be
  implemented soon.
 */
class MassMatrixCorrection : public BoundaryOperator
{
public:
    //! Constructor
    explicit MassMatrixCorrection( const AMP::shared_ptr<MassMatrixCorrectionParameters> &params )
        : BoundaryOperator( params ), d_bSetIdentityOnDiagonal( false )
    {
        d_variable = params->d_variable;
        reset( params );
    }

    //! Destructor
    virtual ~MassMatrixCorrection() {}

    /**
      Set the variable for the vector that will used with this operator.
      */
    void setVariable( const AMP::LinearAlgebra::Variable::shared_ptr &var ) { d_variable = var; }

    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr,
                AMP::LinearAlgebra::Vector::shared_ptr ) override
    {
        // Do Nothing
    }

    void resetBoundaryIds( const AMP::shared_ptr<MassMatrixCorrectionParameters> & );

    /**
      This function modifies the entries of the matrix formed by the volume operator
      in order to impose Dirichlet boundary conditions. This function can also be used
      to change the Dirichlet boundary conditions, if required.
     */
    void reset( const AMP::shared_ptr<OperatorParameters> & ) override;

protected:
    // This must be a simple variable not a dual or multivariable
    AMP::LinearAlgebra::Variable::shared_ptr d_variable;

    bool d_bSetIdentityOnDiagonal;

    std::vector<short int> d_boundaryIds;

    std::vector<std::vector<unsigned int>> d_dofIds;


private:
};
}
}

#endif
