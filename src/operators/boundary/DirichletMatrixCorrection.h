
#ifndef included_AMP_DirichletMatrixCorrection
#define included_AMP_DirichletMatrixCorrection

#include "BoundaryOperator.h"
#include "DirichletMatrixCorrectionParameters.h"
#include "DirichletVectorCorrection.h"

#include "vectors/Variable.h"

namespace AMP {
namespace Operator {

/**
 *  A class used to impose Dirichlet boundary conditions for a linear operator. For a linear
 * operator, imposing
 *  Dirichlet boundary conditions involves the following steps:
 *  1) Modify the entries of the matrix appropriately.
 *  2) Add a vector of corrections to the RHS vector
 *  3) Set the Dirichlet values at the appropriate locations in the RHS vector.
 */
class DirichletMatrixCorrection : public BoundaryOperator
{
public:
    //! Constructor
    explicit DirichletMatrixCorrection(
        const AMP::shared_ptr<DirichletMatrixCorrectionParameters> &params );

    //! Destructor
    virtual ~DirichletMatrixCorrection() {}

    //! Set the variable for the vector that will used with this operator.
    void setVariable( const AMP::LinearAlgebra::Variable::shared_ptr &var ) { d_variable = var; }

    AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() override { return d_variable; }

    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() override { return d_variable; }

    virtual void apply( AMP::LinearAlgebra::Vector::const_shared_ptr,
                        AMP::LinearAlgebra::Vector::shared_ptr ) override
    {
        // Do Nothing
    }

    void parseParams( const AMP::shared_ptr<DirichletMatrixCorrectionParameters> & );

    /**
      This function modifies the entries of the matrix formed by the volume operator
      in order to impose Dirichlet boundary conditions. This function can also be used
      to change the Dirichlet boundary conditions, if required.
      */
    void reset( const AMP::shared_ptr<OperatorParameters> & ) override;

    /**
      Adds a vector to the RHS vector. This is one of the steps for imposing Dirichlet boundary
      conditions.
      This step can be skipped if the Dirichlet boundary conditons are homogeneous.
      */
    void addRHScorrection( AMP::LinearAlgebra::Vector::shared_ptr rhs ) override;

    /**
      Sets the Dirichlet values at the appropriate locations in the RHS vector. This is one
      of the steps for imposing Dirichlet boundary conditions.
      */
    void setRHScorrection( AMP::LinearAlgebra::Vector::shared_ptr rhs ) override;

    std::vector<short int> getBoundaryIds() { return d_boundaryIds; }

    std::vector<std::vector<unsigned int>> getDofIds() { return d_dofIds; }

protected:
    // This must be a simple variable not a dual or multivariable
    AMP::LinearAlgebra::Variable::shared_ptr d_variable;

    std::vector<short int> d_boundaryIds;

    std::vector<std::vector<double>> d_dirichletValues;

    std::vector<std::vector<unsigned int>> d_dofIds;

    AMP::LinearAlgebra::Vector::shared_ptr d_rhsCorrectionAdd;

    AMP::LinearAlgebra::Vector::shared_ptr d_dispVals;

    AMP::shared_ptr<DirichletVectorCorrection> d_rhsCorrectionSet;

    bool d_symmetricCorrection;

    bool d_zeroDirichletBlock;

    bool d_skipRHSaddCorrection;

    bool d_skipRHSsetCorrection;

    bool d_computedAddRHScorrection;

    void initRhsCorrectionSet();

    void initRhsCorrectionAdd( AMP::LinearAlgebra::Vector::shared_ptr rhs );

    void applyMatrixCorrection();

    bool d_applyMatrixCorrectionWasCalled;

    AMP::LinearAlgebra::Matrix::shared_ptr d_inputMatrix;

private:
};
}
}

#endif
