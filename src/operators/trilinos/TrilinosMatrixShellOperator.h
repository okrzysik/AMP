#ifndef included_TrilinosMatrixShellOperator
#define included_TrilinosMatrixShellOperator

#include "operators/LinearOperator.h"

#include "ml_include.h"


namespace AMP {
namespace Operator {


class TrilinosMatrixShellOperator : public LinearOperator
{
public:
    explicit TrilinosMatrixShellOperator( const AMP::shared_ptr<OperatorParameters> &params );

    virtual ~TrilinosMatrixShellOperator() {}

    void setOperator( AMP::shared_ptr<Operator> op );

    void setNodalDofMap( AMP::shared_ptr<AMP::Discretization::DOFManager> dofMap );

    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr f ) override;
    /**
     * Column specific implementation of the residual: f-L(u)
     * \param f: shared pointer to const vector rhs
     * \param u: shared pointer to const vector u
     * \param r: shared pointer to vector residual
     */
    virtual void residual( AMP::LinearAlgebra::Vector::const_shared_ptr f,
                           AMP::LinearAlgebra::Vector::const_shared_ptr u,
                           AMP::LinearAlgebra::Vector::shared_ptr r ) override;

    void reset( const AMP::shared_ptr<OperatorParameters> &params ) override;

    AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() override;

    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() override;

    static int
    matVec( ML_Operator *data, int in_length, double in[], int out_length, double out[] );

    static int getRow( ML_Operator *data,
                       int N_requested_rows,
                       int requested_rows[],
                       int allocated_space,
                       int columns[],
                       double values[],
                       int row_lengths[] );

    void setGetRow( void ( *func )(
        void *object, int row, std::vector<size_t> &cols, std::vector<double> &values ) );

    void getColumn( int column, std::vector<size_t> &rows, std::vector<double> &values );

    size_t getMatrixSize();

private:
    AMP::shared_ptr<AMP::Discretization::DOFManager> d_nodalDofMap;

    AMP::shared_ptr<Operator> d_operator;

    void ( *d_getRow )( void *object,
                        int row,
                        std::vector<size_t> &cols,
                        std::vector<double> &values );
};
} // namespace Operator
} // namespace AMP

#endif
