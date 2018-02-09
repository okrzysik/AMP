#ifndef included_PetscMatrixShellOperator
#define included_PetscMatrixShellOperator

#include "AMP/operators/LinearOperator.h"
#include "AMP/utils/AMP_MPI.h"


#ifdef MPICH_SKIP_MPICXX
#define _FIX_FOR_PETSC_MPI_CXX
#undef MPICH_SKIP_MPICXX
#endif

#ifdef OMPI_SKIP_MPICXX
#define _FIX_FOR_PETSC_OMPI_CXX
#undef OMPI_SKIP_MPICXX
#endif

#include "petscmat.h"

#ifdef _FIX_FOR_PETSC_MPI_CXX
#ifndef MPICH_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#endif
#endif

#ifdef _FIX_FOR_PETSC_OMPI_CXX
#ifndef OMPI_SKIP_MPICXX
#define OMPI_SKIP_MPICXX
#endif
#endif


namespace AMP {
namespace Operator {


class PetscMatrixShellOperator : public LinearOperator
{
public:
    PetscMatrixShellOperator( const AMP::shared_ptr<OperatorParameters> &params );

    virtual ~PetscMatrixShellOperator() {}

    static PetscErrorCode mult( Mat, Vec, Vec );

    void setMatLocalRowSize( int val );

    void setMatLocalColumnSize( int val );

    void setComm( AMP_MPI comm );

    void setOperator( AMP::shared_ptr<Operator> op );

    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    void reset( const AMP::shared_ptr<OperatorParameters> &params ) override;

    AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() override;

    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() override;

private:
    AMP::shared_ptr<Operator> d_operator;
    int d_iMatLocalRowSize;
    int d_iMatLocalColumnSize;
    Mat d_mat;
    AMP_MPI d_comm;
};
} // namespace Operator
} // namespace AMP

#endif
