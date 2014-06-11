#include "solvers/BandedSolver.h"
#include "operators/LinearOperator.h"
#include "utils/Utilities.h"
#include "utils/ProfilerApp.h"
#include "utils/LapackWrappers.h"


namespace AMP {
namespace Solver {


BandedSolver::BandedSolver( boost::shared_ptr<SolverStrategyParameters> parameters ):
    SolverStrategy(parameters),
    AB(NULL), IPIV(NULL), KL(0), KU(0)
{
    KL = parameters->d_db->getInteger("KL");
    KU = parameters->d_db->getInteger("KU");
    reset( parameters );
}


BandedSolver::~BandedSolver()
{
    delete [] AB;
    delete [] IPIV;
}


void BandedSolver::reset( boost::shared_ptr<SolverStrategyParameters> parameters ) 
{
    PROFILE_START("reset");

    // Get the properties from the database
    if ( parameters != NULL ) {
        KL = parameters->d_db->getIntegerWithDefault("KL",KL);
        KU = parameters->d_db->getIntegerWithDefault("KU",KU);
    }
    int K = 2*KL+KU+1;

    // Get the matrix
    boost::shared_ptr<AMP::Operator::LinearOperator> linear_op = 
        boost::dynamic_pointer_cast<AMP::Operator::LinearOperator>(d_pOperator);
    AMP_INSIST(linear_op.get()!=NULL, "ERROR: BandedSolver requires a linear operator");
    AMP::LinearAlgebra::Matrix::shared_ptr matrix = linear_op->getMatrix();
    AMP_INSIST(linear_op.get()!=NULL, "ERROR: BandedSolver requires a matrix");
    rightDOF = matrix->getRightDOFManager();
    leftDOF = matrix->getLeftDOFManager();
    M = static_cast<int>(matrix->numLocalRows());
    N = static_cast<int>(matrix->numGlobalColumns());

    // Allocate space 
    delete [] AB;
    delete [] IPIV;
    int N2 = std::max(N,M);
    AB = new double[K*N2];
    IPIV = new int[N2];

    // Copy the matrix to a banded diagonal form
    memset(AB,0,K*N2*sizeof(double));
    size_t row_begin = leftDOF->beginDOF();
    size_t row_end   = leftDOF->endDOF();
    size_t col_begin = rightDOF->beginDOF();
    size_t col_end   = rightDOF->endDOF();
    std::vector<unsigned int> cols;
    std::vector<double> values;
    for (int i=0; i<M; i++) {
        size_t row = i + row_begin;
        matrix->getRowByGlobalID( row, cols, values );
        for (size_t k=0; k<cols.size(); k++) {
            if ( values[k]==0 )
                continue;
            if ( cols[k]<col_begin || cols[k]>=col_end )
                AMP_ERROR("Matrix has entries that are non-local");
            int j = cols[k] - col_begin;
            if ( j<i-KL || j>i+KU || j<0 || j>=N ) {
                char tmp[100];
                sprintf(tmp,"Banded entry is out of bounds (%i,%i,%e)",i,j,values[k]);
                AMP_ERROR(tmp);
            }
            AB[KL+KU+i-j+j*K] = values[k];
        }
        if ( AB[KL+KU+i*K]==0 ) {
            char msg[100];
            sprintf(msg,"Error diagonal entry M(%i,%i) = 0",i+1,i+1);
            AMP_ERROR(msg);
        }
    }
    
    // Factor the matrix
    int error = 0;
    Lapack::dgbtrf( M, N, KL, KU, AB, K, IPIV, error );
    if ( error != 0 ) {
        char msg[100];
        sprintf(msg,"Error factoring matrix (%i)",error);
        AMP_ERROR(msg);
    }


    PROFILE_STOP("reset");
}


void BandedSolver::resetOperator(const boost::shared_ptr<AMP::Operator::OperatorParameters> params) 
{
    PROFILE_START("resetOperator");
    AMP_INSIST((d_pOperator.get() != NULL), "ERROR: BandedSolver::resetOperator() operator cannot be NULL");
    d_pOperator->reset(params);
    reset( boost::shared_ptr<SolverStrategyParameters>() );
    PROFILE_STOP("resetOperator");
}


void BandedSolver::solve(boost::shared_ptr<const AMP::LinearAlgebra::Vector> f,
    boost::shared_ptr<AMP::LinearAlgebra::Vector> u )
{
    PROFILE_START("solve");
    
    // Copy f 
    AMP_ASSERT(*rightDOF==*(f->getDOFManager()));
    double *B = new double[N];
    f->copyOutRawData(B);

    // Solve the 
    int error = 0;
    Lapack::dgbtrs( 'N', N, KL, KU, 1, AB, 2*KL+KU+1, IPIV, B, N, error );
    d_iNumberIterations = 1;
    if ( error != 0 ) {
        char msg[100];
        sprintf(msg,"Error solving matrix (%i)",error);
        AMP_ERROR(msg);
    }

    // Copy the solution
    AMP_ASSERT(*rightDOF==*(u->getDOFManager()));
    u->putRawData(B);
    delete [] B;

    PROFILE_STOP("solve");
}



} // namespace Solver
} // namespace AMP
