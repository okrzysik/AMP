#include "AMP/solvers/BandedSolver.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/operators/LinearOperator.h"
#include "AMP/utils/Utilities.h"
#include "ProfilerApp.h"

// External includes
#ifdef AMP_USE_LAPACK_WRAPPERS
    #include "LapackWrappers.h"
#endif


namespace AMP::Solver {


using AMP::Utilities::stringf;


BandedSolver::BandedSolver( std::shared_ptr<SolverStrategyParameters> parameters )
    : SolverStrategy( parameters ), N( 0 ), M( 0 ), KL( 0 ), KU( 0 ), AB( nullptr ), IPIV( nullptr )
{
    reset( parameters );
}


BandedSolver::~BandedSolver()
{
    delete[] AB;
    delete[] IPIV;
}


void BandedSolver::reset( std::shared_ptr<SolverStrategyParameters> parameters )
{
    PROFILE_START( "reset" );

    // Reset the parameters
    rightDOF.reset();
    leftDOF.reset();
    delete[] AB;
    delete[] IPIV;
    AB   = nullptr;
    IPIV = nullptr;
    if ( parameters != nullptr ) {
        KL = parameters->d_db->getScalar<int>( "KL" );
        KU = parameters->d_db->getScalar<int>( "KU" );
    }

    // Get the linear operator
    auto linear_op = std::dynamic_pointer_cast<AMP::Operator::LinearOperator>( d_pOperator );
    AMP_INSIST( linear_op, "ERROR: BandedSolver requires a linear operator" );

    // Get the matrix
    auto matrix = linear_op->getMatrix();
    AMP_INSIST( matrix, "ERROR: BandedSolver requires a matrix" );
    rightDOF = matrix->getRightDOFManager();
    leftDOF  = matrix->getLeftDOFManager();
    M        = static_cast<int>( matrix->numLocalRows() );
    N        = static_cast<int>( matrix->numLocalColumns() );

    // Allocate space
    int K  = 2 * KL + KU + 1;
    int N2 = std::max( N, M );
    AB     = new double[K * N2];
    IPIV   = new int[N2];

    // Copy the matrix to a banded diagonal form
    memset( AB, 0, K * N2 * sizeof( double ) );
    size_t row_begin = leftDOF->beginDOF();
    size_t col_begin = rightDOF->beginDOF();
    size_t col_end   = rightDOF->endDOF();
    std::vector<size_t> cols;
    std::vector<double> values;
    for ( int i = 0; i < M; i++ ) {
        size_t row = i + row_begin;
        matrix->getRowByGlobalID( row, cols, values );
        for ( size_t k = 0; k < cols.size(); k++ ) {
            if ( values[k] == 0 )
                continue;
            if ( cols[k] < col_begin || cols[k] >= col_end )
                AMP_ERROR( "Matrix has entries that are non-local" );
            int j = cols[k] - col_begin;
            if ( j < i - KL || j > i + KU || j < 0 || j >= N ) {
                auto msg = stringf( "Banded entry is out of bounds (%i,%i,%e)", i, j, values[k] );
                AMP_ERROR( msg );
            }
            AB[KL + KU + i - j + j * K] = values[k];
        }
        if ( AB[KL + KU + i * K] == 0 ) {
            auto msg = stringf( "Error diagonal entry M(%i,%i) = 0", i + 1, i + 1 );
            AMP_ERROR( msg );
        }
    }

    // Factor the matrix
#ifdef AMP_USE_LAPACK_WRAPPERS
    int error = 0;
    Lapack<double>::gbtrf( M, N, KL, KU, AB, K, IPIV, error );
    if ( error != 0 ) {
        auto msg = stringf( "Error factoring matrix (%i)", error );
        AMP_ERROR( msg );
    }
#else
    AMP_ERROR( "Lapack Required" );
#endif

    PROFILE_STOP( "reset" );
}


void BandedSolver::resetOperator( std::shared_ptr<const AMP::Operator::OperatorParameters> params )
{
    PROFILE_START( "resetOperator" );
    AMP_INSIST( ( d_pOperator ), "ERROR: BandedSolver::resetOperator() operator cannot be NULL" );
    d_pOperator->reset( params );
    reset( std::shared_ptr<SolverStrategyParameters>() );
    PROFILE_STOP( "resetOperator" );
}


void BandedSolver::apply( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                          std::shared_ptr<AMP::LinearAlgebra::Vector> u )
{
    PROFILE_START( "solve" );

    // Copy f
    AMP_ASSERT( *rightDOF == *( f->getDOFManager() ) );
    auto B = new double[N];
    f->getRawData( B );

    // Solve the
#ifdef AMP_USE_LAPACK_WRAPPERS
    int error = 0;
    Lapack<double>::gbtrs( 'N', N, KL, KU, 1, AB, 2 * KL + KU + 1, IPIV, B, N, error );
    d_iNumberIterations = 1;
    if ( error != 0 ) {
        auto msg = stringf( "Error solving matrix (%i)", error );
        AMP_ERROR( msg );
    }
#else
    AMP_ERROR( "Lapack Required" );
#endif

    // Copy the solution
    AMP_ASSERT( *rightDOF == *( u->getDOFManager() ) );
    u->putRawData( B );
    delete[] B;

    PROFILE_STOP( "solve" );
}


} // namespace AMP::Solver
