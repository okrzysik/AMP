#include "AMP/matrices/Matrix.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/ParameterBase.h"
#include <iomanip>

namespace AMP::LinearAlgebra {


/********************************************************
 * Constructors                                          *
 ********************************************************/
Matrix::Matrix( const Matrix &rhs ) : d_comm( rhs.d_comm )
{
    AMPManager::incrementResource( "Matrix" );
}
Matrix::Matrix() { AMPManager::incrementResource( "Matrix" ); }
Matrix::Matrix( std::shared_ptr<MatrixParameters> params ) : d_comm( params->getComm() )
{
    AMP_ASSERT( !d_comm.isNull() );
    AMPManager::incrementResource( "Matrix" );
}
Matrix::~Matrix() { AMPManager::decrementResource( "Matrix" ); }


/********************************************************
 * Get the number of rows/columns in the matrix          *
 ********************************************************/
size_t Matrix::numLocalRows() const
{
    auto DOF = getLeftDOFManager();
    return DOF->numLocalDOF();
}
size_t Matrix::numGlobalRows() const
{
    auto DOF = getLeftDOFManager();
    return DOF->numGlobalDOF();
}
size_t Matrix::numLocalColumns() const
{
    auto DOF = getRightDOFManager();
    return DOF->numLocalDOF();
}
size_t Matrix::numGlobalColumns() const
{
    auto DOF = getRightDOFManager();
    return DOF->numGlobalDOF();
}


/********************************************************
 * Get iterators                                         *
 ********************************************************/
size_t Matrix::beginRow() const
{
    auto DOF = getRightDOFManager();
    return DOF->beginDOF();
}
size_t Matrix::endRow() const
{
    auto DOF = getRightDOFManager();
    return DOF->endDOF();
}


/********************************************************
 * multiply                                             *
 ********************************************************/
std::shared_ptr<Matrix> Matrix::matMultiply( shared_ptr A, shared_ptr B )
{
    if ( A->numGlobalColumns() != B->numGlobalRows() )
        AMP_ERROR( "Inner matrix dimensions must agree" );
    shared_ptr retVal;
    A->multiply( B, retVal );
    return retVal;
}


/********************************************************
 * axpy                                                  *
 ********************************************************/
void Matrix::axpy( double alpha, std::shared_ptr<const Matrix> x )
{
    AMP_ASSERT( x );
    size_t N1 = x->numGlobalColumns();
    size_t N2 = this->numGlobalRows();
    if ( N1 != N2 )
        AMP_ERROR( "Matrix sizes are not compatible" );
    axpy( alpha, *x );
}


/********************************************************
 * Print the matrix to a IO stream                       *
 ********************************************************/
std::ostream &operator<<( std::ostream &out, const Matrix &M_in )
{
    auto *M = (Matrix *) &M_in;
    // Print the matrix type (not supported yet)
    /*out << "Vector type: " << v.type() << "\n";
    if ( v.getVariable() )
    {
      out << "Variable name: " << v.getName() << "\n";
    }*/
    // Print the rank
    auto leftDOF   = M->getLeftDOFManager();
    auto rightDOF  = M->getRightDOFManager();
    auto leftComm  = leftDOF->getComm();
    auto rightComm = rightDOF->getComm();
    if ( leftComm == rightComm ) {
        int rank = leftComm.getRank();
        out << "Processor: " << rank << "\n";
    } else {
        int leftRank  = leftComm.getRank();
        int rightRank = rightComm.getRank();
        out << "Processor (left comm):  " << leftRank << "\n";
        out << "Processor (right comm): " << rightRank << "\n";
    }
    // Print some basic matrix info
    out << "\n"
        << "Global number of rows: " << M->numGlobalRows() << "\n"
        << "Global number of colums: " << M->numGlobalColumns() << "\n"
        << "Local number of rows: " << M->numLocalRows() << "\n"
        << "Local number of colums: " << M->numLocalColumns() << "\n";
    // Loop through each local row
    std::vector<size_t> cols;
    std::vector<double> values;
    out << "Compressed Matix: " << std::endl;
    for ( size_t row = leftDOF->beginDOF(); row < leftDOF->endDOF(); row++ ) {
        M->getRowByGlobalID( row, cols, values );
        out << "Row " << row << " (" << cols.size() << " entries):"
            << "\n";
        for ( size_t i = 0; i < cols.size(); i++ )
            out << "    M(" << row << "," << cols[i] << ") = " << values[i] << "\n";
    }
    return out;
}
} // namespace AMP::LinearAlgebra
