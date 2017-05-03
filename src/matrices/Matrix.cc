
#include "utils/ParameterBase.h"
#include "Matrix.h"
#include <iomanip>

namespace AMP {
namespace LinearAlgebra {

Matrix::shared_ptr Matrix::matMultiply( shared_ptr A, shared_ptr B )
{
    if ( A->numGlobalColumns() != B->numGlobalRows() )
        AMP_ERROR( "Inner matrix dimensions must agree" );
    shared_ptr retVal;
    A->multiply( B, retVal );
    return retVal;
}


// Get the number of local rows in the matrix
size_t Matrix::numLocalRows() const
{
    auto DOF = getLeftDOFManager();
    return DOF->numLocalDOF();
}


// Get the number of global rows in the matrix
size_t Matrix::numGlobalRows() const
{
    auto DOF = getLeftDOFManager();
    return DOF->numGlobalDOF();
}


// Get the number of local columns in the matrix
size_t Matrix::numLocalColumns() const
{
    auto DOF = getRightDOFManager();
    return DOF->numLocalDOF();
}


// Get the number of global columns in the matrix
size_t Matrix::numGlobalColumns() const
{
    auto DOF = getRightDOFManager();
    return DOF->numGlobalDOF();
}

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

// Print the matrix to a IO stream
std::ostream &operator<<( std::ostream &out, const Matrix &M_in )
{
    Matrix *M = (Matrix *) &M_in;
    // Print the matrix type (not supported yet)
    /*out << "Vector type: " << v.type() << "\n";
    if ( v.getVariable() )
    {
      out << "Variable name: " << v.getVariable()->getName() << "\n";
    }*/
    // Print the rank
    Discretization::DOFManager::shared_ptr leftDOF  = M->getLeftDOFManager();
    Discretization::DOFManager::shared_ptr rightDOF = M->getRightDOFManager();
    AMP_MPI leftComm                                = leftDOF->getComm();
    AMP_MPI rightComm                               = rightDOF->getComm();
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
    /*
        out << "Full Matix: " << std::endl;
        out << std::setprecision(15);
        for (size_t row=0; row<leftDOF->numGlobalDOF(); row++) {
            M->getRowByGlobalID( row, cols, values );
            std::vector<double> A(M->numGlobalColumns(),0.);
            for (size_t i=0; i<cols.size(); i++)
              A[cols[i]]=values[i];
            for (size_t i=0; i<A.size(); i++) out<< A[i]<<"  ";
            out<<std::endl;
        }
        out.unsetf(std::ios::floatfield);
    */
    return out;
}
}
}
