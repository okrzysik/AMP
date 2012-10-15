#include "matrices/MatrixParameters.h"
#include "utils/Utilities.h"
#include "utils/AMP_MPI.h"


namespace AMP {
namespace LinearAlgebra {



MatrixParameters::MatrixParameters( AMP::Discretization::DOFManager::shared_ptr left, AMP::Discretization::DOFManager::shared_ptr right, AMP_MPI comm )
{
    d_comm = comm;
    d_DOFManagerLeft = left;
    d_DOFManagerRight = right;
    if ( d_DOFManagerLeft!=NULL ) {
        AMP_ASSERT(d_comm>=d_DOFManagerLeft->getComm());
        numLocalRows = d_DOFManagerLeft->numLocalDOF();
    }
    if ( d_DOFManagerRight!=NULL ) {
        AMP_ASSERT(d_comm>=d_DOFManagerRight->getComm());
        numLocalColumns = d_DOFManagerRight->numLocalDOF();
    }
    numGlobalRows = comm.sumReduce(numLocalRows);
    numGlobalColumns = comm.sumReduce(numLocalColumns);
}



}
}

