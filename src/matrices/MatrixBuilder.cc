#ifdef USE_AMP_VECTORS

#include "MatrixBuilder.h"

#include "discretization/DOF_Manager.h"
#include "vectors/petsc/ManagedPetscVector.h"
#include "vectors/trilinos/EpetraVectorEngine.h"
#include "matrices/petsc/ManagedPetscMatrix.h"


namespace AMP {
namespace LinearAlgebra {


/********************************************************
* Matrix builder                                        *
********************************************************/
AMP::LinearAlgebra::Matrix::shared_ptr  createMatrix( 
    AMP::LinearAlgebra::Vector::shared_ptr operandVec, 
    AMP::LinearAlgebra::Vector::shared_ptr resultVec )
{
    // Get the DOFs
    AMP::Discretization::DOFManager::shared_ptr operandDOF = operandVec->getDOFManager();
    AMP::Discretization::DOFManager::shared_ptr resultDOF = resultVec->getDOFManager();
    if ( operandDOF->getComm().compare(resultVec->getComm()) == 0 )
        AMP_ERROR("operandDOF and resultDOF on different comm groups is NOT tested, and needs to be fixed");
    AMP_MPI comm = operandDOF->getComm();

    // Create the matrix parameters
    boost::shared_ptr<AMP::LinearAlgebra::ManagedPetscMatrixParameters> params( 
        new AMP::LinearAlgebra::ManagedPetscMatrixParameters ( resultDOF->numLocalDOF(),
        resultDOF->numGlobalDOF(), 0, operandDOF->numGlobalDOF(), 0, comm ) );

    // Add the rows to the matrix parameters
    AMP::Mesh::MeshIterator cur_elem = resultDOF->getIterator();
    AMP::Mesh::MeshIterator end_elem = cur_elem.end();
    int columns[1000];   // Preallocate for the columns for speed
    while ( cur_elem != end_elem) {
        AMP::Mesh::MeshElement obj = *cur_elem;
        // Get the result DOFs associated with the given element
        std::vector<size_t> ids;
        resultDOF->getDOFs(obj,ids);
        // Get the operand DOFs associated with the given element
        std::vector<size_t> row = operandDOF->getRowDOFs(obj);
        size_t nnz = row.size();
        for (size_t i=0; i<row.size(); i++)
            columns[i] = (int) row[i];
        // Add the rows
        for (size_t i=0; i<ids.size(); i++) {
            int globalRowID = ids[i];
            int localRowID = globalRowID - resultDOF->beginDOF();
            params->addMapping( localRowID, globalRowID );
            params->setEntriesInRow( localRowID, nnz );
        }
        // Add the columns
        params->addColumns( nnz, columns );
        // Increment the iterator (pre-increment for speed)
        ++cur_elem;
    }

    // Get the communication lists for the vectors
    params->d_CommListLeft = resultVec->getCommunicationList();
    params->d_CommListRight = operandVec->getCommunicationList();
    params->d_DOFManagerLeft = resultDOF;
    params->d_DOFManagerRight = operandDOF;

    // Create the matrix
    boost::shared_ptr<AMP::LinearAlgebra::ManagedPetscMatrix>  newMatrix( new AMP::LinearAlgebra::ManagedPetscMatrix(params) );

    // Initialize the matrix
    cur_elem = resultDOF->getIterator();
    end_elem = cur_elem.end();
    double  values[1000];
    for (size_t i=0 ; i<1000; i++) { values[i] = 0.0; }
    while ( cur_elem != end_elem) {
        AMP::Mesh::MeshElement obj = *cur_elem;
        // Get the result DOFs associated with the given element
        std::vector<size_t> ids;
        resultDOF->getDOFs(obj,ids);
        // Get the operand DOFs associated with the given element
        std::vector<size_t> row = operandDOF->getRowDOFs(obj);
        size_t nnz = row.size();
        for (size_t i=0; i<row.size(); i++)
            columns[i] = (int) row[i];
        // Add the rows
        for (size_t i=0; i<ids.size(); i++) {
            int globalRowID = ids[i];
            newMatrix->createValuesByGlobalID( 1, nnz, &globalRowID, columns, values );
        }
        ++cur_elem;
    }
    newMatrix->castTo<AMP::LinearAlgebra::EpetraMatrix>().setEpetraMaps( resultVec, operandVec );
    newMatrix->makeConsistent();

    return newMatrix;
}


}
}

#endif

