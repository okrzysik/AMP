#ifdef USE_AMP_VECTORS

#include "matrices/MatrixBuilder.h"
#include "matrices/DenseSerialMatrix.h"

#include "discretization/DOF_Manager.h"

#ifdef USE_EXT_PETSC
    #include "vectors/petsc/ManagedPetscVector.h"
    #include "matrices/petsc/ManagedPetscMatrix.h"
#endif
#ifdef USE_EXT_TRILINOS
    #include "vectors/trilinos/EpetraVectorEngine.h"
#endif


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
    if ( comm.getSize()==1 )
        comm = AMP_MPI(AMP_COMM_SELF);

#if defined(USE_EXT_PETSC) && defined(USE_EXT_TRILINOS)
    // Create the matrix parameters
    boost::shared_ptr<AMP::LinearAlgebra::ManagedPetscMatrixParameters> params( 
        new AMP::LinearAlgebra::ManagedPetscMatrixParameters ( resultDOF, operandDOF, comm ) );
    params->d_CommListLeft = resultVec->getCommunicationList();
    params->d_CommListRight = operandVec->getCommunicationList();
    params->d_VariableLeft = resultVec->getVariable();
    params->d_VariableRight = operandVec->getVariable();

    // Add the rows to the matrix parameters
    AMP::Mesh::MeshIterator cur_elem = resultDOF->getIterator();
    AMP::Mesh::MeshIterator end_elem = cur_elem.end();
    int columns[1000];   // Preallocate for the columns for speed
    for (size_t i=0 ; i<1000; i++) { columns[i] = 0.0; }
    while ( cur_elem != end_elem) {
        AMP::Mesh::MeshElement obj = *cur_elem;
        // Get the result DOFs associated with the given element
        std::vector<size_t> ids;
        resultDOF->getDOFs(obj.globalID(),ids);
        // Get the operand DOFs associated with the given element
        std::vector<size_t> row = operandDOF->getRowDOFs(obj);
        size_t nnz = row.size();
        for (size_t i=0; i<row.size(); i++)
            columns[i] = (int) row[i];
        // Add the rows
        for (size_t i=0; i<ids.size(); i++) {
            int globalRowID = ids[i];
            int localRowID = globalRowID - resultDOF->beginDOF();
            params->setEntriesInRow( localRowID, nnz );
        }
        // Add the columns
        params->addColumns( nnz, columns );
        // Increment the iterator (pre-increment for speed)
        ++cur_elem;
    }

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
        resultDOF->getDOFs(obj.globalID(),ids);
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
    newMatrix->zero();
    newMatrix->makeConsistent();

    return newMatrix;
#else
    if ( comm.getSize()!=1 )
        AMP_ERROR("The only native matrix is a serial dense matrix");
    // Create the matrix parameters
    boost::shared_ptr<AMP::LinearAlgebra::MatrixParameters> params( 
        new AMP::LinearAlgebra::MatrixParameters( resultDOF, operandDOF, comm ) );
    // Create the matrix
    boost::shared_ptr<AMP::LinearAlgebra::DenseSerialMatrix>  newMatrix( new AMP::LinearAlgebra::DenseSerialMatrix(params) );
    // Initialize the matrix
    newMatrix->zero();
    newMatrix->makeConsistent();
    return newMatrix;
#endif
}


}
}

#endif

