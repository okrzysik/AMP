#ifndef included_MeshMatrixTests
#define included_MeshMatrixTests

#include "ampmesh/Mesh.h"
#include "matrices/Matrix.h"


template <int DOF_PER_NODE>
void VerifyGetMatrixTrivialTest( AMP::UnitTest *utils, AMP::Mesh::Mesh::shared_ptr mesh ) 
{
    // Create the DOF_Manager
    AMP::Discretization::DOFManagerParameters::shared_ptr DOFparams( new AMP::Discretization::DOFManagerParameters(mesh) );
    boost::shared_ptr<AMP::Discretization::simpleDOFManager> DOFs( new AMP::Discretization::simpleDOFManager(mesh,AMP::Mesh::Vertex,1,DOF_PER_NODE) );

    // Create a nodal variable 
    AMP::LinearAlgebra::Variable::shared_ptr variable( new AMP::Discretization::NodalVariable(DOF_PER_NODE,"test vector") );

    // Create the matrix and vectors
    AMP::LinearAlgebra::Matrix::shared_ptr matrixa = DOFs->createMatrix ( variable, variable );
    AMP::LinearAlgebra::Vector::shared_ptr vectorb = DOFs->createVector ( variable );
    AMP::LinearAlgebra::Vector::shared_ptr vectorc = DOFs->createVector ( variable );

    // Run some tests
    vectorb->setRandomValues ();
    matrixa->makeConsistent ();
    matrixa->mult ( vectorb , vectorc );
    if ( vectorc->L1Norm() < 0.00000001 )
        utils->passes ( "obtained 0 matrix from mesh" );
    else
        utils->failure ( "did not obtain 0 matrix from mesh" );

    // Need to get another matrix to store data due to Epetra insert/replace idiom.  Matrixa is fixed with no entires.
    AMP::LinearAlgebra::Matrix::shared_ptr matrixb = DOFs->createMatrix ( variable, variable );

    vectorc->setToScalar ( 1. );
    matrixb->makeConsistent ();
    matrixb->setDiagonal ( vectorc );
    matrixb->mult ( vectorb , vectorc );
    vectorb->subtract ( vectorb , vectorc );

    if ( vectorb->L1Norm() < 0.0000001 )
        utils->passes ( "created identity matrix from mesh" );
    else
        utils->failure ( "created identity matrix from mesh" );
}


#endif

