#ifndef included_MeshMatrixTests
#define included_MeshMatrixTests

#include "ampmesh/Mesh.h"
#include "matrices/Matrix.h"
#include "vectors/VectorBuilder.h"
#include "matrices/MatrixBuilder.h"

template <int DOF_PER_NODE, bool SPLIT>
void VerifyGetMatrixTrivialTest( AMP::UnitTest *utils, AMP::Mesh::Mesh::shared_ptr mesh ) 
{
    // Create the DOF_Manager
    AMP::Discretization::DOFManagerParameters::shared_ptr DOFparams( new AMP::Discretization::DOFManagerParameters(mesh) );
    AMP::Discretization::DOFManager::shared_ptr DOFs = AMP::Discretization::simpleDOFManager::create(mesh,AMP::Mesh::Vertex,1,DOF_PER_NODE);

    // Create a nodal variable 
    AMP::LinearAlgebra::Variable::shared_ptr variable( new AMP::LinearAlgebra::Variable("test vector") );

    // Create the matrix and vectors
    AMP::LinearAlgebra::Vector::shared_ptr vector1 = AMP::LinearAlgebra::createVector ( DOFs, variable, SPLIT );
    AMP::LinearAlgebra::Vector::shared_ptr vector2 = AMP::LinearAlgebra::createVector ( DOFs, variable, SPLIT );
    AMP::LinearAlgebra::Matrix::shared_ptr matrixa = AMP::LinearAlgebra::createMatrix ( vector1, vector2 );

    // Run some tests
    vector1->setRandomValues ();
    matrixa->makeConsistent ();
    bool isMultiVector = vector1->isA<AMP::LinearAlgebra::MultiVector>();
    matrixa->mult ( vector1 , vector2 );
    if ( vector2->L1Norm() < 0.00000001 )
        utils->passes ( "obtained 0 matrix from mesh" );
    else
        utils->failure ( "did not obtain 0 matrix from mesh" );

    // Need to get another matrix to store data due to Epetra insert/replace idiom.  Matrixa is fixed with no entires.
    AMP::LinearAlgebra::Matrix::shared_ptr matrixb = AMP::LinearAlgebra::createMatrix ( vector1, vector2 );

    vector2->setToScalar ( 1. );
    matrixb->makeConsistent ();
    matrixb->setDiagonal ( vector2 );
    matrixb->mult ( vector1 , vector2 );
    vector1->subtract ( vector1 , vector2 );

    if ( vector1->L1Norm() < 0.0000001 )
        utils->passes ( "created identity matrix from mesh" );
    else
        utils->failure ( "created identity matrix from mesh" );
}


template <int DOF_PER_NODE, bool SPLIT>
void GhostWriteTest( AMP::UnitTest *utils, AMP::Mesh::Mesh::shared_ptr mesh ) 
{
    // Create the DOF_Manager
    AMP::Discretization::DOFManagerParameters::shared_ptr DOFparams( new AMP::Discretization::DOFManagerParameters(mesh) );
    AMP::Discretization::DOFManager::shared_ptr DOFs = AMP::Discretization::simpleDOFManager::create(mesh,AMP::Mesh::Vertex,1,DOF_PER_NODE);

    // Create a nodal variable 
    AMP::LinearAlgebra::Variable::shared_ptr variable( new AMP::LinearAlgebra::Variable("test vector") );

    // Create the matrix and vectors
    AMP::LinearAlgebra::Vector::shared_ptr vector1 = AMP::LinearAlgebra::createVector ( DOFs, variable, SPLIT );
    AMP::LinearAlgebra::Vector::shared_ptr vector2 = AMP::LinearAlgebra::createVector ( DOFs, variable, SPLIT );
    AMP::LinearAlgebra::Matrix::shared_ptr matrix = AMP::LinearAlgebra::createMatrix ( vector1, vector2 );

    try {
        // Loop through the owned and ghost nodes
        AMP::Mesh::MeshIterator el = mesh->getIterator(AMP::Mesh::Volume,0);
        AMP::Mesh::MeshIterator end_el = el.end();
        for( ; el != end_el; ++el) {
            // Get the DOFs for all nodes
            std::vector<size_t> dofs;
            std::vector<size_t> dofIndices;
            std::vector<AMP::Mesh::MeshElement> elements = el->getElements(AMP::Mesh::Vertex);
            for (size_t i=0; i<elements.size(); i++) {
                DOFs->getDOFs( elements[i].globalID(), dofs );
                for (size_t j=0; j<dofs.size(); j++) {
                    dofIndices.push_back(dofs[j]);
                }
            }
            for (size_t j=0; j<dofIndices.size(); j++) {
                for (size_t i=0; i<dofIndices.size(); i++) {
                    matrix->setValueByGlobalID ( dofIndices[j], dofIndices[i], 1.0 );
                }//end for i
            }//end for j
        }//end for el

        // Apply make consistent
        matrix->makeConsistent();

        utils->passes ( "Able to write to ghost entries in matrix" );
    } catch (...) {
        utils->failure ( "Able to write to ghost entries in matrix" );
    }
}


#endif

