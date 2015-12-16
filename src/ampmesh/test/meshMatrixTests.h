#ifndef included_MeshMatrixTests
#define included_MeshMatrixTests

#include "ampmesh/Mesh.h"
#include "matrices/Matrix.h"
#include "matrices/MatrixBuilder.h"
#include "vectors/VectorBuilder.h"

#include "meshTests.h"


template <int DOF_PER_NODE, bool SPLIT>
void VerifyGetMatrixTrivialTest( AMP::UnitTest *utils, AMP::Mesh::Mesh::shared_ptr mesh )
{
    // Create the DOF_Manager
    AMP::Discretization::DOFManagerParameters::shared_ptr DOFparams(
        new AMP::Discretization::DOFManagerParameters( mesh ) );
    AMP::Discretization::DOFManager::shared_ptr DOFs =
        AMP::Discretization::simpleDOFManager::create( mesh, AMP::Mesh::Vertex, 1, DOF_PER_NODE );

    // Create a nodal variable
    AMP::LinearAlgebra::Variable::shared_ptr variable(
        new AMP::LinearAlgebra::Variable( "test vector" ) );

    // Create the matrix and vectors
    AMP::LinearAlgebra::Vector::shared_ptr vector1 =
        AMP::LinearAlgebra::createVector( DOFs, variable, SPLIT );
    AMP::LinearAlgebra::Vector::shared_ptr vector2 =
        AMP::LinearAlgebra::createVector( DOFs, variable, SPLIT );
    AMP::LinearAlgebra::Matrix::shared_ptr matrixa =
        AMP::LinearAlgebra::createMatrix( vector1, vector2 );

    // Currently there is a bug with multivectors
    bool isMultiVector = vector1->isA<AMP::LinearAlgebra::MultiVector>();
    if ( isMultiVector ) {
        utils->expected_failure( "VerifyGetMatrixTrivialTest with split=true" );
        return;
    }

    // Run some tests
    vector1->setRandomValues();
    matrixa->makeConsistent();
    matrixa->mult( vector1, vector2 );
    if ( vector2->L1Norm() < 0.00000001 )
        utils->passes( "obtained 0 matrix from mesh" );
    else
        utils->failure( "did not obtain 0 matrix from mesh" );

    // Need to get another matrix to store data due to Epetra insert/replace idiom.  Matrixa is
    // fixed with no entires.
    AMP::LinearAlgebra::Matrix::shared_ptr matrixb =
        AMP::LinearAlgebra::createMatrix( vector1, vector2 );

    vector2->setToScalar( 1. );
    matrixb->makeConsistent();
    matrixb->setDiagonal( vector2 );
    matrixb->mult( vector1, vector2 );
    vector1->subtract( vector1, vector2 );

    if ( vector1->L1Norm() < 0.0000001 )
        utils->passes( "created identity matrix from mesh" );
    else
        utils->failure( "created identity matrix from mesh" );
}


template <int DOF_PER_NODE, bool SPLIT>
void GhostWriteTest( AMP::UnitTest *utils, AMP::Mesh::Mesh::shared_ptr mesh )
{
    // Create the DOF_Manager
    AMP::Discretization::DOFManagerParameters::shared_ptr DOFparams(
        new AMP::Discretization::DOFManagerParameters( mesh ) );
    AMP::Discretization::DOFManager::shared_ptr DOFs =
        AMP::Discretization::simpleDOFManager::create( mesh, AMP::Mesh::Vertex, 1, DOF_PER_NODE );

    // Create a nodal variable
    AMP::LinearAlgebra::Variable::shared_ptr variable(
        new AMP::LinearAlgebra::Variable( "test vector" ) );

    // Create the matrix and vectors
    AMP::LinearAlgebra::Vector::shared_ptr vector1 =
        AMP::LinearAlgebra::createVector( DOFs, variable, SPLIT );
    AMP::LinearAlgebra::Vector::shared_ptr vector2 =
        AMP::LinearAlgebra::createVector( DOFs, variable, SPLIT );
    AMP::LinearAlgebra::Matrix::shared_ptr matrix =
        AMP::LinearAlgebra::createMatrix( vector1, vector2 );

    // For each mesh, get a mapping of it's processor id's to the comm of the mesh
    std::map<AMP::Mesh::MeshID, std::vector<int>> proc_map = createRankMap( mesh );
    NULL_USE( proc_map );

    // For each processor, make sure it can write to all entries
    AMP::AMP_MPI comm = mesh->getComm();
    for ( int p = 0; p < comm.getSize(); p++ ) {
        matrix->setScalar( -1.0 );
        matrix->makeConsistent();
        // Processor p should fill the values
        if ( p == comm.getRank() ) {
            try {
                double proc = mesh->getComm().getRank();
                bool passes = true;
                // Loop through the owned nodes
                AMP::Mesh::MeshIterator iterator = mesh->getIterator( AMP::Mesh::Vertex, 0 );
                for ( size_t i = 0; i < iterator.size(); i++ ) {
                    // Get the DOFs for the node and it's neighbors
                    std::vector<size_t> localDOFs;
                    DOFs->getDOFs( iterator->globalID(), localDOFs );
                    std::vector<size_t> neighborDOFs, dofs;
                    std::vector<AMP::Mesh::MeshElement::shared_ptr> elements =
                        iterator->getNeighbors();
                    for ( size_t i = 0; i < elements.size(); i++ ) {
                        DOFs->getDOFs( elements[i]->globalID(), dofs );
                        for ( size_t j = 0; j < dofs.size(); j++ ) {
                            neighborDOFs.push_back( dofs[j] );
                        }
                    }
                    // For each local DOF, set all matrix elements involving the current DOF
                    for ( size_t j = 0; j < localDOFs.size(); j++ ) {
                        for ( size_t i = 0; i < localDOFs.size(); i++ ) {
                            matrix->setValueByGlobalID( localDOFs[i], localDOFs[j], proc );
                            matrix->setValueByGlobalID( localDOFs[j], localDOFs[i], proc );
                        }
                        for ( size_t i = 0; i < neighborDOFs.size(); i++ ) {
                            matrix->setValueByGlobalID( localDOFs[j], neighborDOFs[i], proc );
                            matrix->setValueByGlobalID( neighborDOFs[i], localDOFs[j], proc );
                        }
                        std::vector<unsigned int> cols;
                        std::vector<double> values;
                        matrix->getRowByGlobalID( localDOFs[j], cols, values );
                        for ( size_t i1 = 0; i1 < cols.size(); i1++ ) {
                            for ( size_t i2 = 0; i2 < localDOFs.size(); i2++ ) {
                                if ( cols[i1] == localDOFs[i2] ) {
                                    if ( values[i1] != proc )
                                        passes = false;
                                }
                            }
                            for ( size_t i2 = 0; i2 < neighborDOFs.size(); i2++ ) {
                                if ( cols[i1] == neighborDOFs[i2] ) {
                                    if ( values[i1] != proc )
                                        passes = false;
                                }
                            }
                        }
                    }
                    ++iterator;
                }
                if ( passes )
                    utils->passes( "Able to write to ghost entries in matrix" );
                else
                    utils->failure( "Able to write to ghost entries in matrix" );
            } catch ( ... ) {
                utils->failure( "Able to write to ghost entries in matrix (exception)" );
            }
        }

        // Apply make consistent
        matrix->makeConsistent();

        /*// Test that each matrix entry has the proper value
        bool passes = true;
        // Get a list of all nodes owned by the given processor p
        std::set<AMP::Mesh::MeshElementID> nodes_p;
        AMP::Mesh::MeshIterator iterator = mesh->getIterator(AMP::Mesh::Vertex,1);
        for (size_t i=0; i<iterator.size(); i++) {
            AMP::Mesh::MeshElementID id = iterator->globalID();
            const std::vector<int> &map = proc_map.find(id.meshID())->second;
            int rank = map[id.owner_rank()];
            if ( rank == p )
                nodes_p.insert( iterator->globalID() );
            ++iterator;
        }
        // Get a list of all DOFs associated with the given nodes
        std::vector<size_t> dofs_p;
        std::vector<AMP::Mesh::MeshElementID> tmp(nodes_p.begin(),nodes_p.end());
        DOFs->getDOFs( tmp, dofs_p );
        if ( dofs_p.size()==0 )
            dofs_p.push_back( (size_t)-1 );
        AMP::Utilities::quicksort( dofs_p );
        iterator = mesh->getIterator(AMP::Mesh::Vertex,0);
        double proc = p;
        for (size_t i=0; i<iterator.size(); i++) {
            std::vector<size_t> dofs;
            DOFs->getDOFs( iterator->globalID(), dofs );
            for (size_t j=0; j<dofs.size(); j++) {
                size_t row = dofs[j];
                size_t index = AMP::Utilities::findfirst( dofs_p, row );
                if ( index==dofs_p.size() ) { index--; }
                bool row_found = row==dofs_p[index];
                std::vector<unsigned int> cols;
                std::vector<double> values;
                matrix->getRowByGlobalID( row, cols, values );
                for (size_t k=0; k<cols.size(); k++) {
                    index = AMP::Utilities::findfirst( dofs_p, (size_t) cols[k] );
                    if ( index==dofs_p.size() ) { index--; }
                    bool col_found = cols[k]==dofs_p[index];
                    if ( row_found || col_found ) {
                        if ( values[k] != proc )
                            passes = false;
                    } else {
                        if ( values[k] != -1.0 )
                            passes = false;
                    }
                }
            }
            ++iterator;
        }
        char msg[100];
        sprintf(msg,"Matrix entries set by processor %i read correctly on processor
        %i",p,comm.getRank());
        if ( passes )
            utils->passes( msg );
        else
            utils->failure( msg );*/
    }
}


#endif
