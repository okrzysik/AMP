#ifndef included_MeshMatrixTests
#define included_MeshMatrixTests


template <int DOF_PER_NODE>
class VerifyGetMatrixTrivialTest
{
public:
    static const char * get_test_name () { return "Verify matrix interface trivially in MeshAdapter"; }

    static  void run_test ( AMP::UnitTest *utils, AMP::Mesh::MeshAdapter::shared_ptr mesh ) {

          AMP::LinearAlgebra::Variable::shared_ptr variable ( new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable,DOF_PER_NODE> ( "test vector" ) );
          AMP::LinearAlgebra::Matrix::shared_ptr matrixa = mesh->createMatrix ( variable );  // Generates new mtrix
          AMP::LinearAlgebra::Vector::shared_ptr vectorb = mesh->createVector ( variable );  // Gets vector from the cached copy
          AMP::LinearAlgebra::Vector::shared_ptr vectorc = mesh->createVector ( variable );  // Gets vector from the cached copy

          vectorb->setRandomValues ();
          matrixa->makeConsistent ();
          matrixa->mult ( vectorb , vectorc );
          if ( vectorc->L1Norm() < 0.00000001 )
            utils->passes ( "obtained 0 matrix from mesh" );
          else
            utils->failure ( "did not obtain 0 matrix from mesh" );

          AMP::LinearAlgebra::Matrix::shared_ptr matrixb = mesh->createMatrix ( variable );   // Need to get another matrix to
                                                                           // store data due to Epetra insert/replace
                                                                           // idiom.  Matrixa is fixed with no entires.
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
};


#endif

