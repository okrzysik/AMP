#ifndef included_test_MatrixTests
#define included_test_MatrixTests

#include "utils/UnitTest.h"
#include "discretization/DOF_Manager.h"

#include "../../vectors/test/test_VectorLoops.h"

#ifdef USE_PETSC
    #include "matrices/petsc/PetscMatrix.h"
#endif


namespace AMP {
namespace unit_test {

#if defined(USE_PETSC) && defined(USE_PETSC)

template <class MATRIX_FACTORY>
void  fillWithPseudoLaplacian ( AMP::LinearAlgebra::Matrix::shared_ptr matrix )
{
    AMP::Discretization::DOFManager::shared_ptr  dofmap = MATRIX_FACTORY::getDOFMap();
    for ( size_t i = dofmap->beginDOF() ; i != dofmap->endDOF() ; i++ )
    {
        std::vector<unsigned int>  cols;
        std::vector<double>        vals;
        matrix->getRowByGlobalID ( i , cols , vals );
        for ( size_t j = 0 ; j != cols.size() ; j++ )
        {
            if ( cols[j] == i )
                vals[j] = 6;
            else
                vals[j] = -1;
        }
        if ( cols.size() )
            matrix->setValuesByGlobalID ( 1 , (int)cols.size() , (int *)&i , (int *)&(cols[0]) , (double *)&(vals[0]) );
    }
    matrix->makeConsistent ();
}


template <typename FACTORY>
class InstantiateMatrix
{
    public:
      static const char * get_test_name () { return "instantiate matrix"; }

      static  void run_test ( AMP::UnitTest *utils )
      {
        AMP::LinearAlgebra::Matrix::shared_ptr  matrix = FACTORY::getMatrix();
        if ( matrix )
            utils->passes ( "created" );
        else
            utils->failure ( "created" );
      }
};


template <typename MATRIX_FACTORY>
class AmpInterfaceLeftVectorFactory
{
    public:
      typedef  AMP::LinearAlgebra::ManagedPetscVector        vector;

      static AMP::LinearAlgebra::Variable::shared_ptr   getVariable()
      {
        return AMP::LinearAlgebra::Variable::shared_ptr( new AMP::LinearAlgebra::Variable("left") );
      }

      static AMP::LinearAlgebra::Vector::shared_ptr     getVector()
      {
        AMP::LinearAlgebra::Matrix::shared_ptr  matrix = MATRIX_FACTORY::getMatrix();
        AMP::LinearAlgebra::Vector::shared_ptr  vector = matrix->getLeftVector();
        vector->setVariable( getVariable() );
        return vector;
      }
};


template <typename MATRIX_FACTORY>
class AmpInterfaceRightVectorFactory
{
    public:
      typedef  AMP::LinearAlgebra::ManagedPetscVector        vector;

      static AMP::LinearAlgebra::Variable::shared_ptr   getVariable()
      {
        return AMP::LinearAlgebra::Variable::shared_ptr( new AMP::LinearAlgebra::Variable("right") );
      }

      static AMP::LinearAlgebra::Vector::shared_ptr     getVector()
      {
        AMP::LinearAlgebra::Matrix::shared_ptr  matrix = MATRIX_FACTORY::getMatrix ();
        AMP::LinearAlgebra::Vector::shared_ptr  vector = matrix->getRightVector();
        vector->setVariable( getVariable() );
        return vector;

      }
};


template <typename MATRIX_FACTORY>
class PETScInterfaceLeftVectorFactory
{
    public:
      typedef  AMP::LinearAlgebra::ManagedPetscVector        vector;

      static AMP::LinearAlgebra::Variable::shared_ptr   getVariable()
      {
        return AMP::LinearAlgebra::Variable::shared_ptr( new AMP::LinearAlgebra::Variable("petsc_left") );
      }

      static AMP::LinearAlgebra::Vector::shared_ptr     getVector()
      {
        AMP::LinearAlgebra::Matrix::shared_ptr  matrix = MATRIX_FACTORY::getMatrix();
        ::Mat m = matrix->castTo<AMP::LinearAlgebra::PetscMatrix>().getMat ();
        ::Vec v;
        MatGetVecs ( m , &v , 0 );

        boost::shared_ptr<AMP::LinearAlgebra::NativePetscVectorParameters> p( 
            new AMP::LinearAlgebra::NativePetscVectorParameters( v, true ) );
        AMP::LinearAlgebra::Vector::shared_ptr  vector( new AMP::LinearAlgebra::NativePetscVector ( p ) );
        vector->setVariable( getVariable() );
        return vector;
      }

      static AMP::LinearAlgebra::Vector::shared_ptr   getNativeVector()
      {
        return MATRIX_FACTORY::getMatrix()->getLeftVector();
      }

      static AMP::LinearAlgebra::Vector::shared_ptr   getManagedVector()
      {
        AMP::LinearAlgebra::Matrix::shared_ptr  matrix = MATRIX_FACTORY::getMatrix();
        return AMP::LinearAlgebra::PetscVector::view ( matrix->getLeftVector () );
      }
};


template <typename MATRIX_FACTORY>
class PETScInterfaceRightVectorFactory
{
    public:
      typedef  AMP::LinearAlgebra::ManagedPetscVector        vector;

      static AMP::LinearAlgebra::Variable::shared_ptr   getVariable()
      {
        return AMP::LinearAlgebra::Variable::shared_ptr( new AMP::LinearAlgebra::Variable("petsc_right") );
      }

      static AMP::LinearAlgebra::Vector::shared_ptr     getVector()
      {
        AMP::LinearAlgebra::Matrix::shared_ptr  matrix = MATRIX_FACTORY::getMatrix();
        ::Mat m = matrix->castTo<AMP::LinearAlgebra::PetscMatrix>().getMat ();
        ::Vec v;
        MatGetVecs ( m , &v , 0 );

        boost::shared_ptr<AMP::LinearAlgebra::NativePetscVectorParameters> p( 
            new AMP::LinearAlgebra::NativePetscVectorParameters( v, true ) );
        AMP::LinearAlgebra::Vector::shared_ptr  vector( new AMP::LinearAlgebra::NativePetscVector ( p ) );
        vector->setVariable( getVariable() );
        return vector;
      }

      static AMP::LinearAlgebra::Vector::shared_ptr   getNativeVector()
      {
        return getVector();
      }

      static AMP::LinearAlgebra::Vector::shared_ptr   getManagedVector()
      {
        AMP::LinearAlgebra::Matrix::shared_ptr  matrix = MATRIX_FACTORY::getMatrix();
        return AMP::LinearAlgebra::PetscVector::view ( matrix->getRightVector () );
      }
};


template <typename FACTORY>
class VerifyGetLeftRightVector
{
    public:
      static const char * get_test_name () { return "verify getRightVector"; }

      static  void run_test ( AMP::UnitTest *utils )
      {
        test_managed_vectors_loop <AmpInterfaceRightVectorFactory <FACTORY> > ( utils );
        test_managed_vectors_loop <PETScInterfaceRightVectorFactory <FACTORY> > ( utils );
        // test_managed_vectors_loop <AmpInterfaceLeftVectorFactory <FACTORY> > ( utils );
        // test_managed_vectors_loop <PETScInterfaceLeftVectorFactory <FACTORY> > ( utils );
      }

};


template <typename FACTORY>
class VerifyGetSetValuesMatrix
{
    public:
      static const char * get_test_name () { return "verify get and set"; }

      static  void run_test ( AMP::UnitTest *utils )
      {
        AMP::LinearAlgebra::Matrix::shared_ptr  matrix = FACTORY::getMatrix();
        AMP::Discretization::DOFManager::shared_ptr  dofmap = FACTORY::getDOFMap();

        matrix->makeConsistent();
        fillWithPseudoLaplacian<FACTORY>( matrix ); // puts 6 on the diagonal and -1 on allocated off-diagonals
        for ( size_t i = dofmap->beginDOF() ; i != dofmap->endDOF() ; i++ )
        {
          std::vector<unsigned int>  cols;
          std::vector<double>        vals;
          for ( size_t j = 0 ; j != cols.size() ; j++ )
          {
            double ans = ( i == cols[j] ) ? 6. : -1.;
            if ( vals[j] != ans )
            {
              utils->failure ( "bad value in matrix" );
              return;
            }
          }
        }
        utils->passes ("verify get and set");
      }
};


template <typename FACTORY>
class VerifyAXPYMatrix
{
    public:
      static const char * get_test_name () { return "verify AXPY"; }

      static  void run_test ( AMP::UnitTest *utils )
      {
        AMP::LinearAlgebra::Matrix::shared_ptr  matrix1 = FACTORY::getMatrix();
        AMP::LinearAlgebra::Matrix::shared_ptr  matrix2 = FACTORY::getMatrix();
        AMP::LinearAlgebra::Vector::shared_ptr  vector1lhs = matrix1->getRightVector();
        AMP::LinearAlgebra::Vector::shared_ptr  vector2lhs = matrix2->getRightVector();
        AMP::LinearAlgebra::Vector::shared_ptr  vector1rhs = matrix1->getRightVector();
        AMP::LinearAlgebra::Vector::shared_ptr  vector2rhs = matrix2->getRightVector();
        AMP::LinearAlgebra::Vector::shared_ptr  vectorresult = matrix2->getRightVector();

        fillWithPseudoLaplacian<FACTORY>( matrix1 );
        fillWithPseudoLaplacian<FACTORY>( matrix2 );

        matrix2->axpy ( -2. , matrix1 );   // matrix2 = -matrix1
        vector1lhs->setRandomValues ();
        vector2lhs->copyVector ( vector1lhs );

        matrix1->mult ( vector1lhs , vector1rhs );
        matrix2->mult ( vector2lhs , vector2rhs );  // vector2rhs = - vector1rhs

        vectorresult->add ( vector1rhs , vector2rhs );
        if ( vectorresult->L1Norm() < 0.000001 )
          utils->passes ( "matrices are opposite" );
        else
          utils->failure ( "matrices are not opposite" );
        if ( vector1rhs->L1Norm() > 0.00001 )
          utils->passes ( "non-trivial vector" );
        else
          utils->passes ( "trivial vector" );
      }
};


template <typename FACTORY>
class VerifyScaleMatrix
{
    public:
      static const char * get_test_name () { return "verify scale"; }

      static  void run_test ( AMP::UnitTest *utils )
      {
        AMP::LinearAlgebra::Matrix::shared_ptr  matrix1 = FACTORY::getMatrix();
        AMP::LinearAlgebra::Matrix::shared_ptr  matrix2 = FACTORY::getMatrix();
        AMP::LinearAlgebra::Vector::shared_ptr  vector1lhs = matrix1->getRightVector();
        AMP::LinearAlgebra::Vector::shared_ptr  vector2lhs = matrix2->getRightVector();
        AMP::LinearAlgebra::Vector::shared_ptr  vector1rhs = matrix1->getRightVector();
        AMP::LinearAlgebra::Vector::shared_ptr  vector2rhs = matrix2->getRightVector();
        AMP::LinearAlgebra::Vector::shared_ptr  vectorresult = matrix2->getRightVector();

        fillWithPseudoLaplacian<FACTORY>( matrix1 );
        fillWithPseudoLaplacian<FACTORY>( matrix2 );

        matrix2->scale ( 1.234567 );            // matrix2 = matrix1

        vector1lhs->setRandomValues ();
        vector2lhs->copyVector ( vector1lhs );

        matrix1->mult ( vector1lhs , vector1rhs );
        matrix2->mult ( vector2lhs , vector2rhs );  // vector2rhs = 1.234567vector1rhs

        vector1rhs->scale ( 1.234567 );
        vectorresult->subtract ( vector1rhs , vector2rhs );

        if ( vectorresult->L1Norm() < 0.000001 )
          utils->passes ( "matrices are equally scaled" );
        else
          utils->failure ( "matrices are equally scaled" );
        if ( vector1rhs->L1Norm() > 0.00001 )
          utils->passes ( "non-trivial vector" );
        else
          utils->passes ( "trivial vector" );
      }
};


template <typename FACTORY>
class VerifyMultMatrix
{
    public:
      static const char * get_test_name () { return "verify mult"; }

      static  void run_test ( AMP::UnitTest *utils )
      {
        AMP::LinearAlgebra::Matrix::shared_ptr  matrix = FACTORY::getMatrix();
        AMP::LinearAlgebra::Vector::shared_ptr  vectorlhs = matrix->getRightVector();
        AMP::LinearAlgebra::Vector::shared_ptr  vectorrhs = matrix->getRightVector();
        double normlhs , normrhs;

        /* Verify 0 matrix from factory */
        vectorlhs->setRandomValues ();
        matrix->mult ( vectorlhs , vectorrhs );
        normlhs = vectorlhs->L2Norm();
        normrhs = vectorrhs->L2Norm();
        if ( (normlhs > 0) && (normrhs < 0.0000001) )
          utils->passes ( "mult by 0 matrix" );
        else
          utils->failure ( "mult by 0 matrix" );

        /* Verify mult with identity */
        vectorlhs->setToScalar ( 1.0 );
        matrix->setDiagonal ( vectorlhs );
        
        vectorlhs->setRandomValues ();
        matrix->mult ( vectorlhs , vectorrhs );
        normlhs = vectorlhs->L2Norm();
        vectorrhs->subtract ( vectorlhs , vectorrhs );
        normrhs = vectorrhs->L2Norm();
        if ( (normlhs > 0) && (normrhs < 0.0000001) )
          utils->passes ( "mult by I matrix" );
        else
          utils->failure ( "mult by I matrix" );

        /* Try the non-trivial matrix */
        fillWithPseudoLaplacian<FACTORY>( matrix );
        vectorlhs->setRandomValues ();
        matrix->mult ( vectorlhs , vectorrhs );
        normlhs = vectorlhs->L2Norm();
        normrhs = vectorrhs->L2Norm();
        if ( (normlhs > 0) && (normrhs > 0) )
          utils->passes ( "mult by other matrix" );
        else
          utils->failure ( "mult by other matrix" );

      }
};


template <typename FACTORY>
class VerifyExtractDiagonal
{
    public:
      static const char * get_test_name () { return "Verify extractDiagonal"; }

      static  void run_test ( AMP::UnitTest *utils )
      {
        AMP::LinearAlgebra::Matrix::shared_ptr   matrix = FACTORY::getMatrix();
        AMP::LinearAlgebra::Vector::shared_ptr   vector = matrix->getRightVector();
        size_t  firstRow = vector->getCommunicationList()->getStartGID();
        size_t  maxCols = matrix->numGlobalColumns();
        for ( size_t i = 0 ; i != vector->getCommunicationList()->numLocalRows() ; i++ )
        {
          int row = static_cast<int> ( i + firstRow );
          if ( row >= static_cast<int> ( maxCols ) )
          {
            break;
          }
          matrix->setValueByGlobalID ( row , row , static_cast<double> ( row+1 ) );
        }
        AMP::LinearAlgebra::Vector::shared_ptr   diag = matrix->extractDiagonal ();
        double l1norm = diag->L1Norm();
        double numRows = static_cast<double> ( matrix->numGlobalRows() );
        double numCols = static_cast<double> ( matrix->numGlobalColumns() );
        double compareVal = std::min ( numRows , numCols );
        compareVal = compareVal * (1. + compareVal)/2.;
        if ( fabs ( l1norm - compareVal ) < 0.000001 )
          utils->passes ("Verify extractDiagonal");
        else
          utils->failure ("Verify extractDiagonal");
      }
};

#endif

}
}


#endif
