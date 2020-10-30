#ifndef included_AMP_test_MatrixTests
#define included_AMP_test_MatrixTests

#include "AMP/discretization/DOF_Manager.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/testHelpers/VectorTests.h"

#include "test_MatrixVectorFactory.h"

#if defined( USE_EXT_PETSC ) && defined( USE_EXT_TRILINOS )
#include "AMP/matrices/petsc/ManagedPetscMatrix.h"
#endif


namespace AMP {
namespace LinearAlgebra {


template<class MATRIX_FACTORY>
void fillWithPseudoLaplacian( AMP::LinearAlgebra::Matrix::shared_ptr matrix )
{
    auto dofmap = MATRIX_FACTORY::getDOFMap();
    for ( size_t i = dofmap->beginDOF(); i != dofmap->endDOF(); i++ ) {
        std::vector<size_t> cols;
        std::vector<double> vals;
        matrix->getRowByGlobalID( i, cols, vals );
        for ( size_t j = 0; j != cols.size(); j++ ) {
            if ( cols[j] == i )
                vals[j] = 6;
            else
                vals[j] = -1;
        }
        if ( cols.size() )
            matrix->setValuesByGlobalID(
                1, cols.size(), &i, &( cols[0] ), (double *) &( vals[0] ) );
    }
    matrix->makeConsistent();
}


template<typename FACTORY>
class InstantiateMatrix
{
public:
    static const char *get_test_name() { return "instantiate matrix"; }

    static void run_test( AMP::UnitTest *utils )
    {
        PROFILE_START( "InstantiateMatrix" );
        auto matrix = FACTORY::getMatrix();
        if ( matrix )
            utils->passes( "created" );
        else
            utils->failure( "created" );
        PROFILE_STOP( "InstantiateMatrix" );
    }
};


template<typename FACTORY>
class VerifyGetLeftRightVector
{
public:
    static const char *get_test_name() { return "verify getRightVector"; }

    static void run_test( AMP::UnitTest *utils )
    {
        PROFILE_START( "VerifyGetLeftRightVector" );
        global_cached_matrix = FACTORY::getMatrix();
        auto factory1        = std::make_shared<AmpInterfaceRightVectorFactory>();
        auto factory2        = std::make_shared<AmpInterfaceLeftVectorFactory>();
        VectorTests tests1( factory1 );
        VectorTests tests2( factory2 );
        tests1.testBasicVector( utils );
        tests2.testBasicVector( utils );
#if defined( USE_EXT_PETSC ) && defined( USE_EXT_TRILINOS )
        if ( std::dynamic_pointer_cast<AMP::LinearAlgebra::ManagedPetscMatrix>(
                 global_cached_matrix ) ) {
            auto factory3 = std::make_shared<PETScInterfaceRightVectorFactory>();
            auto factory4 = std::make_shared<PETScInterfaceLeftVectorFactory>();
            VectorTests tests3( factory3 );
            VectorTests tests4( factory4 );
            tests3.testPetsc( utils );
            tests4.testPetsc( utils );
        } else {
            utils->expected_failure(
                "PetscMatrix::createView is not ready for arbitrary matricies" );
        }
#endif
        global_cached_matrix.reset();
        PROFILE_STOP( "VerifyGetLeftRightVector" );
    }
};


template<typename FACTORY>
class VerifyGetSetValuesMatrix
{
public:
    static const char *get_test_name() { return "verify get and set"; }

    static void run_test( AMP::UnitTest *utils )
    {
        PROFILE_START( "VerifyGetSetValuesMatrix" );
        auto matrix = FACTORY::getMatrix();
        auto dofmap = FACTORY::getDOFMap();

        matrix->makeConsistent();
        fillWithPseudoLaplacian<FACTORY>(
            matrix ); // puts 6 on the diagonal and -1 on allocated off-diagonals
        for ( size_t i = dofmap->beginDOF(); i != dofmap->endDOF(); i++ ) {
            std::vector<size_t> cols;
            std::vector<double> vals;
            matrix->getRowByGlobalID( i, cols, vals );
            for ( size_t j = 0; j != cols.size(); j++ ) {
                double ans   = ( i == cols[j] ) ? 6. : -1.;
                double value = matrix->getValueByGlobalID( i, cols[j] );
                if ( vals[j] != ans || value != vals[j] ) {
                    utils->failure( "bad value in matrix" );
                    return;
                }
            }
        }
        utils->passes( "verify get and set" );
        PROFILE_STOP( "VerifyGetSetValuesMatrix" );
    }
};


template<typename FACTORY>
class VerifyAXPYMatrix
{
public:
    static const char *get_test_name() { return "verify AXPY"; }

    static void run_test( AMP::UnitTest *utils )
    {
        PROFILE_START( "VerifyAXPYMatrix" );

        // Create vectors/matricies from the factory
        auto matrix1      = FACTORY::getMatrix();
        auto matrix2      = FACTORY::getMatrix();
        auto vector1lhs   = matrix1->getRightVector();
        auto vector2lhs   = matrix2->getRightVector();
        auto vector1rhs   = matrix1->getRightVector();
        auto vector2rhs   = matrix2->getRightVector();
        auto vectorresult = matrix2->getRightVector();
        fillWithPseudoLaplacian<FACTORY>( matrix1 );
        fillWithPseudoLaplacian<FACTORY>( matrix2 );

        // Test axpy
        matrix2->axpy( -2., matrix1 ); // matrix2 = -matrix1
        vector1lhs->setRandomValues();
        vector2lhs->copyVector( vector1lhs );
        matrix1->mult( vector1lhs, vector1rhs );
        matrix2->mult( vector2lhs, vector2rhs ); // vector2rhs = - vector1rhs
        vectorresult->add( *vector1rhs, *vector2rhs );
        if ( vectorresult->L1Norm() < 0.000001 )
            utils->passes( "matrices are opposite" );
        else
            utils->failure( "matrices are not opposite" );
        if ( vector1rhs->L1Norm() > 0.00001 )
            utils->passes( "non-trivial vector" );
        else
            utils->passes( "trivial vector" );

        // Test that axpy failes with different sized matricies
        std::vector<size_t> row( 7 );
        for ( size_t i = 0; i < row.size(); i++ )
            row[i] = i;
        auto smallVec =
            AMP::LinearAlgebra::createSimpleVector<double>( 7, vector1lhs->getVariable() );
        auto smallMat = AMP::LinearAlgebra::createMatrix(
            smallVec, smallVec, FACTORY::type(), [row]( size_t ) { return row; } );
        try {
            matrix2->axpy( -2., smallMat ); // matrix2 = -matrix1
            utils->failure( "axpy did not crash with different sized matrices" );
        } catch ( std::exception &err ) {
            utils->passes( "axpy correctly fails with different sized matrices" );
        } catch ( ... ) {
            utils->failure( "axpy fails with different sized matrices (unknown failure)" );
        }
        PROFILE_STOP( "VerifyAXPYMatrix" );
    }
};


template<typename FACTORY>
class VerifyScaleMatrix
{
public:
    static const char *get_test_name() { return "verify scale"; }

    static void run_test( AMP::UnitTest *utils )
    {
        PROFILE_START( "VerifyScaleMatrix" );
        auto matrix1      = FACTORY::getMatrix();
        auto matrix2      = FACTORY::getMatrix();
        auto vector1lhs   = matrix1->getRightVector();
        auto vector2lhs   = matrix2->getRightVector();
        auto vector1rhs   = matrix1->getRightVector();
        auto vector2rhs   = matrix2->getRightVector();
        auto vectorresult = matrix2->getRightVector();

        fillWithPseudoLaplacian<FACTORY>( matrix1 );
        fillWithPseudoLaplacian<FACTORY>( matrix2 );

        matrix2->scale( 1.234567 ); // matrix2 = matrix1

        vector1lhs->setRandomValues();
        vector2lhs->copyVector( vector1lhs );

        matrix1->mult( vector1lhs, vector1rhs );
        matrix2->mult( vector2lhs, vector2rhs ); // vector2rhs = 1.234567vector1rhs

        vector1rhs->scale( 1.234567 );
        vectorresult->subtract( *vector1rhs, *vector2rhs );

        if ( vectorresult->L1Norm() < 0.000001 )
            utils->passes( "matrices are equally scaled" );
        else
            utils->failure( "matrices are equally scaled" );
        if ( vector1rhs->L1Norm() > 0.00001 )
            utils->passes( "non-trivial vector" );
        else
            utils->passes( "trivial vector" );
        PROFILE_STOP( "VerifyScaleMatrix" );
    }
};


template<typename FACTORY>
class VerifyExtractDiagonal
{
public:
    static const char *get_test_name() { return "Verify extractDiagonal"; }

    static void run_test( AMP::UnitTest *utils )
    {
        PROFILE_START( "VerifyExtractDiagonal" );
        auto matrix     = FACTORY::getMatrix();
        auto vector     = matrix->getRightVector();
        size_t firstRow = vector->getCommunicationList()->getStartGID();
        size_t maxCols  = matrix->numGlobalColumns();
        for ( size_t i = 0; i != vector->getCommunicationList()->numLocalRows(); i++ ) {
            int row = static_cast<int>( i + firstRow );
            if ( row >= static_cast<int>( maxCols ) )
                break;
            matrix->setValueByGlobalID( row, row, static_cast<double>( row + 1 ) );
        }
        auto diag         = matrix->extractDiagonal();
        double l1norm     = static_cast<double>( diag->L1Norm() );
        double numRows    = static_cast<double>( matrix->numGlobalRows() );
        double numCols    = static_cast<double>( matrix->numGlobalColumns() );
        double compareVal = std::min( numRows, numCols );
        compareVal        = compareVal * ( 1. + compareVal ) / 2.;
        if ( fabs( l1norm - compareVal ) < 0.000001 )
            utils->passes( "Verify extractDiagonal" );
        else
            utils->failure( "Verify extractDiagonal" );
        PROFILE_STOP( "VerifyExtractDiagonal" );
    }
};


template<typename FACTORY>
class VerifyMultMatrix
{
public:
    static const char *get_test_name() { return "verify mult"; }

    static void run_test( AMP::UnitTest *utils )
    {
        PROFILE_START( "VerifyMultMatrix" );
        auto matrix    = FACTORY::getMatrix();
        auto vectorlhs = matrix->getRightVector();
        auto vectorrhs = matrix->getRightVector();
        double normlhs, normrhs;

        // Verify 0 matrix from factory
        if ( matrix->L1Norm() == 0.0 )
            utils->passes( "Factory returns 0 matrix" );
        else
            utils->failure( "Factory returns 0 matrix" );

        // Verify mult with 0 matrix
        matrix->zero();
        vectorlhs->setRandomValues();
        matrix->mult( vectorlhs, vectorrhs );
        normlhs = static_cast<double>( vectorlhs->L2Norm() );
        normrhs = static_cast<double>( vectorrhs->L2Norm() );
        if ( ( normlhs > 0 ) && ( normrhs < 0.0000001 ) )
            utils->passes( "mult by 0 matrix" );
        else
            utils->failure( "mult by 0 matrix" );

        // Verify mult with identity
        vectorlhs->setToScalar( 1.0 );
        matrix->setDiagonal( vectorlhs );
        vectorlhs->setRandomValues();
        matrix->mult( vectorlhs, vectorrhs );
        normlhs = static_cast<double>( vectorlhs->L2Norm() );
        vectorrhs->subtract( *vectorlhs, *vectorrhs );
        normrhs = static_cast<double>( vectorrhs->L2Norm() );
        if ( ( normlhs > 0 ) && ( normrhs < 0.0000001 ) )
            utils->passes( "mult by I matrix" );
        else
            utils->failure( "mult by I matrix" );

        // Try the non-trivial matrix
        fillWithPseudoLaplacian<FACTORY>( matrix );
        vectorlhs->setRandomValues();
        matrix->mult( vectorlhs, vectorrhs );
        normlhs = static_cast<double>( vectorlhs->L2Norm() );
        normrhs = static_cast<double>( vectorrhs->L2Norm() );
        if ( ( normlhs > 0 ) && ( normrhs > 0 ) )
            utils->passes( "mult by other matrix" );
        else
            utils->failure( "mult by other matrix" );
        PROFILE_STOP( "VerifyMultMatrix" );
    }
};


// Test matrix-matrix multiplication (this tests takes a long time for large matrices)
template<typename FACTORY>
class VerifyMatMultMatrix
{
public:
    static const char *get_test_name() { return "verify mult"; }

    static void run_test( AMP::UnitTest *utils )
    {
        PROFILE_START( "VerifyMatMultMatrix" );
        auto matZero   = FACTORY::getMatrix();
        auto matIdent  = FACTORY::getMatrix();
        auto matLaplac = FACTORY::getMatrix();
        auto matSol    = FACTORY::getMatrix();
        auto vector1   = matZero->getRightVector();
        auto vector2   = matZero->getRightVector();
        auto vector3   = matZero->getRightVector();

        if ( vector1->getGlobalSize() > 1000 ) {
            // Matrix-matrix multiplies take a long time (skip it)
            PROFILE_STOP2( "VerifyMatMultMatrix" );
            return;
        }

        // Create the matrices and vectors of interest
        matZero->zero();
        matIdent->zero();
        vector2->setToScalar( 1.0 );
        matIdent->setDiagonal( vector2 );
        fillWithPseudoLaplacian<FACTORY>( matLaplac );
        vector1->setRandomValues();
        double ans1, ans2, ans3;

        // Verify matMultiply with 0 matrix
        matSol = AMP::LinearAlgebra::Matrix::matMultiply( matZero, matLaplac );
        if ( matSol->L1Norm() == 0.0 )
            utils->passes( "matMultiply with 0 matrix" );
        else
            utils->failure( "matMultiply with 0 matrix" );

        // Verify mult with identity
        matLaplac->mult( vector1, vector2 );
        ans1   = static_cast<double>( vector2->L2Norm() );
        matSol = AMP::LinearAlgebra::Matrix::matMultiply( matIdent, matLaplac );
        matSol->mult( vector1, vector2 );
        ans2   = static_cast<double>( vector2->L2Norm() );
        matSol = AMP::LinearAlgebra::Matrix::matMultiply( matLaplac, matIdent );
        matSol->mult( vector1, vector2 );
        ans3 = static_cast<double>( vector2->L2Norm() );
        if ( AMP::Utilities::approx_equal( ans1, ans2 ) &&
             AMP::Utilities::approx_equal( ans1, ans3 ) && ans1 != 0.0 )
            utils->passes( "matMultiply with identity matrix" );
        else
            utils->failure( "matMultiply with identity matrix" );

        // Verify mult with two trival matrices
        matLaplac->mult( vector1, vector2 );
        matLaplac->mult( vector2, vector3 );
        ans1   = static_cast<double>( vector3->L2Norm() );
        matSol = AMP::LinearAlgebra::Matrix::matMultiply( matLaplac, matLaplac );
        matSol->mult( vector1, vector2 );
        ans2 = static_cast<double>( vector2->L2Norm() );
        if ( AMP::Utilities::approx_equal( ans1, ans2 ) && ans1 != 0.0 )
            utils->passes( "matMultiply with trival matrix" );
        else
            utils->failure( "matMultiply with trival matrix" );

        PROFILE_STOP( "VerifyMatMultMatrix" );
    }
};


template<typename FACTORY>
class VerifyAddElementNode
{
public:
    static const char *get_test_name() { return "verify set nodes by element"; }

    static void run_test( AMP::UnitTest *utils )
    {
        PROFILE_START( "VerifySetElementNode" );
        auto mesh   = FACTORY::getMesh();
        auto dofmap = FACTORY::getDOFMap();
        auto matrix = FACTORY::getMatrix();
        matrix->zero();

        // Fill all the node-node entries
        auto it  = mesh->getIterator( AMP::Mesh::GeomType::Volume, 0 );
        auto end = it.end();
        std::vector<size_t> dofs;
        dofs.reserve( 24 );
        while ( it != end ) {
            auto nodes = it->getElements( AMP::Mesh::GeomType::Vertex );
            dofs.clear();
            for ( auto &node : nodes ) {
                std::vector<size_t> dofsNode;
                dofmap->getDOFs( node.globalID(), dofsNode );
                for ( auto &elem : dofsNode )
                    dofs.push_back( elem );
            }
            for ( size_t r = 0; r < dofs.size(); r++ ) {
                for ( size_t c = 0; c < dofs.size(); c++ ) {
                    double val = -1.0;
                    if ( r == c )
                        val = dofs.size() - 1;
                    matrix->addValueByGlobalID( dofs[r], dofs[c], val );
                }
            }
            ++it;
        }
        matrix->makeConsistent();

        // Call makeConsistent a second time
        // This can illustrate a bug where the fill pattern of remote data has changed
        //   and epetra maintains the list of remote rows, but updates the columns
        //   resulting in an access error using the std::vector
        // Another example of this bug can be found in extra_tests/test_Epetra_FECrsMatrix_bug
        // Note: there is no point in catching this bug with a try catch since a failure
        //   will cause asymettric behavior that create a deadlock with one process waiting
        //   for the failed process
        // The current workaround is to disable the GLIBCXX_DEBUG flags?
        matrix->makeConsistent();

        // Check the values
        bool pass = true;
        it        = mesh->getIterator( AMP::Mesh::GeomType::Vertex, 0 );
        end       = it.end();
        std::vector<size_t> cols;
        std::vector<double> values;
        while ( it != end ) {
            dofmap->getDOFs( it->globalID(), dofs );
            for ( auto &dof : dofs ) {
                matrix->getRowByGlobalID( dof, cols, values );
                double sum = 0.0;
                for ( auto &value : values )
                    sum += value;
                if ( fabs( sum ) > 1e-14 || cols.empty() )
                    pass = false;
            }
            ++it;
        }
        if ( pass )
            utils->passes( "VerifySetElementNode" );
        else
            utils->failure( "VerifySetElementNode" );

        PROFILE_STOP( "VerifySetElementNode" );
    }
};
} // namespace LinearAlgebra
} // namespace AMP

#endif
