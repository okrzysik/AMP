#include "AMP/matrices/testHelpers/MatrixTests.h"
#include "AMP/matrices/Matrix.h"
#include "AMP/matrices/MatrixBuilder.h"
#include "AMP/matrices/testHelpers/test_MatrixVectorFactory.h"
#include "AMP/utils/Utilities.h"

#include "ProfilerApp.h"


namespace AMP::LinearAlgebra {


void fillWithPseudoLaplacian( std::shared_ptr<AMP::LinearAlgebra::Matrix> matrix,
                              std::shared_ptr<AMP::Discretization::DOFManager> dofmap )
{
    if ( matrix->type() == "NativePetscMatrix" ) {
        std::map<size_t, std::vector<size_t>> allCols;
        std::map<size_t, std::vector<double>> allVals;
        for ( size_t i = dofmap->beginDOF(); i != dofmap->endDOF(); i++ ) {
            const auto cols  = matrix->getColumnIDs( i );
            const auto ncols = cols.size();
            std::vector<double> vals( ncols );
            for ( size_t j = 0; j != ncols; j++ ) {
                if ( cols[j] == i )
                    vals[j] = static_cast<double>( ncols );
                else
                    vals[j] = -1;
            }
            allVals[i] = vals;
            allCols[i] = cols;
        }
        for ( size_t i = dofmap->beginDOF(); i != dofmap->endDOF(); i++ ) {
            auto &cols = allCols[i];
            auto &vals = allVals[i];
            matrix->setValuesByGlobalID( 1, cols.size(), &i, cols.data(), vals.data() );
        }
    } else {

        for ( size_t i = dofmap->beginDOF(); i != dofmap->endDOF(); i++ ) {
            auto cols        = matrix->getColumnIDs( i );
            const auto ncols = cols.size();
            std::vector<double> vals( ncols );
            for ( size_t j = 0; j != ncols; j++ ) {
                if ( cols[j] == i )
                    vals[j] = static_cast<double>( ncols );
                else
                    vals[j] = -1;
            }
            if ( ncols ) {
                matrix->setValuesByGlobalID<double>( 1, ncols, &i, cols.data(), vals.data() );
            }
        }
    }

    matrix->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD );
}

static void fillWithPseudoLaplacian( std::shared_ptr<AMP::LinearAlgebra::Matrix> matrix,
                                     std::shared_ptr<const MatrixFactory> factory )
{
    auto dofmap = factory->getDOFMap();
    fillWithPseudoLaplacian( matrix, dofmap );
}

void MatrixTests::InstantiateMatrix( AMP::UnitTest *utils )
{
    PROFILE_START( "InstantiateMatrix" );
    auto matrix = d_factory->getMatrix();
    if ( matrix )
        utils->passes( "created " + matrix->type() );
    else
        utils->failure( "created " );
    PROFILE_STOP( "InstantiateMatrix" );
}


void MatrixTests::VerifyGetLeftRightVector( AMP::UnitTest *utils )
{
    PROFILE_START( "VerifyGetLeftRightVector" );
    auto matrix   = d_factory->getMatrix();
    auto factory1 = std::make_shared<AmpInterfaceRightVectorFactory>( matrix );
    auto factory2 = std::make_shared<AmpInterfaceLeftVectorFactory>( matrix );
    VectorTests tests1( factory1 );
    VectorTests tests2( factory2 );
    tests1.testBasicVector( utils );
    tests2.testBasicVector( utils );
#if defined( AMP_USE_PETSC )
    auto factory3 = std::make_shared<PETScInterfaceRightVectorFactory>( matrix );
    auto factory4 = std::make_shared<PETScInterfaceLeftVectorFactory>( matrix );
    VectorTests tests3( factory3 );
    VectorTests tests4( factory4 );
    tests3.testPetsc( utils );
    tests4.testPetsc( utils );
#endif
    PROFILE_STOP( "VerifyGetLeftRightVector" );
}


void MatrixTests::VerifyGetSetValuesMatrix( AMP::UnitTest *utils )
{
    PROFILE_START( "VerifyGetSetValuesMatrix" );
    auto matrix = d_factory->getMatrix();
    auto dofmap = d_factory->getDOFMap();

    //    matrix->makeConsistent();
    fillWithPseudoLaplacian( matrix, d_factory );
    for ( size_t i = dofmap->beginDOF(); i != dofmap->endDOF(); i++ ) {
        std::vector<size_t> cols;
        std::vector<double> vals;
        matrix->getRowByGlobalID( i, cols, vals );
        for ( size_t j = 0; j != cols.size(); j++ ) {
            double ans   = ( i == cols[j] ) ? cols.size() : -1.;
            double value = matrix->getValueByGlobalID( i, cols[j] );
            if ( vals[j] != ans || value != vals[j] ) {
                utils->failure( "bad value in matrix " + matrix->type() );
                return;
            }
        }
    }
    utils->passes( "verify get and set" + matrix->type() );
    PROFILE_STOP( "VerifyGetSetValuesMatrix" );
}


void MatrixTests::VerifyAXPYMatrix( AMP::UnitTest *utils )
{
    PROFILE_START( "VerifyAXPYMatrix" );

    // Create vectors/matricies from the factory
    auto matrix1      = d_factory->getMatrix();
    auto matrix2      = d_factory->getMatrix();
    auto vector1lhs   = matrix1->getRightVector();
    auto vector2lhs   = matrix2->getRightVector();
    auto vector1rhs   = matrix1->getRightVector();
    auto vector2rhs   = matrix2->getRightVector();
    auto vectorresult = matrix2->getRightVector();
    fillWithPseudoLaplacian( matrix1, d_factory );
    fillWithPseudoLaplacian( matrix2, d_factory );

    // Test axpy
    matrix2->axpy( -2., matrix1 ); // matrix2 = -matrix1
    vector1lhs->setRandomValues();
    vector2lhs->copyVector( vector1lhs );
    matrix1->mult( vector1lhs, vector1rhs );
    matrix2->mult( vector2lhs, vector2rhs ); // vector2rhs = - vector1rhs
    vectorresult->add( *vector1rhs, *vector2rhs );
    if ( vectorresult->L1Norm() < 0.000001 )
        utils->passes( "matrices are opposite " + matrix1->type() );
    else
        utils->failure( "matrices are not opposite " + matrix1->type() );
    if ( vector1rhs->L1Norm() > 0.00001 )
        utils->passes( "non-trivial vector " + matrix1->type() );
    else
        utils->passes( "trivial vector " + matrix1->type() );

        // currently can't seem to correctly catch signal that PETSc throws
#if !defined( AMP_USE_PETSC )
    // Test that axpy failes with different sized matricies
    std::vector<size_t> row( 7 );
    for ( size_t i = 0; i < row.size(); i++ )
        row[i] = i;
    auto smallVec = AMP::LinearAlgebra::createSimpleVector<double>( 7, vector1lhs->getVariable() );
    auto smallMat = AMP::LinearAlgebra::createMatrix(
        smallVec, smallVec, d_factory->type(), [row]( size_t ) { return row; } );
    try {
        matrix2->axpy( -2., smallMat ); // matrix2 = -matrix1
        utils->failure( "axpy did not crash with different sized matrices " + matrix1->type() );
    } catch ( std::exception &err ) {
        utils->passes( "axpy correctly fails with different sized matrices " + matrix1->type() );
    } catch ( ... ) {
        utils->failure( "axpy fails with different sized matrices (unknown failure) " +
                        matrix1->type() );
    }
#endif
    PROFILE_STOP( "VerifyAXPYMatrix" );
}


void MatrixTests::VerifyScaleMatrix( AMP::UnitTest *utils )
{
    PROFILE_START( "VerifyScaleMatrix" );
    auto matrix1      = d_factory->getMatrix();
    auto matrix2      = d_factory->getMatrix();
    auto vector1lhs   = matrix1->getRightVector();
    auto vector2lhs   = matrix2->getRightVector();
    auto vector1rhs   = matrix1->getRightVector();
    auto vector2rhs   = matrix2->getRightVector();
    auto vectorresult = matrix2->getRightVector();

    fillWithPseudoLaplacian( matrix1, d_factory );
    fillWithPseudoLaplacian( matrix2, d_factory );

    matrix2->scale( 1.234567 ); // matrix2 = matrix1

    vector1lhs->setRandomValues();
    vector2lhs->copyVector( vector1lhs );

    matrix1->mult( vector1lhs, vector1rhs );
    matrix2->mult( vector2lhs, vector2rhs ); // vector2rhs = 1.234567vector1rhs

    vector1rhs->scale( 1.234567 );
    vectorresult->subtract( *vector1rhs, *vector2rhs );

    if ( vectorresult->L1Norm() < 0.000001 )
        utils->passes( "matrices are equally scaled " + matrix1->type() );
    else
        utils->failure( "matrices are equally scaled " + matrix1->type() );
    if ( vector1rhs->L1Norm() > 0.00001 )
        utils->passes( "non-trivial vector " + matrix1->type() );
    else
        utils->passes( "trivial vector " + matrix1->type() );
    PROFILE_STOP( "VerifyScaleMatrix" );
}


void MatrixTests::VerifyExtractDiagonal( AMP::UnitTest *utils )
{
    PROFILE_START( "VerifyExtractDiagonal" );
    auto matrix = d_factory->getMatrix();
    //    matrix->makeConsistent(); // required by PETSc
    auto vector     = matrix->getRightVector();
    size_t firstRow = vector->getCommunicationList()->getStartGID();
    size_t maxCols  = matrix->numGlobalColumns();
    for ( size_t i = 0; i != vector->getCommunicationList()->numLocalRows(); i++ ) {
        int row = static_cast<int>( i + firstRow );
        if ( row >= static_cast<int>( maxCols ) )
            break;
        matrix->setValueByGlobalID( row, row, static_cast<double>( row + 1 ) );
    }
    matrix->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD ); // required by PETSc
    auto diag         = matrix->extractDiagonal();
    double l1norm     = static_cast<double>( diag->L1Norm() );
    auto numRows      = static_cast<double>( matrix->numGlobalRows() );
    auto numCols      = static_cast<double>( matrix->numGlobalColumns() );
    double compareVal = std::min( numRows, numCols );
    compareVal        = compareVal * ( 1. + compareVal ) / 2.;
    if ( fabs( l1norm - compareVal ) < 0.000001 )
        utils->passes( "Verify extractDiagonal " + matrix->type() );
    else
        utils->failure( "Verify extractDiagonal " + matrix->type() );
    PROFILE_STOP( "VerifyExtractDiagonal" );
}


void MatrixTests::VerifyMultMatrix( AMP::UnitTest *utils )
{
    PROFILE_START( "VerifyMultMatrix" );
    auto matrix    = d_factory->getMatrix();
    auto vectorlhs = matrix->getRightVector();
    auto vectorrhs = matrix->getRightVector();
    double normlhs, normrhs;

    // Verify 0 matrix from factory
    if ( matrix->L1Norm() == 0.0 )
        utils->passes( "Factory returns 0 matrix" + matrix->type() );
    else
        utils->failure( "Factory returns 0 matrix" + matrix->type() );

    // Verify mult with 0 matrix
    matrix->zero();
    vectorlhs->setRandomValues();
    matrix->mult( vectorlhs, vectorrhs );
    normlhs = static_cast<double>( vectorlhs->L2Norm() );
    normrhs = static_cast<double>( vectorrhs->L2Norm() );
    if ( ( normlhs > 0 ) && ( normrhs < 0.0000001 ) )
        utils->passes( "mult by 0 matrix " + matrix->type() );
    else
        utils->failure( "mult by 0 matrix " + matrix->type() );

    // Verify mult with identity
    vectorlhs->setToScalar( 1.0 );
    matrix->setDiagonal( vectorlhs );
    vectorlhs->setRandomValues();
    matrix->mult( vectorlhs, vectorrhs );
    normlhs = static_cast<double>( vectorlhs->L2Norm() );
    vectorrhs->subtract( *vectorlhs, *vectorrhs );
    normrhs = static_cast<double>( vectorrhs->L2Norm() );
    if ( ( normlhs > 0 ) && ( normrhs < 0.0000001 ) )
        utils->passes( "mult by I matrix " + matrix->type() );
    else
        utils->failure( "mult by I matrix " + matrix->type() );

    // Try the non-trivial matrix
    fillWithPseudoLaplacian( matrix, d_factory );
    vectorlhs->setRandomValues();
    matrix->mult( vectorlhs, vectorrhs );
    normlhs = static_cast<double>( vectorlhs->L2Norm() );
    normrhs = static_cast<double>( vectorrhs->L2Norm() );
    if ( ( normlhs > 0 ) && ( normrhs > 0 ) )
        utils->passes( "mult by other matrix " + matrix->type() );
    else
        utils->failure( "mult by other matrix " + matrix->type() );
    PROFILE_STOP( "VerifyMultMatrix" );
}


// Test matrix-matrix multiplication (this tests takes a long time for large matrices)
void MatrixTests::VerifyMatMultMatrix( AMP::UnitTest *utils )
{
    PROFILE_START( "VerifyMatMultMatrix" );
    auto matZero   = d_factory->getMatrix();
    auto matIdent  = d_factory->getMatrix();
    auto matLaplac = d_factory->getMatrix();
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
    fillWithPseudoLaplacian( matLaplac, d_factory );
    vector1->setRandomValues();
    vector1->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
    vector2->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
    AMP_ASSERT( static_cast<double>( vector1->L2Norm() ) > 0.0 );

    // Verify matMultiply with 0 matrix
    auto matSol = AMP::LinearAlgebra::Matrix::matMultiply( matZero, matLaplac );
    if ( matSol->L1Norm() == 0.0 )
        utils->passes( "matMultiply with 0 matrix " + matZero->type() );
    else
        utils->failure( "matMultiply with 0 matrix " + matZero->type() );

    // Verify mult with identity
    matLaplac->mult( vector1, vector2 );
    if ( vector2->getUpdateStatus() != AMP::LinearAlgebra::UpdateState::UNCHANGED )
        utils->failure( "matMultiply leaves vector in an inconsistent state" );
    auto ans1 = static_cast<double>( vector2->L2Norm() );
    matSol    = AMP::LinearAlgebra::Matrix::matMultiply( matIdent, matLaplac );
    matSol->mult( vector1, vector2 );
    auto ans2 = static_cast<double>( vector2->L2Norm() );
    matSol    = AMP::LinearAlgebra::Matrix::matMultiply( matLaplac, matIdent );
    matSol->mult( vector1, vector2 );
    auto ans3 = static_cast<double>( vector2->L2Norm() );
    if ( AMP::Utilities::approx_equal( ans1, ans2 ) && AMP::Utilities::approx_equal( ans1, ans3 ) &&
         ans1 != 0.0 )
        utils->passes( "matMultiply with identity matrix " + matSol->type() );
    else
        utils->failure( "matMultiply with identity matrix " + matSol->type() );

    // Verify mult with two trivial matrices
    matLaplac->mult( vector1, vector2 );
    matLaplac->mult( vector2, vector3 );
    ans1   = static_cast<double>( vector3->L2Norm() );
    matSol = AMP::LinearAlgebra::Matrix::matMultiply( matLaplac, matLaplac );
    matSol->mult( vector1, vector2 );
    ans2 = static_cast<double>( vector2->L2Norm() );
    if ( AMP::Utilities::approx_equal( ans1, ans2 ) && ans1 != 0.0 )
        utils->passes( "matMultiply with trivial matrix " + matSol->type() );
    else
        utils->failure( "matMultiply with trivial matrix " + matSol->type() );

    PROFILE_STOP( "VerifyMatMultMatrix" );
}


void MatrixTests::VerifyAddElementNode( AMP::UnitTest *utils )
{
    PROFILE_START( "VerifySetElementNode" );
    auto mesh   = d_factory->getMesh();
    auto dofmap = d_factory->getDOFMap();
    auto matrix = d_factory->getMatrix();
    matrix->zero();

    // Fill all the node-node entries
    auto it  = mesh->getIterator( AMP::Mesh::GeomType::Cell, 0 );
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
    matrix->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD );

    // Call makeConsistent a second time
    // This can illustrate a bug where the fill pattern of remote data has changed
    //   and epetra maintains the list of remote rows, but updates the columns
    //   resulting in an access error using the std::vector
    // Another example of this bug can be found in extra_tests/test_Epetra_FECrsMatrix_bug
    // Note: there is no point in catching this bug with a try catch since a failure
    //   will cause asymettric behavior that create a deadlock with one process waiting
    //   for the failed process
    // The current workaround is to disable the GLIBCXX_DEBUG flags?
    matrix->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD );

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
        utils->passes( "VerifyAddElementNode " + matrix->type() );
    else
        utils->failure( "VerifyAddElementNode " + matrix->type() );

    PROFILE_STOP( "VerifySetElementNode" );
}


} // namespace AMP::LinearAlgebra
