
#include "DirichletMatrixCorrection.h"
#include "utils/InputDatabase.h"
#include "utils/Utilities.h"

namespace AMP {
namespace Operator {


/****************************************************************
* Constructors                                                  *
****************************************************************/
DirichletMatrixCorrection::DirichletMatrixCorrection(
    const AMP::shared_ptr<DirichletMatrixCorrectionParameters> &params )
    : BoundaryOperator( params )
{
    d_variable                 = params->d_variable;
    d_computedAddRHScorrection = false;
    d_symmetricCorrection  = ( params->d_db )->getBoolWithDefault( "symmetric_correction", true );
    d_zeroDirichletBlock   = ( params->d_db )->getBoolWithDefault( "zero_dirichlet_block", false );
    d_skipRHSsetCorrection = ( params->d_db )->getBoolWithDefault( "skip_rhs_correction", true );
    d_skipRHSaddCorrection =
        ( params->d_db )->getBoolWithDefault( "skip_rhs_add_correction", d_skipRHSsetCorrection );
    d_applyMatrixCorrectionWasCalled = false;

    reset( params );
}


/****************************************************************
* Reset                                                         *
****************************************************************/
void DirichletMatrixCorrection::reset( const AMP::shared_ptr<OperatorParameters> &params )
{
    AMP::shared_ptr<DirichletMatrixCorrectionParameters> myParams =
        AMP::dynamic_pointer_cast<DirichletMatrixCorrectionParameters>( params );

    AMP_INSIST( ( ( myParams.get() ) != nullptr ), "NULL parameters" );

    parseParams( myParams );

    d_inputMatrix = myParams->d_inputMatrix;
    AMP_INSIST( ( ( d_inputMatrix.get() ) != nullptr ), "NULL matrix" );
    d_inputMatrix->makeConsistent(); // Check that we can call makeConsistent

    if ( d_skipRHSsetCorrection ) {
        AMP_ASSERT( d_skipRHSaddCorrection );
    }
    if ( !d_skipRHSaddCorrection ) {
        AMP_ASSERT( !d_skipRHSsetCorrection );
    }

    d_applyMatrixCorrectionWasCalled = false;

    initRhsCorrectionSet();

    if ( d_skipRHSaddCorrection ) {
        applyMatrixCorrection();
    } // end if
}
void DirichletMatrixCorrection::parseParams(
    const AMP::shared_ptr<DirichletMatrixCorrectionParameters> &params )
{
    AMP_INSIST( ( ( ( params->d_db ).get() ) != nullptr ), "NULL database" );
    bool skipParams = ( params->d_db )->getBoolWithDefault( "skip_params", false );

    if ( !skipParams ) {
        d_symmetricCorrection =
            ( params->d_db )->getBoolWithDefault( "symmetric_correction", true );
        d_zeroDirichletBlock =
            ( params->d_db )->getBoolWithDefault( "zero_dirichlet_block", false );

        d_skipRHSsetCorrection =
            ( params->d_db )->getBoolWithDefault( "skip_rhs_correction", true );
        d_skipRHSaddCorrection =
            ( params->d_db )
                ->getBoolWithDefault( "skip_rhs_add_correction", d_skipRHSsetCorrection );

        if ( d_symmetricCorrection == false ) {
            d_skipRHSaddCorrection = true;
        }

        AMP_INSIST( ( params->d_db )->keyExists( "number_of_ids" ),
                    "Key ''number_of_ids'' is missing!" );
        int numIds = ( params->d_db )->getInteger( "number_of_ids" );

        d_boundaryIds.resize( numIds );
        d_dofIds.resize( numIds );

        char key[100];
        for ( int j = 0; j < numIds; ++j ) {
            sprintf( key, "id_%d", j );
            AMP_INSIST( ( params->d_db )->keyExists( key ), "Key is missing!" );
            d_boundaryIds[j] = ( params->d_db )->getInteger( key );

            sprintf( key, "number_of_dofs_%d", j );
            AMP_INSIST( ( params->d_db )->keyExists( key ), "Key is missing!" );
            int numDofIds = ( params->d_db )->getInteger( key );

            d_dofIds[j].resize( numDofIds );
            for ( int i = 0; i < numDofIds; ++i ) {
                sprintf( key, "dof_%d_%d", j, i );
                AMP_INSIST( ( params->d_db )->keyExists( key ), "Key is missing!" );
                d_dofIds[j][i] = ( params->d_db )->getInteger( key );
            } // end for i
        }     // end for j

        if ( !d_skipRHSsetCorrection ) {
            d_dirichletValues.resize( numIds );
            for ( int j = 0; j < numIds; ++j ) {
                int numDofIds = d_dofIds[j].size();
                d_dirichletValues[j].resize( numDofIds );

                for ( int i = 0; i < numDofIds; ++i ) {
                    sprintf( key, "value_%d_%d", j, i );
                    d_dirichletValues[j][i] = ( params->d_db )->getDoubleWithDefault( key, 0.0 );
                } // end for i
            }     // end for j
        }
    }
}


/****************************************************************
* applyMatrixCorrection                                         *
****************************************************************/
void DirichletMatrixCorrection::applyMatrixCorrection()
{
    AMP_ASSERT( !d_applyMatrixCorrectionWasCalled );
    d_applyMatrixCorrectionWasCalled = true;

    AMP::LinearAlgebra::Vector::shared_ptr inVec        = d_inputMatrix->getRightVector();
    AMP::Discretization::DOFManager::shared_ptr dof_map = inVec->getDOFManager();
    AMP_ASSERT( ( *dof_map ) == ( *d_inputMatrix->getLeftDOFManager() ) );
    AMP_ASSERT( ( *dof_map ) == ( *d_inputMatrix->getRightDOFManager() ) );

    for ( size_t k = 0; k < d_boundaryIds.size(); ++k ) {
        AMP::Mesh::MeshIterator bnd =
            d_Mesh->getBoundaryIDIterator( AMP::Mesh::Vertex, d_boundaryIds[k], 0 );
        AMP::Mesh::MeshIterator end_bnd = bnd.end();

        for ( ; bnd != end_bnd; ++bnd ) {
            std::vector<size_t> bndDofIds;
            dof_map->getDOFs( bnd->globalID(), bndDofIds );

            // Get neighbors does not include the calling node (bnd) itself.
            // Get neighbors also returns remote neighbors
            // The calling node (bnd) must be owned locally.
            std::vector<AMP::Mesh::MeshElement::shared_ptr> neighbors = bnd->getNeighbors();
            for ( unsigned int i = 0; i < neighbors.size(); ++i ) {
                AMP_ASSERT( ( *( neighbors[i] ) ) != ( *bnd ) );
            } // end for i

            for ( unsigned int j = 0; j < d_dofIds[k].size(); ++j ) {
                for ( unsigned int i = 0; i < bndDofIds.size(); ++i ) {
                    if ( d_dofIds[k][j] == i ) {
                        if ( d_zeroDirichletBlock ) {
                            d_inputMatrix->setValueByGlobalID( bndDofIds[i], bndDofIds[i], 0.0 );
                        } else {
                            d_inputMatrix->setValueByGlobalID( bndDofIds[i], bndDofIds[i], 1.0 );
                        }
                    } else {
                        d_inputMatrix->setValueByGlobalID(
                            bndDofIds[d_dofIds[k][j]], bndDofIds[i], 0.0 );
                        if ( d_symmetricCorrection ) {
                            d_inputMatrix->setValueByGlobalID(
                                bndDofIds[i], bndDofIds[d_dofIds[k][j]], 0.0 );
                        }
                    }
                } // end for i
                for ( size_t n = 0; n < neighbors.size(); ++n ) {
                    std::vector<size_t> nhDofIds;
                    dof_map->getDOFs( neighbors[n]->globalID(), nhDofIds );
                    for ( unsigned int i = 0; i < nhDofIds.size(); ++i ) {
                        d_inputMatrix->setValueByGlobalID(
                            bndDofIds[d_dofIds[k][j]], nhDofIds[i], 0.0 );
                        if ( d_symmetricCorrection ) {
                            d_inputMatrix->setValueByGlobalID(
                                nhDofIds[i], bndDofIds[d_dofIds[k][j]], 0.0 );
                        }
                    } // end for i
                }     // end for n
            }         // end for j
        }             // end for bnd
    }                 // end for k

    // This does consistent for both "Sum-into" and "set".
    d_inputMatrix->makeConsistent();
}


void DirichletMatrixCorrection::initRhsCorrectionSet()
{
    if ( !d_skipRHSsetCorrection ) {
        int numIds = d_dofIds.size();
        char key[100];
        AMP::shared_ptr<AMP::InputDatabase> tmp_db( new AMP::InputDatabase( "Dummy" ) );
        tmp_db->putBool( "skip_params", false );
        tmp_db->putBool( "isAttachedToVolumeOperator", false );
        tmp_db->putInteger( "number_of_ids", numIds );
        tmp_db->putInteger( "print_info_level", d_iDebugPrintInfoLevel );
        for ( int j = 0; j < numIds; j++ ) {
            int numDofIds = d_dofIds[j].size();

            sprintf( key, "id_%d", j );
            tmp_db->putInteger( key, d_boundaryIds[j] );

            sprintf( key, "number_of_dofs_%d", j );
            tmp_db->putInteger( key, numDofIds );

            for ( int i = 0; i < numDofIds; i++ ) {
                sprintf( key, "dof_%d_%d", j, i );
                tmp_db->putInteger( key, d_dofIds[j][i] );

                sprintf( key, "value_%d_%d", j, i );
                tmp_db->putDouble( key, d_dirichletValues[j][i] );
            } // end for i
        }     // end for j

        AMP::shared_ptr<DirichletVectorCorrectionParameters> setDispOpParams(
            new DirichletVectorCorrectionParameters( tmp_db ) );
        setDispOpParams->d_variable = d_variable;
        setDispOpParams->d_Mesh     = d_Mesh;

        if ( d_rhsCorrectionSet.get() == nullptr ) {
            d_rhsCorrectionSet.reset( new DirichletVectorCorrection( setDispOpParams ) );
        } else {
            d_rhsCorrectionSet->reset( setDispOpParams );
        }
    }
}


void DirichletMatrixCorrection::initRhsCorrectionAdd( AMP::LinearAlgebra::Vector::shared_ptr rhs )
{
    AMP_ASSERT( !d_applyMatrixCorrectionWasCalled );

    if ( !d_skipRHSsetCorrection ) {
        if ( !d_skipRHSaddCorrection ) {
            if ( d_dispVals.get() == nullptr ) {
                d_dispVals = ( subsetOutputVector( rhs ) )->cloneVector();
                AMP_ASSERT( ( *( d_dispVals->getVariable() ) ) == ( *d_variable ) );
            }

            d_dispVals->zero();

            AMP::LinearAlgebra::Vector::shared_ptr emptyVec;
            d_rhsCorrectionSet->apply( emptyVec, d_dispVals );

            if ( d_rhsCorrectionAdd.get() == nullptr ) {
                d_rhsCorrectionAdd = d_dispVals->cloneVector();
            }

            d_inputMatrix->mult( d_dispVals, d_rhsCorrectionAdd );

            d_rhsCorrectionAdd->scale( -1.0 );
        }
    }
}


void DirichletMatrixCorrection::addRHScorrection( AMP::LinearAlgebra::Vector::shared_ptr rhs )
{
    if ( !d_skipRHSaddCorrection ) {
        if ( !d_applyMatrixCorrectionWasCalled ) {
            initRhsCorrectionAdd( rhs );
            applyMatrixCorrection();
        } // end if
        AMP::LinearAlgebra::Vector::shared_ptr myRhs = subsetOutputVector( rhs );
        myRhs->add( myRhs, d_rhsCorrectionAdd );
    }
}


void DirichletMatrixCorrection::setRHScorrection( AMP::LinearAlgebra::Vector::shared_ptr rhs )
{
    if ( !d_skipRHSsetCorrection ) {
        AMP::LinearAlgebra::Vector::shared_ptr emptyVec;
        d_rhsCorrectionSet->apply( emptyVec, rhs );
    }
}
}
}
