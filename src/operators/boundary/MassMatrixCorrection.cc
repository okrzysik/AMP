#include "MassMatrixCorrection.h"
#include "utils/InputDatabase.h"
#include "utils/Utilities.h"

namespace AMP {
namespace Operator {

void MassMatrixCorrection::resetBoundaryIds(
    const AMP::shared_ptr<MassMatrixCorrectionParameters> &params )
{
    AMP_INSIST( ( ( ( params->d_db ).get() ) != NULL ), "NULL database" );
    bool skipParams = ( params->d_db )->getBoolWithDefault( "skip_params", true );
    d_bSetIdentityOnDiagonal =
        ( params->d_db )->getBoolWithDefault( "setIdentityOnDiagonal", false );

    if ( !skipParams ) {
        AMP_INSIST( ( params->d_db )->keyExists( "number_of_ids" ),
                    "Key ''number_of_ids'' is missing!" );
        int numIds = ( params->d_db )->getInteger( "number_of_ids" );

        d_boundaryIds.resize( numIds );
        d_dofIds.resize( numIds );

        char key[100];
        for ( int j = 0; j < numIds; j++ ) {
            sprintf( key, "id_%d", j );
            AMP_INSIST( ( params->d_db )->keyExists( key ), "Key is missing!" );
            d_boundaryIds[j] = ( params->d_db )->getInteger( key );

            sprintf( key, "number_of_dofs_%d", j );
            AMP_INSIST( ( params->d_db )->keyExists( key ), "Key is missing!" );
            int numDofIds = ( params->d_db )->getInteger( key );

            d_dofIds[j].resize( numDofIds );
            for ( int i = 0; i < numDofIds; i++ ) {
                sprintf( key, "dof_%d_%d", j, i );
                AMP_INSIST( ( params->d_db )->keyExists( key ), "Key is missing!" );
                d_dofIds[j][i] = ( params->d_db )->getInteger( key );
            } // end for i
        }     // end for j
    }
}

void MassMatrixCorrection::reset( const AMP::shared_ptr<OperatorParameters> &params )
{


    AMP::shared_ptr<MassMatrixCorrectionParameters> myParams =
        AMP::dynamic_pointer_cast<MassMatrixCorrectionParameters>( params );

    AMP_INSIST( ( ( myParams.get() ) != NULL ), "NULL parameters" );

    resetBoundaryIds( myParams );

    double diagVal = d_bSetIdentityOnDiagonal ? 1.0 : 0.0;

    AMP::LinearAlgebra::Matrix::shared_ptr inputMatrix = myParams->d_inputMatrix;
    AMP_INSIST( ( ( inputMatrix.get() ) != NULL ), "NULL matrix" );

    AMP::LinearAlgebra::Vector::shared_ptr inVec        = inputMatrix->getRightVector();
    AMP::Discretization::DOFManager::shared_ptr dof_map = inVec->getDOFManager();

    unsigned int numIds = d_boundaryIds.size();

    for ( unsigned int k = 0; k < numIds; k++ ) {

        AMP::Mesh::MeshIterator bnd =
            d_Mesh->getBoundaryIDIterator( AMP::Mesh::Vertex, d_boundaryIds[k], 0 );
        AMP::Mesh::MeshIterator end_bnd = bnd.end();

        for ( ; bnd != end_bnd; ++bnd ) {
            std::vector<size_t> bndGlobalIds;
            dof_map->getDOFs( bnd->globalID(), bndGlobalIds );

            std::vector<AMP::Mesh::MeshElement::shared_ptr> neighbors = bnd->getNeighbors();

            for ( size_t i = 0; i < neighbors.size(); ++i ) {
                AMP_ASSERT( ( *( neighbors[i] ) ) != ( *bnd ) );
            } // end for el

            for ( unsigned int j = 0; j < d_dofIds[k].size(); ++j ) {
                for ( unsigned int i = 0; i < bndGlobalIds.size(); ++i ) {

                    if ( d_dofIds[k][j] == i ) {
                        inputMatrix->setValueByGlobalID(
                            bndGlobalIds[j], bndGlobalIds[j], diagVal );
                    } else {
                        inputMatrix->setValueByGlobalID(
                            bndGlobalIds[d_dofIds[k][j]], bndGlobalIds[i], 0.0 );
                        inputMatrix->setValueByGlobalID(
                            bndGlobalIds[i], bndGlobalIds[d_dofIds[k][j]], 0.0 );
                    }
                } // end for i
                for ( size_t n = 0; n < neighbors.size(); ++n ) {
                    std::vector<size_t> nhDofIds;
                    dof_map->getDOFs( neighbors[n]->globalID(), nhDofIds );
                    for ( unsigned int i = 0; i < nhDofIds.size(); ++i ) {
                        inputMatrix->setValueByGlobalID(
                            bndGlobalIds[d_dofIds[k][j]], nhDofIds[i], 0.0 );
                        inputMatrix->setValueByGlobalID(
                            nhDofIds[i], bndGlobalIds[d_dofIds[k][j]], 0.0 );
                    } // end for i
                }     // end for n
            }         // end for j
        }             // end for bnd
    }                 // end for k

    // This does consistent for both "Sum-into" and "set".
    inputMatrix->makeConsistent();
}
}
}
