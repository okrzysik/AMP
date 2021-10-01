
#include "DirichletVectorCorrection.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"
#include "DirichletMatrixCorrectionParameters.h"

namespace AMP {
namespace Operator {

void DirichletVectorCorrection::reset( std::shared_ptr<const OperatorParameters> tmpParams )
{
    auto params = std::dynamic_pointer_cast<const DirichletVectorCorrectionParameters>( tmpParams );

    AMP_INSIST( params, "NULL parameters" );
    AMP_INSIST( params->d_db, "NULL database" );

    bool skipParams = params->d_db->getWithDefault<bool>( "skip_params", false );

    if ( !skipParams ) {
        d_scalingFactor = params->d_db->getWithDefault<double>( "SCALING_FACTOR", 1.0 );
        d_setResidual   = params->d_db->getWithDefault<bool>( "setResidual", false );
        d_isAttachedToVolumeOperator =
            params->d_db->getWithDefault<bool>( "isAttachedToVolumeOperator", false );
        d_valuesType = params->d_db->getWithDefault<int>( "valuesType", 1 );
        AMP_INSIST( ( ( d_valuesType == 1 ) || ( d_valuesType == 2 ) ), "Wrong value." );

        AMP_INSIST( params->d_db->keyExists( "number_of_ids" ),
                    "Key ''number_of_ids'' is missing!" );
        int numIds = params->d_db->getScalar<int>( "number_of_ids" );

        d_boundaryIds.resize( numIds );
        d_dofIds.resize( numIds );

        if ( d_valuesType == 1 ) {
            d_dirichletValues1.resize( numIds );
        }

        char key[100];
        for ( int j = 0; j < numIds; j++ ) {
            sprintf( key, "id_%d", j );
            AMP_INSIST( params->d_db->keyExists( key ), "Key is missing!" );
            d_boundaryIds[j] = params->d_db->getScalar<int>( key );

            sprintf( key, "number_of_dofs_%d", j );
            AMP_INSIST( params->d_db->keyExists( key ), "Key is missing!" );
            int numDofIds = params->d_db->getScalar<int>( key );

            d_dofIds[j].resize( numDofIds );

            if ( d_valuesType == 1 ) {
                d_dirichletValues1[j].resize( numDofIds );
            }
            for ( int i = 0; i < numDofIds; i++ ) {
                sprintf( key, "dof_%d_%d", j, i );
                AMP_INSIST( params->d_db->keyExists( key ), "Key is missing!" );
                d_dofIds[j][i] = params->d_db->getScalar<int>( key );

                if ( d_valuesType == 1 ) {
                    sprintf( key, "value_%d_%d", j, i );
                    d_dirichletValues1[j][i] = params->d_db->getWithDefault<double>( key, 0.0 );
                }
            } // end for i
        }     // end for j
    }
}

// This is an in-place apply
void DirichletVectorCorrection::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                                       AMP::LinearAlgebra::Vector::shared_ptr r )
{
    auto rInternal = mySubsetVector( r, d_variable );

    if ( d_iDebugPrintInfoLevel > 3 ) {
        AMP::pout << "L2 Norm of rInternal entering DirichletVectorCorrection::apply is : "
                  << rInternal->L2Norm() << std::endl;
    }

    if ( d_setResidual ) {
        this->applyResidual( u, rInternal );
    } else if ( d_isAttachedToVolumeOperator ) {
        this->applyZeroValues( rInternal );
    } else {
        this->applyNonZeroValues( rInternal );
    }

    rInternal->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );

    if ( d_iDebugPrintInfoLevel > 3 ) {
        AMP::pout << "L2 Norm of rInternal leaving DirichletVectorCorrection::apply is : "
                  << rInternal->L2Norm() << std::endl;
    }
}

void DirichletVectorCorrection::applyZeroValues( AMP::LinearAlgebra::Vector::shared_ptr r )
{
    auto rInternal = mySubsetVector( r, d_variable );
    auto dof_map   = rInternal->getDOFManager();
    size_t numIds  = d_boundaryIds.size();
    for ( size_t j = 0; j < numIds; j++ ) {
        auto bnd =
            d_Mesh->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, d_boundaryIds[j], 0 );
        auto end_bnd = bnd.end();

        for ( ; bnd != end_bnd; ++bnd ) {
            std::vector<size_t> bndGlobalIds;
            dof_map->getDOFs( bnd->globalID(), bndGlobalIds );
            const double val = 0.0;
            for ( auto &elem : d_dofIds[j] ) {
                rInternal->setLocalValuesByGlobalID( 1, &bndGlobalIds[elem], &val );
            } // end for i
        }     // end for bnd
    }         // end for j
    rInternal->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
}

void DirichletVectorCorrection::applyNonZeroValues( AMP::LinearAlgebra::Vector::shared_ptr r )
{
    AMP::LinearAlgebra::Vector::shared_ptr rInternal         = mySubsetVector( r, d_variable );
    std::shared_ptr<AMP::Discretization::DOFManager> dof_map = rInternal->getDOFManager();
    size_t numIds                                            = d_boundaryIds.size();
    for ( size_t j = 0; j < numIds; j++ ) {
        auto bnd =
            d_Mesh->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, d_boundaryIds[j], 0 );
        auto end_bnd = bnd.end();

        for ( ; bnd != end_bnd; ++bnd ) {
            std::vector<size_t> bndGlobalIds;
            dof_map->getDOFs( bnd->globalID(), bndGlobalIds );

            for ( size_t i = 0; i < d_dofIds[j].size(); ++i ) {
                double dVal =
                    ( d_valuesType == 1 ) ?
                        d_dirichletValues1[j][i] :
                        d_dirichletValues2->getLocalValueByGlobalID( bndGlobalIds[d_dofIds[j][i]] );
                dVal *= d_scalingFactor;
                rInternal->setLocalValuesByGlobalID( 1, &bndGlobalIds[d_dofIds[j][i]], &dVal );
            } // end for i
        }     // end for bnd
    }         // end for j
    rInternal->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
}

void DirichletVectorCorrection::applyResidual( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                                               AMP::LinearAlgebra::Vector::shared_ptr r )
{
    auto uInternal = mySubsetVector( u, d_variable );
    auto dof_map   = uInternal->getDOFManager();
    size_t numIds  = d_boundaryIds.size();
    for ( size_t j = 0; j < numIds; j++ ) {
        auto bnd =
            d_Mesh->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, d_boundaryIds[j], 0 );
        auto end_bnd = bnd.end();

        for ( ; bnd != end_bnd; ++bnd ) {
            std::vector<size_t> bndGlobalIds;
            dof_map->getDOFs( bnd->globalID(), bndGlobalIds );
            for ( size_t i = 0; i < d_dofIds[j].size(); i++ ) {
                double uVal = uInternal->getLocalValueByGlobalID( bndGlobalIds[d_dofIds[j][i]] );
                double dVal =
                    ( d_valuesType == 1 ) ?
                        d_dirichletValues1[j][i] :
                        d_dirichletValues2->getLocalValueByGlobalID( bndGlobalIds[d_dofIds[j][i]] );
                dVal = d_scalingFactor * ( uVal - dVal );
                r->setLocalValuesByGlobalID( 1, &bndGlobalIds[d_dofIds[j][i]], &dVal );
            } // end for i
        }     // end for bnd
    }         // end for j
    r->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
}

std::shared_ptr<OperatorParameters>
    DirichletVectorCorrection::getJacobianParameters( AMP::LinearAlgebra::Vector::const_shared_ptr )
{
    auto tmp_db = std::make_shared<AMP::Database>( "Dummy" );
    tmp_db->putScalar( "skip_params", true );

    auto outParams = std::make_shared<DirichletMatrixCorrectionParameters>( tmp_db );

    return outParams;
}

AMP::LinearAlgebra::Vector::shared_ptr
DirichletVectorCorrection::mySubsetVector( AMP::LinearAlgebra::Vector::shared_ptr vec,
                                           AMP::LinearAlgebra::Variable::shared_ptr var )
{
    if ( d_Mesh ) {
        AMP::LinearAlgebra::VS_Mesh meshSelector( d_Mesh );
        auto meshSubsetVec = vec->select( meshSelector, ( vec->getVariable() )->getName() );
        auto varSubsetVec  = meshSubsetVec->subsetVectorForVariable( var );
        return varSubsetVec;
    } else {
        return vec->subsetVectorForVariable( var );
    }
}

AMP::LinearAlgebra::Vector::const_shared_ptr
DirichletVectorCorrection::mySubsetVector( AMP::LinearAlgebra::Vector::const_shared_ptr vec,
                                           AMP::LinearAlgebra::Variable::shared_ptr var )
{
    if ( d_Mesh ) {
        AMP::LinearAlgebra::VS_Mesh meshSelector( d_Mesh );
        auto meshSubsetVec = vec->constSelect( meshSelector, ( vec->getVariable() )->getName() );
        auto varSubsetVec  = meshSubsetVec->constSubsetVectorForVariable( var );
        return varSubsetVec;
    } else {
        return vec->constSubsetVectorForVariable( var );
    }
}
} // namespace Operator
} // namespace AMP
