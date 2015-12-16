
#include "LinearFEOperator.h"
#include "ProfilerApp.h"
#include "matrices/MatrixBuilder.h"
#include "utils/Utilities.h"
#include "vectors/VectorBuilder.h"

#include "libmesh/cell_hex8.h"
#include "libmesh/node.h"

namespace AMP {
namespace Operator {


LinearFEOperator::LinearFEOperator( const AMP::shared_ptr<LinearFEOperatorParameters> &params )
    : LinearOperator( params ), d_currElemPtr( nullptr )
{
    d_elemOp    = ( params->d_elemOp );
    d_inDofMap  = ( params->d_inDofMap );
    d_outDofMap = ( params->d_outDofMap );
    if ( d_inDofMap == nullptr ) {
        d_inDofMap = d_outDofMap;
    }
    if ( d_outDofMap == nullptr ) {
        d_outDofMap = d_inDofMap;
    }
}


void LinearFEOperator::reset( const AMP::shared_ptr<OperatorParameters> &params )
{
    PROFILE_START( "reset" );
    AMP_INSIST( ( params != nullptr ), "NULL parameter" );
    AMP_INSIST( ( ( params->d_db ) != nullptr ), "NULL database" );

    const bool reuse_matrix = ( params->d_db )->getBoolWithDefault( "reset_reuses_matrix", true );

    if ( ( d_matrix.get() == nullptr ) || ( !reuse_matrix ) ) {
        AMP::LinearAlgebra::Vector::shared_ptr inVec =
            AMP::LinearAlgebra::createVector( d_inDofMap, getInputVariable(), true );
        AMP::LinearAlgebra::Vector::shared_ptr outVec =
            AMP::LinearAlgebra::createVector( d_outDofMap, getOutputVariable(), true );
        d_matrix = AMP::LinearAlgebra::createMatrix( inVec, outVec );
        d_matrix->zero();
        d_matrix->makeConsistent();
        d_matrix->makeConsistent();
    }
    AMP_ASSERT( ( *d_inDofMap ) == ( *d_matrix->getLeftDOFManager() ) );
    AMP_ASSERT( ( *d_inDofMap ) == ( *d_matrix->getRightDOFManager() ) );

    AMP::Mesh::MeshIterator el           = d_Mesh->getIterator( AMP::Mesh::Volume, 0 );
    const AMP::Mesh::MeshIterator end_el = el.end();

    this->preAssembly( params );

    for ( ; el != end_el; ++el ) {
        this->preElementOperation( *el );
        d_elemOp->apply();
        this->postElementOperation();
    } // end for el

    this->postAssembly();

    PROFILE_STOP( "reset" );
}


void LinearFEOperator::createCurrentLibMeshElement()
{
    d_currElemPtr = new ::Hex8;
    for ( size_t j = 0; j < d_currNodes.size(); ++j ) {
        std::vector<double> pt       = d_currNodes[j].coord();
        d_currElemPtr->set_node( j ) = new ::Node( pt[0], pt[1], pt[2], j );
    } // end for j
}


void LinearFEOperator::destroyCurrentLibMeshElement()
{
    for ( size_t j = 0; j < d_currElemPtr->n_nodes(); ++j ) {
        delete ( d_currElemPtr->get_node( j ) );
        d_currElemPtr->set_node( j ) = nullptr;
    } // end for j
    delete d_currElemPtr;
    d_currElemPtr = nullptr;
}
}
} // end namespace
