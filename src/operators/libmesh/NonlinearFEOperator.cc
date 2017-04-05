
#include "NonlinearFEOperator.h"
#include "ProfilerApp.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/node.h"
#include "utils/Utilities.h"

namespace AMP {
namespace Operator {


NonlinearFEOperator::NonlinearFEOperator( const AMP::shared_ptr<FEOperatorParameters> &params )
    : Operator( params )
{
    d_elemOp = ( params->d_elemOp );
    createLibMeshElementList();
    d_currElemIdx = static_cast<unsigned int>( -1 );
}


NonlinearFEOperator::~NonlinearFEOperator() { destroyLibMeshElementList(); }


void NonlinearFEOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                                 AMP::LinearAlgebra::Vector::shared_ptr r )
{
    PROFILE_START( "apply" );

    AMP_INSIST( ( r != nullptr ), "NULL Residual/Output Vector" );
    AMP::LinearAlgebra::Vector::shared_ptr rInternal = this->subsetOutputVector( r );
    AMP_INSIST( ( rInternal != nullptr ), "NULL Residual/Output Vector" );

    if ( u.get() != nullptr )
        AMP_ASSERT( u->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::UNCHANGED );

    d_currElemIdx = static_cast<unsigned int>( -1 );
    this->preAssembly( u, rInternal );

    PROFILE_START( "loop over elements" );
    AMP::Mesh::MeshIterator el = d_Mesh->getIterator( AMP::Mesh::GeomType::Volume, 0 );
    for ( d_currElemIdx = 0; d_currElemIdx < el.size(); ++d_currElemIdx, ++el ) {
        this->preElementOperation( *el );
        d_elemOp->apply();
        this->postElementOperation();
    } // end for el
    PROFILE_STOP( "loop over elements" );

    d_currElemIdx = static_cast<unsigned int>( -1 );
    this->postAssembly();

    rInternal->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );

    if ( d_iDebugPrintInfoLevel > 2 )
        AMP::pout << "L2 norm of result of NonlinearFEOperator::apply is: " << rInternal->L2Norm()
                  << std::endl;
    if ( d_iDebugPrintInfoLevel > 5 )
        std::cout << rInternal << std::endl;

    PROFILE_STOP( "apply" );
}


void NonlinearFEOperator::createLibMeshElementList()
{
    AMP::Mesh::MeshIterator el = d_Mesh->getIterator( AMP::Mesh::GeomType::Volume, 0 );
    d_currElemPtrs.resize( d_Mesh->numLocalElements( AMP::Mesh::GeomType::Volume ) );
    for ( size_t i = 0; i < el.size(); ++el, ++i ) {
        std::vector<AMP::Mesh::MeshElement> currNodes = el->getElements( AMP::Mesh::GeomType::Vertex );
        d_currElemPtrs[i]                             = new ::Hex8;
        for ( size_t j = 0; j < currNodes.size(); ++j ) {
            std::vector<double> pt           = currNodes[j].coord();
            d_currElemPtrs[i]->set_node( j ) = new ::Node( pt[0], pt[1], pt[2], j );
        } // end for j
    }     // end for i
}


void NonlinearFEOperator::destroyLibMeshElementList()
{
    for ( auto &elem : d_currElemPtrs ) {
        for ( size_t j = 0; j < elem->n_nodes(); ++j ) {
            delete ( elem->get_node( j ) );
            elem->set_node( j ) = nullptr;
        } // end for j
        delete ( elem );
        elem = nullptr;
    } // end for i
    d_currElemPtrs.clear();
}
}
} // end namespace
