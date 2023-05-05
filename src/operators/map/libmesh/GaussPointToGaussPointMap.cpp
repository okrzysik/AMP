#include "AMP/operators/map/libmesh/GaussPointToGaussPointMap.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/MultiMesh.h"
#include "AMP/vectors/VectorBuilder.h"

// Libmesh files
DISABLE_WARNINGS
#include "libmesh/auto_ptr.h"
#include "libmesh/elem.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/face_quad4.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_type.h"
#include "libmesh/node.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
ENABLE_WARNINGS

#include <string>

namespace AMP::Operator {


// Constructor
GaussPointToGaussPointMap::GaussPointToGaussPointMap(
    std::shared_ptr<const AMP::Operator::OperatorParameters> params )
    : NodeToNodeMap( params )
{
    createIdxMap( params );
    d_useFrozenInputVec = params->d_db->getWithDefault<bool>( "FrozenInput", false );
}


// Apply start
void GaussPointToGaussPointMap::applyStart( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                                            AMP::LinearAlgebra::Vector::shared_ptr f )
{
    auto uInternal = u;
    if ( d_useFrozenInputVec ) {
        uInternal = d_frozenInputVec;
    }
    AMP::Operator::NodeToNodeMap::applyStart( uInternal, f );
}


// Apply finish
void GaussPointToGaussPointMap::applyFinish( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                                             AMP::LinearAlgebra::Vector::shared_ptr f )
{
    AMP::Operator::NodeToNodeMap::applyFinish( u, f );
    correctLocalOrdering();
}


// Check if we have the correct map
bool GaussPointToGaussPointMap::validMapType( const std::string &t )
{
    if ( t == "GaussPointToGaussPoint" )
        return true;
    return false;
}


// Correct the local ordering
void GaussPointToGaussPointMap::correctLocalOrdering()
{
    auto dofMap = d_OutputVector->getDOFManager();
    std::vector<size_t> localDofs( DofsPerObj );
    for ( size_t i = 0; i < d_recvList.size(); ++i ) {
        dofMap->getDOFs( d_recvList[i], localDofs );
        std::vector<double> vals( DofsPerObj );
        for ( int j = 0; j < DofsPerObj; ++j ) {
            vals[j] = d_OutputVector->getLocalValueByGlobalID( localDofs[j] );
        } // end j
        int DofsPerGaussPt = DofsPerObj / ( d_idxMap[i].size() );
        for ( size_t j = 0; j < d_idxMap[i].size(); ++j ) {
            for ( int k = 0; k < DofsPerGaussPt; ++k ) {
                d_OutputVector->setLocalValuesByGlobalID(
                    1,
                    &localDofs[( j * DofsPerGaussPt ) + k],
                    &vals[( ( d_idxMap[i][j] ) * DofsPerGaussPt ) + k] );
            } // end k
        }     // end j
    }         // end i
}


void GaussPointToGaussPointMap::createIdxMap(
    std::shared_ptr<const AMP::Operator::OperatorParameters> params )
{
    auto db              = params->d_db;
    auto feTypeOrderName = db->getWithDefault<std::string>( "FE_ORDER", "FIRST" );
    auto feTypeOrder     = libMesh::Utility::string_to_enum<libMeshEnums::Order>( feTypeOrderName );

    auto feFamilyName = db->getWithDefault<std::string>( "FE_FAMILY", "LAGRANGE" );
    auto feFamily     = libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( feFamilyName );

    auto qruleTypeName = db->getWithDefault<std::string>( "QRULE_TYPE", "QGAUSS" );
    auto qruleType =
        libMesh::Utility::string_to_enum<libMeshEnums::QuadratureType>( qruleTypeName );

    auto qruleOrderName = db->getWithDefault<std::string>( "QRULE_ORDER", "DEFAULT" );

    int faceDim = db->getWithDefault<int>( "DIMENSION", 2 );

    auto feType = std::make_shared<libMesh::FEType>( feTypeOrder, feFamily );

    libMeshEnums::Order qruleOrder;

    if ( qruleOrderName == "DEFAULT" ) {
        qruleOrder = feType->default_quadrature_order();
    } else {
        qruleOrder = libMesh::Utility::string_to_enum<libMeshEnums::Order>( qruleOrderName );
    }

    std::shared_ptr<libMesh::QBase> qrule(
        libMesh::QBase::build( qruleType, faceDim, qruleOrder ).release() );
    qrule->init( libMesh::QUAD4, 0 );

    unsigned int numGaussPtsPerElem = qrule->n_points();

    unsigned int dofsPerElem = ( dim * numGaussPtsPerElem );

    auto variable = std::make_shared<AMP::LinearAlgebra::Variable>( "GaussPoints" );

    std::vector<std::shared_ptr<AMP::Mesh::Mesh>> meshesForMap( 2 );
    meshesForMap[0] = d_mesh1;
    meshesForMap[1] = d_mesh2;
    auto multiMesh = std::make_shared<AMP::Mesh::MultiMesh>( "MultiMesh", d_MapComm, meshesForMap );

    // auto surfIter = multiMesh->getSurfaceIterator(AMP::Mesh::GeomType::Face, 0);
    // auto dofMap = AMP::Discretization::simpleDOFManager::create(multiMesh, surfIter, surfIter,
    // dofsPerElem);
    auto submesh =
        multiMesh->Subset( multiMesh->getSurfaceIterator( AMP::Mesh::GeomType::Face, 0 ) );
    auto dofMap = AMP::Discretization::simpleDOFManager::create(
        submesh, AMP::Mesh::GeomType::Face, 0, dofsPerElem, true );

    auto inVec  = AMP::LinearAlgebra::createVector( dofMap, variable );
    auto outVec = inVec->clone();

    std::vector<size_t> localDofs( dofsPerElem );
    for ( auto &_i : d_sendList ) {
        AMP::Mesh::MeshElement el = multiMesh->getElement( _i );

        auto currNodes = el.getElements( AMP::Mesh::GeomType::Vertex );

        libMesh::Elem *elem = new libMesh::Quad4;
        for ( size_t j = 0; j < currNodes.size(); ++j ) {
            auto pt             = currNodes[j].coord();
            elem->set_node( j ) = new libMesh::Node( pt[0], pt[1], pt[2], j );
        } // end for j

        std::shared_ptr<libMesh::FEBase> fe(
            ( libMesh::FEBase::build( faceDim, ( *feType ) ) ).release() );
        fe->attach_quadrature_rule( qrule.get() );
        fe->reinit( elem );

        const auto &xyz = fe->get_xyz();

        dofMap->getDOFs( _i, localDofs );

        for ( unsigned int j = 0; j < numGaussPtsPerElem; ++j ) {
            for ( int k = 0; k < dim; ++k ) {
                inVec->setLocalValuesByGlobalID( 1, &localDofs[( j * dim ) + k], &xyz[j]( k ) );
            } // end for k
        }     // end for j

        for ( unsigned int j = 0; j < elem->n_nodes(); ++j ) {
            delete ( elem->node_ptr( j ) );
            elem->set_node( j ) = nullptr;
        } // end for j
        delete elem;
        elem = nullptr;
    } // end i

    db->putScalar( "DOFsPerObject", dofsPerElem );
    db->putScalar( "VariableName", "GaussPoints" );
    auto n2nMap = std::make_shared<AMP::Operator::NodeToNodeMap>( params );
    n2nMap->setVector( outVec );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    n2nMap->apply( inVec, nullVec );

    d_idxMap.clear();
    for ( auto &_i : d_recvList ) {
        auto el = multiMesh->getElement( _i );

        auto currNodes = el.getElements( AMP::Mesh::GeomType::Vertex );

        libMesh::Elem *elem = new libMesh::Quad4;
        for ( size_t j = 0; j < currNodes.size(); ++j ) {
            auto pt             = currNodes[j].coord();
            elem->set_node( j ) = new libMesh::Node( pt[0], pt[1], pt[2], j );
        } // end for j

        std::shared_ptr<libMesh::FEBase> fe(
            ( libMesh::FEBase::build( faceDim, ( *feType ) ) ).release() );
        fe->get_xyz();
        fe->attach_quadrature_rule( qrule.get() );
        fe->reinit( elem );

        const auto &xyz = fe->get_xyz();

        dofMap->getDOFs( _i, localDofs );

        std::vector<double> vals( dofsPerElem );
        for ( unsigned int j = 0; j < dofsPerElem; ++j ) {
            vals[j] = outVec->getLocalValueByGlobalID( localDofs[j] );
        } // end j

        std::vector<unsigned int> locMap( numGaussPtsPerElem, static_cast<unsigned int>( -1 ) );

        for ( unsigned int j = 0; j < numGaussPtsPerElem; ++j ) {
            for ( unsigned int k = 0; k < numGaussPtsPerElem; ++k ) {
                bool found = true;
                for ( int d = 0; d < dim; ++d ) {
                    if ( fabs( xyz[j]( d ) - vals[( k * dim ) + d] ) > 1.0e-11 ) {
                        found = false;
                        break;
                    }
                } // end d
                if ( found ) {
                    locMap[j] = k;
                    break;
                }
            } // end k
        }     // end j

        d_idxMap.push_back( locMap );

        for ( unsigned int j = 0; j < elem->n_nodes(); ++j ) {
            delete ( elem->node_ptr( j ) );
            elem->set_node( j ) = nullptr;
        } // end for j
        delete elem;
        elem = nullptr;
    } // end for i
}
} // namespace AMP::Operator
