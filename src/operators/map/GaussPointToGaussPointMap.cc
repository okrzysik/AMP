
#include "operators/map/GaussPointToGaussPointMap.h"
#include "ampmesh/MultiMesh.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"

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

namespace AMP {
namespace Operator {


// Constructor
GaussPointToGaussPointMap::GaussPointToGaussPointMap(
    const AMP::shared_ptr<AMP::Operator::OperatorParameters> &params )
    : NodeToNodeMap( params )
{
    createIdxMap( params );
    d_useFrozenInputVec = params->d_db->getBoolWithDefault( "FrozenInput", false );
}


// Apply start
void GaussPointToGaussPointMap::applyStart( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                                            AMP::LinearAlgebra::Vector::shared_ptr f )
{
    AMP::LinearAlgebra::Vector::const_shared_ptr uInternal = u;
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
    AMP::Discretization::DOFManager::shared_ptr dofMap = d_OutputVector->getDOFManager();
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
                d_OutputVector->setLocalValueByGlobalID(
                    localDofs[( j * DofsPerGaussPt ) + k],
                    vals[( ( d_idxMap[i][j] ) * DofsPerGaussPt ) + k] );
            } // end k
        }     // end j
    }         // end i
}


void GaussPointToGaussPointMap::createIdxMap(
    AMP::shared_ptr<AMP::Operator::OperatorParameters> params )
{
    AMP::shared_ptr<AMP::Database> db = params->d_db;
    std::string feTypeOrderName       = db->getStringWithDefault( "FE_ORDER", "FIRST" );
    libMeshEnums::Order feTypeOrder =
        Utility::string_to_enum<libMeshEnums::Order>( feTypeOrderName );

    std::string feFamilyName = db->getStringWithDefault( "FE_FAMILY", "LAGRANGE" );
    libMeshEnums::FEFamily feFamily =
        Utility::string_to_enum<libMeshEnums::FEFamily>( feFamilyName );

    std::string qruleTypeName = db->getStringWithDefault( "QRULE_TYPE", "QGAUSS" );
    libMeshEnums::QuadratureType qruleType =
        Utility::string_to_enum<libMeshEnums::QuadratureType>( qruleTypeName );

    std::string qruleOrderName = db->getStringWithDefault( "QRULE_ORDER", "DEFAULT" );

    int faceDim = db->getIntegerWithDefault( "DIMENSION", 2 );

    AMP::shared_ptr<::FEType> feType( new ::FEType( feTypeOrder, feFamily ) );

    libMeshEnums::Order qruleOrder;

    if ( qruleOrderName == "DEFAULT" ) {
        qruleOrder = feType->default_quadrature_order();
    } else {
        qruleOrder = Utility::string_to_enum<libMeshEnums::Order>( qruleOrderName );
    }

    AMP::shared_ptr<::QBase> qrule( (::QBase::build( qruleType, faceDim, qruleOrder ) ).release() );
    qrule->init( QUAD4, 0 );

    unsigned int numGaussPtsPerElem = qrule->n_points();

    unsigned int dofsPerElem = ( dim * numGaussPtsPerElem );

    AMP::LinearAlgebra::Variable::shared_ptr variable(
        new AMP::LinearAlgebra::Variable( "GaussPoints" ) );

    std::vector<AMP::Mesh::Mesh::shared_ptr> meshesForMap( 2 );
    meshesForMap[0] = d_mesh1;
    meshesForMap[1] = d_mesh2;
    AMP::Mesh::Mesh::shared_ptr multiMesh( new AMP::Mesh::MultiMesh( d_MapComm, meshesForMap ) );

    //      AMP::Mesh::MeshIterator surfIter =
    //      multiMesh->getSurfaceIterator(AMP::Mesh::GeomType::Face, 0);
    //      AMP::Discretization::DOFManager::shared_ptr dofMap =
    //      AMP::Discretization::simpleDOFManager::create(multiMesh,
    //          surfIter, surfIter, dofsPerElem);
    AMP::Mesh::Mesh::shared_ptr submesh =
        multiMesh->Subset( multiMesh->getSurfaceIterator( AMP::Mesh::GeomType::Face, 0 ) );
    AMP::Discretization::DOFManager::shared_ptr dofMap =
        AMP::Discretization::simpleDOFManager::create(
            submesh, AMP::Mesh::GeomType::Face, 0, dofsPerElem, true );

    AMP::LinearAlgebra::Vector::shared_ptr inVec =
        AMP::LinearAlgebra::createVector( dofMap, variable );

    AMP::LinearAlgebra::Vector::shared_ptr outVec = inVec->cloneVector();

    std::vector<size_t> localDofs( dofsPerElem );
    for ( auto &_i : d_sendList ) {
        AMP::Mesh::MeshElement el = multiMesh->getElement( _i );

        std::vector<AMP::Mesh::MeshElement> currNodes =
            el.getElements( AMP::Mesh::GeomType::Vertex );

        ::Elem *elem = new ::Quad4;
        for ( size_t j = 0; j < currNodes.size(); ++j ) {
            std::vector<double> pt = currNodes[j].coord();
            elem->set_node( j )    = new ::Node( pt[0], pt[1], pt[2], j );
        } // end for j

        AMP::shared_ptr<::FEBase> fe( (::FEBase::build( faceDim, ( *feType ) ) ).release() );
        fe->attach_quadrature_rule( qrule.get() );
        fe->reinit( elem );

        const std::vector<::Point> &xyz = fe->get_xyz();

        dofMap->getDOFs( _i, localDofs );

        for ( unsigned int j = 0; j < numGaussPtsPerElem; ++j ) {
            for ( int k = 0; k < dim; ++k ) {
                inVec->setLocalValueByGlobalID( localDofs[( j * dim ) + k], xyz[j]( k ) );
            } // end for k
        }     // end for j

        for ( unsigned int j = 0; j < elem->n_nodes(); ++j ) {
            delete ( elem->get_node( j ) );
            elem->set_node( j ) = nullptr;
        } // end for j
        delete elem;
        elem = nullptr;
    } // end i

    db->putInteger( "DOFsPerObject", dofsPerElem );
    db->putString( "VariableName", "GaussPoints" );
    AMP::shared_ptr<AMP::Operator::NodeToNodeMap> n2nMap(
        new AMP::Operator::NodeToNodeMap( params ) );
    n2nMap->setVector( outVec );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    n2nMap->apply( inVec, nullVec );

    d_idxMap.clear();
    for ( auto &_i : d_recvList ) {
        AMP::Mesh::MeshElement el = multiMesh->getElement( _i );

        std::vector<AMP::Mesh::MeshElement> currNodes =
            el.getElements( AMP::Mesh::GeomType::Vertex );

        ::Elem *elem = new ::Quad4;
        for ( size_t j = 0; j < currNodes.size(); ++j ) {
            std::vector<double> pt = currNodes[j].coord();
            elem->set_node( j )    = new ::Node( pt[0], pt[1], pt[2], j );
        } // end for j

        AMP::shared_ptr<::FEBase> fe( (::FEBase::build( faceDim, ( *feType ) ) ).release() );
        fe->attach_quadrature_rule( qrule.get() );
        fe->reinit( elem );

        const std::vector<::Point> &xyz = fe->get_xyz();

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
            delete ( elem->get_node( j ) );
            elem->set_node( j ) = nullptr;
        } // end for j
        delete elem;
        elem = nullptr;
    } // end for i
}
}
}
