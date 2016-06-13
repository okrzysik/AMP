
#include "operators/map/Map3Dto1D.h"
#include "discretization/DOF_Manager.h"
#include "utils/InputDatabase.h"
#include "utils/Utilities.h"

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


namespace AMP {
namespace Operator {


Map3Dto1D::Map3Dto1D( const AMP::shared_ptr<OperatorParameters> &params ) : MapOperator( params )
{
    AMP::shared_ptr<MapOperatorParameters> myparams =
        AMP::dynamic_pointer_cast<MapOperatorParameters>( params );
    reset( myparams );
}


void Map3Dto1D::reset( const AMP::shared_ptr<OperatorParameters> &params )
{
    AMP::shared_ptr<MapOperatorParameters> myparams =
        AMP::dynamic_pointer_cast<MapOperatorParameters>( params );

    AMP_INSIST( ( ( myparams.get() ) != nullptr ), "NULL parameter" );
    AMP_INSIST( ( ( ( myparams->d_db ).get() ) != nullptr ), "NULL database" );
    AMP_INSIST( !myparams->d_MapComm.isNull(), "NULL communicator" );
    d_Mesh    = myparams->d_Mesh;
    d_MapMesh = myparams->d_MapMesh;
    d_MapComm = myparams->d_MapComm;
    AMP_INSIST( d_MapComm.sumReduce<int>( d_MapMesh.get() != nullptr ? 1 : 0 ) > 0,
                "Somebody must own the mesh" );

    d_useGaussVec = myparams->d_db->getBoolWithDefault( "UseGaussVec", false );

    if ( d_useGaussVec ) {
        AMP::Mesh::MeshIterator iterator =
            d_MapMesh->getBoundaryIDIterator( AMP::Mesh::Face, d_boundaryId, 0 );
        libmeshElements.reinit( iterator );
    }

    AMP_INSIST( myparams->d_db->keyExists( "InputVariable" ), "key not found" );
    std::string inpVar = myparams->d_db->getString( "InputVariable" );
    d_inpVariable.reset( new AMP::LinearAlgebra::Variable( inpVar ) );

    AMP_INSIST( myparams->d_db->keyExists( "OutputVariable" ), "key not found" );
    std::string outVar = myparams->d_db->getString( "OutputVariable" );
    d_outVariable.reset( new AMP::LinearAlgebra::Variable( outVar ) );
}


void Map3Dto1D::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                       AMP::LinearAlgebra::Vector::shared_ptr f )
{

    if ( d_useGaussVec ) {
        apply_Gauss( u, f );
    } else {
        apply_Nodal( u, f );
    }
}

void Map3Dto1D::apply_Gauss( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                             AMP::LinearAlgebra::Vector::shared_ptr )
{
    const unsigned int numPoints = outputVec->getLocalSize();
    std::vector<double> mapValues( numPoints, 0 );
    std::vector<int> numFaceGauss( numPoints, 0 );

    // Get the local contributions to the map
    if ( d_MapMesh != nullptr ) {
        AMP_ASSERT( u != nullptr );
        AMP::LinearAlgebra::Vector::const_shared_ptr inputVec = subsetInputVector( u );
        AMP_ASSERT( inputVec != nullptr );
        AMP_ASSERT( inputVec->getUpdateStatus() == AMP::LinearAlgebra::Vector::UNCHANGED );

        AMP::Discretization::DOFManager::shared_ptr dof_map = inputVec->getDOFManager();

        if ( d_iDebugPrintInfoLevel > 5 ) {
            AMP::pout << "The input to Map3Dto1D " << std::endl;
            AMP::pout << inputVec << std::endl;
        }

        AMP_ASSERT( outputVec != nullptr );

        // Get an iterator over the side elements
        AMP::Mesh::MeshIterator bnd =
            d_MapMesh->getBoundaryIDIterator( AMP::Mesh::Face, d_boundaryId, 0 );
        AMP::Mesh::MeshIterator end_bnd = bnd.end();

        // Iterator for the solid-clad boundary
        for ( ; bnd != end_bnd; ++bnd ) {
            libMeshEnums::Order feTypeOrder =
                Utility::string_to_enum<libMeshEnums::Order>( "FIRST" );
            libMeshEnums::FEFamily feFamily =
                Utility::string_to_enum<libMeshEnums::FEFamily>( "LAGRANGE" );

            AMP::shared_ptr<::FEType> d_feType( new ::FEType( feTypeOrder, feFamily ) );
            AMP::shared_ptr<::FEBase> d_fe( (::FEBase::build( 2, ( *d_feType ) ) ).release() );

            libMeshEnums::Order qruleOrder =
                Utility::string_to_enum<libMeshEnums::Order>( "SECOND" );
            AMP::shared_ptr<::QBase> d_qrule(
                (::QBase::build( "QGAUSS", 2, qruleOrder ) ).release() );

            d_fe->attach_quadrature_rule( d_qrule.get() );

            d_fe->reinit( libmeshElements.getElement( bnd->globalID() ) );

            // Get the current position and DOF
            std::vector<Point> coordinates = d_fe->get_xyz();

            std::vector<size_t> ids;
            dof_map->getDOFs( bnd->globalID(), ids );

            std::vector<double> zcoords;
            for ( size_t i = 0; i < ids.size(); i++ ) {
                zcoords.push_back( coordinates[i]( 2 ) );
            }

            std::sort( zcoords.begin(), zcoords.end() );

            std::vector<double> tmpZcoords = zcoords;
            std::vector<int> tmpIds( zcoords.size() );

            for ( size_t i = 0; i < ids.size(); i++ ) {
                tmpIds[i] = i;
            }

            std::vector<int> originalGaussOrder( zcoords.size() );

            for ( size_t i = 0; i < ids.size(); i++ ) {
                double myZ = coordinates[i]( 2 );
                for ( unsigned int j = 0; j < tmpZcoords.size(); j++ ) {
                    if ( fabs( tmpZcoords[j] - myZ ) <= 1.e-12 ) {
                        originalGaussOrder[tmpIds[j]] = i;
                        tmpZcoords.erase( tmpZcoords.begin() + j );
                        tmpIds.erase( tmpIds.begin() + j );
                        break;
                    }
                }
            }

            std::vector<double> z( 4, 0 );
            std::vector<double> y( 4, 0 );
            std::vector<double> x( 4, 0 );
            for ( int i = 0; i < 4; i++ ) {
                x[i] = coordinates[originalGaussOrder[i]]( 0 );
                y[i] = coordinates[originalGaussOrder[i]]( 1 );
                z[i] = coordinates[originalGaussOrder[i]]( 2 );
            }

            int pickId;
            if ( pow( ( y[0] - y[3] ), 2 ) + pow( ( x[0] - x[3] ), 2 ) <
                 pow( ( y[0] - y[2] ), 2 ) + pow( ( x[0] - x[2] ), 2 ) ) {
                pickId = 3;
            } else {
                pickId = 2;
            }

            // Iterator for the fluid boundary
            for ( unsigned int i = 0; i < numPoints; i++ ) {

                double cur_gauss = d_zLocations[i];

                // Section of the Clad boundary map corresponding to the fluid Element
                if ( cur_gauss >= z[0] && cur_gauss <= z[pickId] ) {
                    mapValues[i] +=
                        ( ( inputVec )->getValueByGlobalID( ids[originalGaussOrder[0]] ) *
                              ( z[pickId] - cur_gauss ) +
                          ( inputVec )->getValueByGlobalID( ids[originalGaussOrder[pickId]] ) *
                              ( cur_gauss - z[0] ) ) /
                        ( z[pickId] - z[0] );
                    numFaceGauss[i] += 1;
                }
            } // end for i
        }     // end for bnd
    }         // end if

    // Gather the results from all processors
    std::vector<double> aggMapValues( numPoints );
    std::vector<int> aggNumFaceGauss( numPoints );
    d_MapComm.sumReduce( (double *) &( mapValues[0] ), (double *) &( aggMapValues[0] ), numPoints );
    d_MapComm.sumReduce( (int *) &( numFaceGauss[0] ), (int *) &( aggNumFaceGauss[0] ), numPoints );

    // Store the results
    for ( unsigned int i = 0; i < numPoints; i++ ) {
        outputVec->setValueByLocalID( i,
                                      aggMapValues[i] / static_cast<double>( aggNumFaceGauss[i] ) );
    }

    if ( d_iDebugPrintInfoLevel > 4 ) {
        AMP::pout << "The output to Map3Dto1D " << std::endl;
        AMP::pout << outputVec << std::endl;
    }
}


void Map3Dto1D::apply_Nodal( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                             AMP::LinearAlgebra::Vector::shared_ptr )
{

    const unsigned int numPoints = outputVec->getLocalSize();
    std::vector<double> mapValues( numPoints, 0 );
    std::vector<int> numFaceNodes( numPoints, 0 );

    // Get the local contributions to the map
    if ( d_MapMesh != nullptr ) {
        AMP_ASSERT( u != nullptr );

        // Subset u for the local vector of interest
        AMP::LinearAlgebra::Vector::const_shared_ptr inputVec = subsetInputVector( u );
        AMP_ASSERT( inputVec != nullptr );
        AMP_ASSERT( inputVec->getUpdateStatus() == AMP::LinearAlgebra::Vector::UNCHANGED );

        AMP::Discretization::DOFManager::shared_ptr dof_map = inputVec->getDOFManager();

        if ( d_iDebugPrintInfoLevel > 5 ) {
            AMP::pout << "The input to Map3Dto1D " << std::endl;
            AMP::pout << inputVec << std::endl;
        }

        AMP_ASSERT( outputVec != nullptr );

        // Get an iterator over the side elements
        AMP::Mesh::MeshIterator bnd =
            d_MapMesh->getBoundaryIDIterator( AMP::Mesh::Face, d_boundaryId, 0 );
        AMP::Mesh::MeshIterator end_bnd = bnd.end();

        // Iterator for the solid-clad boundary
        for ( ; bnd != end_bnd; ++bnd ) {

            AMP::Mesh::MeshElement cur_side           = *bnd;
            std::vector<AMP::Mesh::MeshElement> nodes = cur_side.getElements( AMP::Mesh::Vertex );
            AMP_ASSERT( nodes.size() == 4 );

            std::vector<double> zcoords;
            for ( auto &node : nodes ) {
                std::vector<double> coord = node.coord();
                zcoords.push_back( coord[2] );
            }

            std::sort( zcoords.begin(), zcoords.end() );

            std::vector<double> tmpZcoords = zcoords;
            std::vector<int> tmpIds( zcoords.size() );

            for ( size_t i = 0; i < nodes.size(); i++ ) {
                tmpIds[i] = i;
            }

            std::vector<int> originalNodeOrder( zcoords.size() );

            for ( size_t i = 0; i < nodes.size(); i++ ) {
                std::vector<double> coord = nodes[i].coord();
                double myZ                = coord[2];
                for ( unsigned int j = 0; j < tmpZcoords.size(); j++ ) {
                    if ( fabs( tmpZcoords[j] - myZ ) <= 1.e-12 ) {
                        originalNodeOrder[tmpIds[j]] = i;
                        tmpZcoords.erase( tmpZcoords.begin() + j );
                        tmpIds.erase( tmpIds.begin() + j );
                        break;
                    }
                }
            }

            std::vector<double> z( 4, 0 );
            std::vector<double> y( 4, 0 );
            std::vector<double> x( 4, 0 );
            for ( int i = 0; i < 4; i++ ) {
                std::vector<double> coord = nodes[originalNodeOrder[i]].coord();
                x[i]                      = coord[0];
                y[i]                      = coord[1];
                z[i]                      = coord[2];
            }

            int pickId;
            if ( pow( ( y[0] - y[3] ), 2 ) + pow( ( x[0] - x[3] ), 2 ) <
                 pow( ( y[0] - y[2] ), 2 ) + pow( ( x[0] - x[2] ), 2 ) ) {
                pickId = 3;
            } else {
                pickId = 2;
            }

            // Iterator for the fluid boundary
            for ( unsigned int i = 0; i < numPoints; i++ ) {

                double cur_node = d_zLocations[i];

                // Section of the Clad boundary map corresponding to the fluid Element
                if ( cur_node >= z[0] && cur_node <= z[pickId] ) {
                    std::vector<size_t> dof1;
                    std::vector<size_t> dof2;
                    dof_map->getDOFs( nodes[originalNodeOrder[0]].globalID(), dof1 );
                    dof_map->getDOFs( nodes[originalNodeOrder[pickId]].globalID(), dof2 );
                    AMP_ASSERT( dof1.size() == 1 && dof2.size() == 1 );

                    mapValues[i] +=
                        ( ( inputVec )->getValueByGlobalID( dof1[0] ) * ( z[pickId] - cur_node ) +
                          ( inputVec )->getValueByGlobalID( dof2[0] ) * ( cur_node - z[0] ) ) /
                        ( z[pickId] - z[0] );

                    numFaceNodes[i] += 1;
                }
            } // end for i
        }     // end for bnd
    }

    // Gather the results from all processors
    std::vector<double> aggMapValues( numPoints );
    std::vector<int> aggNumFaceNodes( numPoints );
    d_MapComm.sumReduce( (double *) &( mapValues[0] ), (double *) &( aggMapValues[0] ), numPoints );
    d_MapComm.sumReduce( (int *) &( numFaceNodes[0] ), (int *) &( aggNumFaceNodes[0] ), numPoints );

    // Store the results
    for ( unsigned int i = 0; i < numPoints; i++ ) {
        outputVec->setValueByLocalID( i,
                                      aggMapValues[i] / static_cast<double>( aggNumFaceNodes[i] ) );
    }

    if ( d_iDebugPrintInfoLevel > 4 ) {
        AMP::pout << "The output to Map3Dto1D " << std::endl;
        AMP::pout << outputVec << std::endl;
    }
}
}
}
