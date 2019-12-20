#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/libmesh/libmeshMesh.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/VectorBuilder.h"
#include <memory>

#include "libmesh/auto_ptr.h"
#include "libmesh/boundary_info.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/elem.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_type.h"
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/node.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>


static void myTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string log_file = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    std::shared_ptr<libMesh::Mesh> mesh( new libMesh::Mesh(libMesh::Parallel::Communicator(), 3 ) );
    libMesh::MeshTools::Generation::build_cube(
				      ( *( mesh.get() ) ), 1, 1, 1, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, libMesh::HEX8, false );

    libMesh::Elem *elemPtr = mesh->elem( 0 );

    ( elemPtr->point( 4 ) )( 0 ) -= 0.4;

    ( elemPtr->point( 5 ) )( 0 ) -= 0.4;

    ( elemPtr->point( 6 ) )( 0 ) -= 0.4;

    ( elemPtr->point( 7 ) )( 0 ) -= 0.4;

    auto ampMesh = AMP::Mesh::Mesh::shared_ptr( new AMP::Mesh::libmeshMesh( mesh, "TestMesh" ) );

    auto myVar   = std::make_shared<AMP::LinearAlgebra::Variable>( "myVar" );
    auto dof_map = AMP::Discretization::simpleDOFManager::create(
        ampMesh, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    auto T = AMP::LinearAlgebra::createVector( dof_map, myVar, true );

    FILE *fp = fopen( "InverseJacobian.txt", "w" );

    auto el     = ampMesh->getIterator( AMP::Mesh::GeomType::Volume, 0 );
    auto end_el = el.end();

    while ( el != end_el ) {
        auto nodes = el->getElements( AMP::Mesh::GeomType::Vertex );
        for ( size_t i = 0; i < nodes.size(); i++ ) {
            auto pt = nodes[i].coord();
            fprintf(
                fp, "nd = %d, x = %.15f, y = %.15f, z = %.15f \n", (int) i, pt[0], pt[1], pt[2] );
        }
    }

    auto feTypeOrder = libMesh::Utility::string_to_enum<libMeshEnums::Order>( "FIRST" );
    auto feFamily    = libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( "LAGRANGE" );

    std::shared_ptr<libMesh::FEType> feType( new libMesh::FEType( feTypeOrder, feFamily ) );
    std::shared_ptr<libMesh::FEBase> fe( (libMesh::FEBase::build( 3, ( *feType ) ) ).release() );

    //  const std::vector<std::vector<libMesh::Real> > &phi = fe->get_phi();
    const std::vector<std::vector<libMesh::Real>> &dphidxi   = fe->get_dphidxi();
    const std::vector<std::vector<libMesh::Real>> &dphideta  = fe->get_dphideta();
    const std::vector<std::vector<libMesh::Real>> &dphidzeta = fe->get_dphidzeta();
    const std::vector<std::vector<libMesh::Real>> &dphidx    = fe->get_dphidx();
    const std::vector<std::vector<libMesh::Real>> &dphidy    = fe->get_dphidy();
    const std::vector<std::vector<libMesh::Real>> &dphidz    = fe->get_dphidz();
    //  const std::vector<libMesh::Real> &djxw = fe->get_JxW();

    const std::vector<libMesh::RealGradient> &dxyzdxi   = fe->get_dxyzdxi();
    const std::vector<libMesh::RealGradient> &dxyzdeta  = fe->get_dxyzdeta();
    const std::vector<libMesh::RealGradient> &dxyzdzeta = fe->get_dxyzdzeta();

    auto qruleType = libMesh::Utility::string_to_enum<libMeshEnums::QuadratureType>( "QGAUSS" );
    libMeshEnums::Order qruleOrder = feType->default_quadrature_order();
    std::shared_ptr<libMesh::QBase> qrule( (libMesh::QBase::build( qruleType, 3, qruleOrder ) ).release() );

    fe->attach_quadrature_rule( qrule.get() );

    std::vector<AMP::Mesh::MeshElement> d_currNodes;

    std::vector<size_t> d_dofIndices;
    el          = ampMesh->getIterator( AMP::Mesh::GeomType::Volume, 0 );
    d_currNodes = el->getElements( AMP::Mesh::GeomType::Vertex );

    std::vector<AMP::Mesh::MeshElementID> globalIDs( d_currNodes.size() );

    for ( unsigned int j = 0; j < d_currNodes.size(); j++ ) {
        globalIDs[j] = d_currNodes[j].globalID();
    } // end of j
    dof_map->getDOFs( globalIDs, d_dofIndices );

    libMesh::Elem *currElemPtr;
    currElemPtr = new libMesh::Hex8;
    for ( size_t j = 0; j < d_currNodes.size(); j++ ) {
        auto pt                    = d_currNodes[j].coord();
        currElemPtr->set_node( j ) = new libMesh::Node( pt[0], pt[1], pt[2], j );
    } // end for j

    fe->reinit( currElemPtr );

    for ( unsigned int i = 0; i < d_dofIndices.size(); i++ ) {
        T->setValueByGlobalID( d_dofIndices[i], 300 * ( i + 1 ) );
    }

    fprintf( fp,
             " dx/dxi = %.15f, dydxi = %.15f, dzdxi = %.15f \n",
             dxyzdxi[0]( 0 ),
             dxyzdxi[0]( 1 ),
             dxyzdxi[0]( 2 ) );
    fprintf( fp,
             " dx/deta = %.15f, dydeta = %.15f, dzdeta = %.15f \n",
             dxyzdeta[0]( 0 ),
             dxyzdeta[0]( 1 ),
             dxyzdeta[0]( 2 ) );
    fprintf( fp,
             " dx/dzeta = %.15f, dydzeta = %.15f, dzdzeta = %.15f \n",
             dxyzdzeta[0]( 0 ),
             dxyzdzeta[0]( 1 ),
             dxyzdzeta[0]( 2 ) );

    std::vector<libMesh::Real> Jinv1( 3 );
    std::vector<libMesh::Real> Jinv2( 3 );
    std::vector<libMesh::Real> Jinv3( 3 );
    Jinv1[0] = 2.;
    Jinv1[1] = 0.;
    Jinv1[2] = 0.;

    Jinv2[0] = 0;
    Jinv2[1] = 2;
    Jinv2[2] = 0;

    Jinv3[0] = 0.8;
    Jinv3[1] = 0.;
    Jinv3[2] = 2.;

    libMesh::Real dTdxi = 0, dTdeta = 0, dTdzeta = 0, dTdx = 0;
    libMesh::Real dTdy = 0, dTdz = 0, lib_dTdx = 0, lib_dTdy = 0, lib_dTdz = 0;

    for ( unsigned int i = 0; i < d_dofIndices.size(); i++ ) {
        dTdxi += ( dphidxi[i][0] * T->getValueByGlobalID( d_dofIndices[i] ) );
        dTdeta += ( dphideta[i][0] * T->getValueByGlobalID( d_dofIndices[i] ) );
        dTdzeta += ( dphidzeta[i][0] * T->getValueByGlobalID( d_dofIndices[i] ) );
    }

    dTdx = Jinv1[0] * dTdxi + Jinv1[1] * dTdeta + Jinv1[2] * dTdzeta;
    dTdy = Jinv2[0] * dTdxi + Jinv2[1] * dTdeta + Jinv2[2] * dTdzeta;
    dTdz = Jinv3[0] * dTdxi + Jinv3[1] * dTdeta + Jinv3[2] * dTdzeta;

    for ( unsigned int i = 0; i < d_dofIndices.size(); i++ ) {
        lib_dTdx += ( dphidx[i][0] * T->getValueByGlobalID( d_dofIndices[i] ) );
        lib_dTdy += ( dphidy[i][0] * T->getValueByGlobalID( d_dofIndices[i] ) );
        lib_dTdz += ( dphidz[i][0] * T->getValueByGlobalID( d_dofIndices[i] ) );
    }

    fprintf( fp, " dT/dx = %.15f, dTdy = %.15f, dTdz = %.15f \n", dTdx, dTdy, dTdz );
    fprintf( fp,
             " lib_dT/dx = %.15f, lib_dTdy = %.15f, lib_dTdz = %.15f \n",
             lib_dTdx,
             lib_dTdy,
             lib_dTdz );
    fclose( fp );
    ut->passes( "Ran to completion" );
}

int testInverseJacobianBug( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::string exeName = "testInverseJacobianBug";

    myTest( &ut, exeName );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
