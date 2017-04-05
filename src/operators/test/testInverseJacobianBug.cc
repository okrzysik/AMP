
#include "utils/shared_ptr.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>

#include "utils/AMPManager.h"
#include "utils/InputManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/libmesh/libMesh.h"

#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"

#include "libmesh/boundary_info.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/elem.h"
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/node.h"

#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_type.h"
#include "libmesh/quadrature.h"

#include "libmesh/auto_ptr.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/string_to_enum.h"

void myTest( AMP::UnitTest *ut, std::string exeName )
{
    std::string log_file = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    AMP::shared_ptr<::Mesh> mesh( new ::Mesh( 3 ) );
    MeshTools::Generation::build_cube(
        ( *( mesh.get() ) ), 1, 1, 1, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, HEX8, false );

    ::Elem *elemPtr = mesh->elem( 0 );

    ( elemPtr->point( 4 ) )( 0 ) -= 0.4;

    ( elemPtr->point( 5 ) )( 0 ) -= 0.4;

    ( elemPtr->point( 6 ) )( 0 ) -= 0.4;

    ( elemPtr->point( 7 ) )( 0 ) -= 0.4;

    AMP::Mesh::Mesh::shared_ptr ampMesh =
        AMP::Mesh::Mesh::shared_ptr( new AMP::Mesh::libMesh( mesh, "TestMesh" ) );

    AMP::LinearAlgebra::Variable::shared_ptr myVar( new AMP::LinearAlgebra::Variable( "myVar" ) );
    AMP::Discretization::DOFManager::shared_ptr dof_map =
        AMP::Discretization::simpleDOFManager::create( ampMesh, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    AMP::LinearAlgebra::Vector::shared_ptr T =
        AMP::LinearAlgebra::createVector( dof_map, myVar, true );

    FILE *fp;
    fp = fopen( "InverseJacobian.txt", "w" );

    AMP::Mesh::MeshIterator el     = ampMesh->getIterator( AMP::Mesh::GeomType::Volume, 0 );
    AMP::Mesh::MeshIterator end_el = el.end();

    while ( el != end_el ) {
        std::vector<AMP::Mesh::MeshElement> nodes = el->getElements( AMP::Mesh::GeomType::Vertex );
        for ( size_t i = 0; i < nodes.size(); i++ ) {
            std::vector<double> pt = nodes[i].coord();
            fprintf(
                fp, "nd = %d, x = %.15f, y = %.15f, z = %.15f \n", (int) i, pt[0], pt[1], pt[2] );
        }
    }

    libMeshEnums::Order feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>( "FIRST" );
    libMeshEnums::FEFamily feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>( "LAGRANGE" );

    AMP::shared_ptr<::FEType> feType( new ::FEType( feTypeOrder, feFamily ) );
    AMP::shared_ptr<::FEBase> fe( (::FEBase::build( 3, ( *feType ) ) ).release() );

    //  const std::vector<std::vector<Real> > &phi = fe->get_phi();
    const std::vector<std::vector<Real>> &dphidxi   = fe->get_dphidxi();
    const std::vector<std::vector<Real>> &dphideta  = fe->get_dphideta();
    const std::vector<std::vector<Real>> &dphidzeta = fe->get_dphidzeta();
    const std::vector<std::vector<Real>> &dphidx    = fe->get_dphidx();
    const std::vector<std::vector<Real>> &dphidy    = fe->get_dphidy();
    const std::vector<std::vector<Real>> &dphidz    = fe->get_dphidz();
    //  const std::vector<Real> &djxw = fe->get_JxW();

    const std::vector<RealGradient> &dxyzdxi   = fe->get_dxyzdxi();
    const std::vector<RealGradient> &dxyzdeta  = fe->get_dxyzdeta();
    const std::vector<RealGradient> &dxyzdzeta = fe->get_dxyzdzeta();

    libMeshEnums::QuadratureType qruleType =
        Utility::string_to_enum<libMeshEnums::QuadratureType>( "QGAUSS" );
    libMeshEnums::Order qruleOrder = feType->default_quadrature_order();
    AMP::shared_ptr<::QBase> qrule( (::QBase::build( qruleType, 3, qruleOrder ) ).release() );

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

    ::Elem *currElemPtr;
    currElemPtr = new ::Hex8;
    for ( size_t j = 0; j < d_currNodes.size(); j++ ) {
        std::vector<double> pt     = d_currNodes[j].coord();
        currElemPtr->set_node( j ) = new ::Node( pt[0], pt[1], pt[2], j );
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

    std::vector<Real> Jinv1( 3 );
    std::vector<Real> Jinv2( 3 );
    std::vector<Real> Jinv3( 3 );
    Jinv1[0] = 2.;
    Jinv1[1] = 0.;
    Jinv1[2] = 0.;

    Jinv2[0] = 0;
    Jinv2[1] = 2;
    Jinv2[2] = 0;

    Jinv3[0] = 0.8;
    Jinv3[1] = 0.;
    Jinv3[2] = 2.;

    Real dTdxi = 0, dTdeta = 0, dTdzeta = 0, dTdx = 0;
    Real dTdy = 0, dTdz = 0, lib_dTdx = 0, lib_dTdy = 0, lib_dTdz = 0;

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

int main( int argc, char *argv[] )
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
