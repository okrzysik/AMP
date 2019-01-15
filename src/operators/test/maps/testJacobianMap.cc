#include "AMP/utils/AMPManager.h"
#include "AMP/utils/UnitTest.h"
#include <cstdlib>
#include <iostream>

#include "AMP/utils/ReadTestMesh.h"
#include "AMP/utils/shared_ptr.h"

// Libmesh files
DISABLE_WARNINGS
#include "libmesh/auto_ptr.h"
#include "libmesh/boundary_info.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_type.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/quadrature.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/string_to_enum.h"
ENABLE_WARNINGS


// Using mesh and function calls from testLibmeshGeomType::FaceStuff.cc

static void calculateGrad( AMP::UnitTest *ut )
{
    const unsigned int mesh_dim = 3;
    AMP::shared_ptr<Mesh> mesh( new Mesh( mesh_dim ) );

    AMP::readTestMesh( "distortedElementMesh", mesh );

    MeshCommunication().broadcast( *( mesh.get() ) );

    mesh->prepare_for_use( false );

    EquationSystems equation_systems( *( mesh.get() ) );

    auto &system = equation_systems.add_system<LinearImplicitSystem>( "Poisson" );

    system.add_variable( "V", FIRST );
    equation_systems.init();

    const unsigned int V_var = system.variable_number( "V" );

    FEType fe_type = system.variable_type( V_var );

    QGauss qrule( 3, fe_type.default_quadrature_order() );

    AutoPtr<FEBase> fe_3d( FEBase::build( 3, fe_type ) );
    fe_3d->attach_quadrature_rule( &qrule );
    // const std::vector<Point>& q_point3d = fe_3d->get_xyz();

    const std::vector<std::vector<Real>> &dphi3d = fe_3d->get_dphidx();

    MeshBase::const_element_iterator el           = mesh->local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh->local_elements_end();
    std::cout << "Entering Element Iyerator" << std::endl;
    for ( ; el != end_el; ++el ) {
        const Elem *elem = *el;

        fe_3d->reinit( elem );
        // std::vector<Point> coordinates = fe_3d->get_xyz();
        std::vector<double> computedAtGauss( qrule.n_points(), 0.0 );
        std::cout << "Entering Gauss Point loop : " << qrule.n_points() << std::endl;
        for ( unsigned int qp = 0; qp < qrule.n_points(); qp++ ) {
            std::cout << "dphidx.size = " << dphi3d.size() << std::endl;
            for ( size_t l = 0; l < dphi3d.size(); l++ ) {
                std::cout << "dphidx[" << l << "][" << qp << "]  = " << dphi3d[l][qp] << std::endl;
            }
        }
    }
    ut->passes( "Ran to completion" );
}

int testJacobianMap( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::cout << "Entering main" << std::endl;
    calculateGrad( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
