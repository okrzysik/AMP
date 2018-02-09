#include "AMP/materials/Material.h"
#include "AMP/operators/NeutronicsRhs.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/InputDatabase.h"
#include "AMP/utils/InputManager.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/shared_ptr.h"
#include "AMP/vectors/Variable.h"
#include <string>

#include "AMP/operators/diffusion/DiffusionLinearElement.h"
#include "AMP/operators/diffusion/DiffusionNonlinearElement.h"
#include "AMP/operators/diffusion/DiffusionTransportModel.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/Vector.h"

/* libMesh files */
#include "dense_matrix.h"
#include "dof_map.h"
#include "elem.h"
#include "equation_systems.h"
#include "fe.h"
#include "libmesh.h"
#include "libmesh_config.h"
#include "linear_implicit_system.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "petsc_matrix.h"
#include "petsc_vector.h"
#include "quadrature_gauss.h"
#include "sparse_matrix.h"

#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "AMP/operators/libmesh/MassLinearElement.h"
#include "AMP/operators/libmesh/MassLinearFEOperator.h"

#include "AMP/operators/DirichletMatrixCorrection.h"
#include "AMP/operators/DirichletVectorCorrection.h"
#include "AMP/operators/IsotropicElasticModel.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/MechanicsLinearElement.h"
#include "AMP/operators/MechanicsLinearFEOperator.h"

#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#include "AMP/time_integrators/sundials/IDATimeIntegrator.h"

#define ITFAILS ut.failure( __LINE__ );
#define UNIT_TEST( a ) \
    if ( !( a ) )      \
        ut.failure( __LINE__ );

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//


void IDATimeIntegratorTest( AMP::UnitTest *ut )
{

    AMP::AMP_MPI::initialize();
    AMP::AMPManager::startup();


    AMP::Materials::Material::shared_ptr my_material = AMP::Materials::create( "UO2_MSRZC_09" );
    AMP::Materials::PropertyType my_property =
        AMP::Materials::findProperty( "ThermalConductivity" );

    int num_components = 20;
    std::vector<double> result( num_components );
    std::vector<double> temp( num_components );
    std::vector<double> conc( num_components );
    std::vector<double> burn( num_components );

    for ( int i = 0; i < num_components; i++ ) {
        temp[i] = 600. + i * 10;
        conc[i] = .1;
        burn[i] = 0.;
    }
    my_material->evalv(
        my_property, &( result[0] ), &( temp[0] ), &conc[0], &burn[0], num_components );


    for ( int j = 0; j < num_components; j++ ) {
        cout << "result[" << j << "] = " << result[j] << endl;
    }

    if ( ut.numFails == 0 ) {
        ut.passes( "testIDATimeIntegrator successful" );
    }
}


//---------------------------------------------------------------------------//

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    IDATimeIntegratorTest( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}

//---------------------------------------------------------------------------//
//                        end of SundialsVectorTest.cc
//---------------------------------------------------------------------------//
