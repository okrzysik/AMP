#include "utils/Database.h"
#include "materials/Material.h"
#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"


#include "ampmesh/Mesh.h"
#include "utils/Writer.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/MultiVariable.h"
#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"

#include "operators/mechanics/MechanicsLinearFEOperator.h"
#include "operators/mechanics/MechanicsNonlinearFEOperator.h"
#include "operators/mechanics/PericElastoViscoPlasticModel.h"
#include "operators/mechanics/VonMisesElastoPlasticModel.h"
#include "operators/mechanics/VonMises_IsotropicKinematicHardening.h"

#include "operators/boundary/DirichletVectorCorrection.h"

#include "operators/BVPOperatorParameters.h"
#include "operators/ElementPhysicsModelFactory.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/OperatorBuilder.h"

#include "solvers/petsc/PetscKrylovSolver.h"
#include "solvers/petsc/PetscKrylovSolverParameters.h"
#include "solvers/petsc/PetscSNESSolver.h"
#include "solvers/petsc/PetscSNESSolverParameters.h"
#include "solvers/trilinos/TrilinosMLSolver.h"

#define ITFAILS ut.failure( __LINE__ );
#define UNIT_TEST( a ) \
    if ( !( a ) )      \
        ut.failure( __LINE__ );

void myTest( AMP::UnitTest *ut, std::string exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
    AMP::shared_ptr<AMP::Operator::MechanicsMaterialModel> mechanicsMaterialModel;
    AMP::shared_ptr<AMP::Operator::PericElastoViscoPlasticModel> pevpModel;

    AMP::shared_ptr<AMP::Database> matModelDatabase =
        input_db->getDatabase( "MechanicsMaterialModel" );
    elementPhysicsModel =
        AMP::Operator::ElementPhysicsModelFactory::createElementPhysicsModel( matModelDatabase );
    mechanicsMaterialModel =
        AMP::dynamic_pointer_cast<AMP::Operator::MechanicsMaterialModel>( elementPhysicsModel );
    pevpModel = AMP::dynamic_pointer_cast<AMP::Operator::PericElastoViscoPlasticModel>(
        mechanicsMaterialModel );

    pevpModel->preNonlinearInit( true, true );
    pevpModel->nonlinearInitGaussPointOperation( 500.0 );
    pevpModel->preNonlinearAssembly();

    int max_load = 850;
    // int max_num = 3 * max_load;

    double *stress;
    double *C;
    // double sigy = 0.2, E = 70.0, Ep = 0.135, H_a = 1.5, nu = 0.3333333, diff_norm;
    double eph11[1000], sig11[1000], sig11p[1000], eph11p[1000], slope[1000], slope_p[1000],
        time[1000];
    // double sig11_init[1000];
    std::vector<std::vector<double>> strain( 1 );
    strain[0].push_back( 0.0 );
    strain[0].push_back( 0.0 );
    strain[0].push_back( 0.0 );
    strain[0].push_back( 0.0 );
    strain[0].push_back( 0.0 );
    strain[0].push_back( 0.0 );

    double current_time = 0.0;
    double epsilon_dot  = 0.1;
    // double strain_const = 0.0;
    double time_multiplier = 1.0;
    double previous_time   = 0.0;
    double previous_strain = 0.0;

    for ( int i = 0; i < max_load; i++ ) {
        /*sig11_init[i] = ((double)(i + 1)) * 0.001;
        if(sig11_init[i] <= sigy) {
          strain[0][0] = sig11_init[i] / E;
          strain[0][1] = -nu * (sig11_init[i] / E);
          strain[0][2] = -nu * (sig11_init[i] / E);
        } else {
          strain[0][0] = (sig11_init[i] / E) + ((sig11_init[i] - sigy) / H_a);
          strain[0][1] = (-nu * (sig11_init[i] / E)) - ((1.0 / 2.0) * (sig11_init[i] - sigy) / H_a);
          strain[0][2] = (-nu * (sig11_init[i] / E)) - ((1.0 / 2.0) * (sig11_init[i] - sigy) / H_a);
          //strain[0][0] = (sig11_init[i] / E) + ((sig11_init[i] - sigy) / Ep);
          //strain[0][1] = (-nu * (sig11_init[i] / E)) - ((1.0 / 2.0) * (sig11_init[i] - sigy) /
        Ep);
          //strain[0][2] = (-nu * (sig11_init[i] / E)) - ((1.0 / 2.0) * (sig11_init[i] - sigy) /
        Ep);
        }*/

        if ( i <= 70 ) {
            epsilon_dot     = 0.1;
            time_multiplier = 0.001;
        } else if ( i <= 350 ) {
            epsilon_dot     = 0.01;
            time_multiplier = 0.001;
        } else if ( i <= 450 ) {
            epsilon_dot     = 0.1;
            time_multiplier = 0.001;
        } else if ( i <= 700 ) {
            epsilon_dot     = 0.01;
            time_multiplier = 0.001;
        } else if ( i <= 850 ) {
            epsilon_dot     = 0.1;
            time_multiplier = 0.001;
        }
        /*else if(i <= 600) {
          epsilon_dot = 0.01;
          time_multiplier = 0.001;
        } else if(i<= 650) {
          epsilon_dot = 0.1;
          time_multiplier = 0.001;
        }*/

        current_time = previous_time + time_multiplier;
        time[i]      = current_time;
        strain[0][0] = previous_strain + ( time_multiplier * epsilon_dot );
        mechanicsMaterialModel->updateTime( current_time );

        eph11[i] = strain[0][0];
        eph11p[i] =
            strain[0][0] - ( ( 1.0 / 3.0 ) * ( strain[0][0] + strain[0][1] + strain[0][2] ) );

        pevpModel->getInternalStress( strain, stress );

        sig11[i]  = stress[0];
        sig11p[i] = stress[0] - ( ( 1.0 / 3.0 ) * ( stress[0] + stress[1] + stress[2] ) );
        pevpModel->globalReset();

        if ( i == 0 ) {
            slope[0]   = 0.0;
            slope_p[0] = 0.0;
        } else {
            slope[i]   = ( sig11[i] - sig11[i - 1] ) / ( eph11[i] - eph11[i - 1] );
            slope_p[i] = ( sig11p[i] - sig11p[i - 1] ) / ( eph11p[i] - eph11p[i - 1] );
        }

        previous_time   = current_time;
        previous_strain = strain[0][0];
    }

    FILE *fout;
    fout = fopen( "pevp_stress_strain_results.txt", "w" );
    for ( int i = 0; i < max_load; i++ ) {
        // fprintf(fout,"%15.8lf%15.8lf%15.8lf%15.8lf%15.8lf%15.8lf%15.8lf\n",eph11[i],sig11[i],slope[i],eph11p[i],sig11p[i],slope_p[i],sig11_init[i]);
        fprintf( fout,
                 "%15.8f%15.8f%15.8f%15.8f%15.8f%15.8f%15.8f\n",
                 time[i],
                 eph11[i],
                 sig11[i],
                 slope[i],
                 eph11p[i],
                 sig11p[i],
                 slope_p[i] );
    }
    fclose( fout );

    // for(int i = 0; i < 6; i++) {
    //  printf("stress[%d]=%lf\n",i,stress[i]);
    //}

    // pevpModel->preNonlinearJacobian();
    // pevpModel->nonlinearJacobianGaussPointOperation(strain);
    // pevpModel->postNonlinearJacobianGaussPointOperation();

    pevpModel->preLinearAssembly();
    pevpModel->getConstitutiveMatrix( C );
    pevpModel->postLinearGaussPointOperation();

    // for(int i = 0; i < 6; i++) {
    //  for(int j = 0; j < 6; j++) {
    //    printf("C[%d][%d] = %lf\n",i,j,C[(i * 6) + j]);
    //  }
    //}

    /*diff_norm = 0.0;
    for(int i = 0; i < max_num; i++) {
      diff_norm += ((sig11_init[i] - sig11[i]) * (sig11_init[i] - sig11[i]));
    }
    diff_norm = sqrt(diff_norm) / sig11[max_num - 1];
    if(fabs(diff_norm) > (1.0e-8)) {
      ut->failure("Gauss point test for Linear Isotropic Hardening");
    } else {
      ut->passes("");
    }*/

    // std::cout << "diff_norm = " << fabs(diff_norm) << std::endl;

    ut->passes( exeName );
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.push_back( "testMaterialModels-PericElastoViscoPlasticity-1" );

    for ( auto &exeName : exeNames ) {
        try {
            myTest( &ut, exeName );
        } catch ( std::exception &err ) {
            AMP::pout << "ERROR: While testing " << argv[0] << err.what() << std::endl;
            ut.failure( "ERROR: While testing" );
        } catch ( ... ) {
            AMP::pout << "ERROR: While testing " << argv[0] << "An unknown exception was thrown."
                      << std::endl;
            ut.failure( "ERROR: While testing" );
        }
    } // end for i

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
