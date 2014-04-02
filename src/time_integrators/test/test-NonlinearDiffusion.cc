#include <string>
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/AMP_MPI.h"
#include "materials/Material.h"
#include "boost/shared_ptr.hpp"
#include "utils/InputDatabase.h"
#include "utils/Utilities.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/Database.h"
#include "operators/NeutronicsRhs.h"
#include "vectors/Variable.h"

#include "utils/Writer.h"
#include "vectors/Vector.h"
#include "operators/diffusion/DiffusionLinearElement.h"
#include "operators/diffusion/DiffusionNonlinearElement.h"
#include "operators/diffusion/DiffusionTransportModel.h"

/* libMesh files */
#include "libmesh_config.h"
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "equation_systems.h"
#include "fe.h"
#include "quadrature_gauss.h"
#include "dof_map.h"
#include "sparse_matrix.h"
#include "petsc_matrix.h"
#include "petsc_vector.h"
#include "dense_matrix.h"
#include "linear_implicit_system.h"
#include "elem.h"

#include "operators/libmesh/MassLinearElement.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "operators/libmesh/MassLinearFEOperator.h"

#include "operators/IsotropicElasticModel.h"
#include "operators/MechanicsLinearElement.h"
#include "operators/MechanicsLinearFEOperator.h"
#include "operators/DirichletMatrixCorrection.h"
#include "operators/DirichletVectorCorrection.h"
#include "operators/LinearBVPOperator.h"

#include "solvers/trilinos/TrilinosMLSolver.h"
#include "time_integrators/sundials/IDATimeIntegrator.h"

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//



void IDATimeIntegratorTest(AMP::UnitTest *ut)
{
    
    AMP::AMP_MPI::initialize();
    AMP::AMPManager::startup();
    
    
    AMP::Materials::Material::shared_ptr my_material = AMP::Materials::create("UO2_MSRZC_09");
    AMP::Materials::PropertyType         my_property = AMP::Materials::findProperty("ThermalConductivity");
    
    int num_components = 20;
    std::vector<double> result(num_components);
    std::vector<double> temp(num_components);
    std::vector<double> conc(num_components);
    std::vector<double> burn(num_components);
    
    for (int i=0; i<num_components; i++){
        temp[i] = 600.+i*10;
        conc[i] = .1;
        burn[i]=0.;
    }
    my_material->evalv(my_property, &(result[0]), &(temp[0]), &conc[0], &burn[0], num_components);
    
    
    for (int j=0; j<num_components; j++){
        cout << "result[" << j << "] = " << result[j] << endl;
    }    
        
    if (ut.numFails == 0)
    {
        ut.passes("testIDATimeIntegrator successful");
    }
}


//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;
    
    try {        
        IDATimeIntegratorTest(&ut);        
    } catch (std::exception &err) {
        std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
        ut.failure("ERROR: While testing");
    } catch( ... ) {
        std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
        ut.failure("ERROR: While testing");
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
    
}   

//---------------------------------------------------------------------------//
//                        end of SundialsVectorTest.cc
//---------------------------------------------------------------------------//










