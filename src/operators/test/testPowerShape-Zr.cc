#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <string>
#include "utils/AMPManager.h"
#include "utils/shared_ptr.h"
#include "utils/InputDatabase.h"
#include "utils/Utilities.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/Database.h"
#include "vectors/Variable.h"

#include "utils/Writer.h"
#include "vectors/Vector.h"

#include "ampmesh/Mesh.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"

#include "operators/libmesh/PowerShape.h"


void test_with_shape(AMP::UnitTest *ut )
{

//--------------------------------------------------
//  Read Input File.
//--------------------------------------------------

    AMP::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
    AMP::shared_ptr<AMP::Database> mesh_db = input_db->putDatabase("Mesh");
    mesh_db->putString("FileName","cylinder270.e");
    mesh_db->putString("MeshType","libMesh");
    mesh_db->putString("MeshName","fuel");
    mesh_db->putInteger("dim",3);
    mesh_db->putDouble("x_offset",0.);
    mesh_db->putDouble("y_offset",0.);
    mesh_db->putDouble("z_offset",0.);
//--------------------------------------------------
//   Create the Mesh.
//--------------------------------------------------
    AMP::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(new AMP::Mesh::MeshParameters(mesh_db));
    mgrParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
    AMP::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh(mgrParams);

//--------------------------------------------------
//  Construct PowerShape for a radial only term.
//--------------------------------------------------
    AMP::shared_ptr<AMP::Database> shape_db = input_db->putDatabase("shape_db");
    shape_db->putString("coordinateSystem","cylindrical");
    shape_db->putString("type","zernikeRadial");
    shape_db->putInteger("print_info_level",1);
    
    // Create a DOF manager for a gauss point vector 
    int DOFsPerNode = 8;
    int ghostWidth = 0;
    bool split = true;
    AMP::Discretization::DOFManager::shared_ptr dof_map = AMP::Discretization::simpleDOFManager::create(meshAdapter, AMP::Mesh::Volume, ghostWidth, DOFsPerNode, split);

    #ifdef USE_EXT_SILO
      AMP::LinearAlgebra::Variable::shared_ptr shapeVar(new AMP::LinearAlgebra::Variable("PowerShape"));
      AMP::LinearAlgebra::Vector::shared_ptr shapeVec = AMP::LinearAlgebra::createVector( dof_map, shapeVar, split );
      AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter("Silo");
      siloWriter->registerMesh( meshAdapter );
      siloWriter->registerVector( shapeVec, meshAdapter, AMP::Mesh::Volume, "PowerShape" );
    #endif
      
    for ( int nMoments = 0; nMoments < 3; nMoments++ ) {
      shape_db->putInteger("numMoments", nMoments); 
      if( nMoments > 0) {
        std::vector<double> moments(nMoments,0.);
        moments[nMoments-1] = -1.;
        shape_db->putDoubleArray("Moments", moments);
      } 
      AMP::shared_ptr<AMP::Operator::PowerShapeParameters> shape_params(new AMP::Operator::PowerShapeParameters( shape_db ));
      shape_params->d_Mesh = meshAdapter;
      AMP::shared_ptr<AMP::Operator::PowerShape> shape(new AMP::Operator::PowerShape( shape_params ));
    
      // Create a shared pointer to a Variable - Power - Output because it will be used in the "residual" location of apply. 
      AMP::LinearAlgebra::Variable::shared_ptr SpecificPowerShapeVar(new AMP::LinearAlgebra::Variable("SpecificPowerInWattsPerKg"));
  
      // Create a vector associated with the Variable.
      AMP::LinearAlgebra::Vector::shared_ptr SpecificPowerShapeVec = AMP::LinearAlgebra::createVector( dof_map, SpecificPowerShapeVar, split );
      AMP::LinearAlgebra::Vector::shared_ptr  SpecificPowerMagnitudeVec = SpecificPowerShapeVec->cloneVector( );
      SpecificPowerMagnitudeVec->setToScalar(1.0);
    
      // Set the initial value for all nodes of SpecificPowerVec to zero.
      SpecificPowerShapeVec->setToScalar(0.0);
      AMP::LinearAlgebra::Vector::shared_ptr   nullVec;
      try   { shape->apply(SpecificPowerMagnitudeVec, SpecificPowerShapeVec); }
      catch ( std::exception const & a ) {  
        std::cout << a.what() << std::endl;  
        ut->failure("error");
      }
      if( nMoments == 0 ) { 
        if( !AMP::Utilities::approx_equal(SpecificPowerShapeVec->max(), 1.0, 1e-9) ) {
          ut->failure("flat solution is not really flat (max).");
          printf("This %.9e is not 1.0. \n", SpecificPowerShapeVec->max() ); 
        }
        if( !AMP::Utilities::approx_equal(SpecificPowerShapeVec->min(), 1.0, 1e-9) ) {
          ut->failure("flat solution is not really flat (min).");
          printf("This %.9e is not 1.0. \n", SpecificPowerShapeVec->min() ); 
        }
      }
      
      #ifdef USE_EXT_SILO
        shapeVec->copyVector(SpecificPowerShapeVec);
        // SpecificPowerShapeVec->copyVector(shapeVec);
        siloWriter->writeFile( "PowerShape-Zr" , nMoments );
      #endif
      
    }//end for elements

    ut->passes("PowerShape produces a non-negative power shape.");

//--------------------------------------------------
//  Construct PowerShape for a full Zernike basis. -
//--------------------------------------------------
    shape_db->putString("type","zernike");
    int i=0;
    for ( int nMoments = 0; nMoments < 9; nMoments++ ) {
//    for ( int nMoments = 0; nMoments < 5; nMoments++ ) {
    for ( int n = -nMoments; n<= nMoments; i++ ) {
      shape_db->putInteger("numMoments", nMoments); 
      int nMNmoments = (nMoments+2)*(nMoments+1)/2-1; 
      if( nMoments > 0) {
        std::vector<double> moments(nMNmoments,0.);
        moments[i-1] = -1.;
        shape_db->putDoubleArray("Moments", moments);
      } 
      AMP::shared_ptr<AMP::Operator::PowerShapeParameters> shape_params(new AMP::Operator::PowerShapeParameters( shape_db ));
      shape_params->d_Mesh = meshAdapter;
      AMP::shared_ptr<AMP::Operator::PowerShape> shape(new AMP::Operator::PowerShape( shape_params ));
      AMP::LinearAlgebra::Variable::shared_ptr SpecificPowerShapeVar(new AMP::LinearAlgebra::Variable("SpecificPowerShape") );
    
      // Create a vector associated with the Variable.
      AMP::LinearAlgebra::Vector::shared_ptr SpecificPowerShapeVec = AMP::LinearAlgebra::createVector( dof_map, SpecificPowerShapeVar, split );
      AMP::LinearAlgebra::Vector::shared_ptr  SpecificPowerMagnitudeVec = SpecificPowerShapeVec->cloneVector( );
      SpecificPowerMagnitudeVec->setToScalar(1.0);
    
      // Set the initial value for all nodes of SpecificPowerVec to zero.
      SpecificPowerShapeVec->setToScalar(0.0);
      AMP::LinearAlgebra::Vector::shared_ptr   nullVec;
      try   { shape->apply( SpecificPowerMagnitudeVec, SpecificPowerShapeVec); }
      catch ( std::exception const & a ) {  
        std::cout << a.what() << std::endl;  
        ut->failure("PowerShape error");
      }
      if( nMoments == 0 ) { 
        if( !AMP::Utilities::approx_equal(SpecificPowerShapeVec->max(), 1.0, 1e-9) ) {
          ut->failure("flat solution is not flat (max).");
          printf("This %.9e is not 1.0. \n", SpecificPowerShapeVec->max() ); 
        }
        if( !AMP::Utilities::approx_equal(SpecificPowerShapeVec->min(), 1.0, 1e-9) ) {
          ut->failure("flat solution is not flat (min).");
          printf("This %.9e is not 1.0. \n", SpecificPowerShapeVec->min() ); 
        }
      }
      
      #ifdef USE_EXT_SILO
        shapeVec->copyVector(SpecificPowerShapeVec);
        // SpecificPowerShapeVec->copyVector(shapeVec);
        siloWriter->writeFile( "PowerShape-Zr" , i+3 );
      #endif

      n+=2; 
    }
    }

    ut->passes("PowerShape produces a non-negative power shape.");

}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    try {
        test_with_shape(&ut);
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



