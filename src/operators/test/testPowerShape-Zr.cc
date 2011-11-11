#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <string>
#include "utils/AMPManager.h"
#include "ampmesh/MeshManager.h"
#include "boost/shared_ptr.hpp"
#include "utils/InputDatabase.h"
#include "utils/Utilities.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/Database.h"
#include "ampmesh/MeshAdapter.h"
#include "vectors/Variable.h"

#include "ampmesh/SiloIO.h"
#include "vectors/Vector.h"


#include "../PowerShape.h"


void test_with_shape(AMP::UnitTest *ut )
{

//--------------------------------------------------
//  Read Input File.
//--------------------------------------------------

    boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
    input_db->putInteger("NumberOfMeshes",1);
    boost::shared_ptr<AMP::Database> mesh_db = input_db->putDatabase("Mesh_1");
    mesh_db->putString("Filename","cylinder270.e");
    mesh_db->putString("MeshName","fuel");
    mesh_db->putDouble("x_offset",0.);
    mesh_db->putDouble("y_offset",0.);
    mesh_db->putDouble("z_offset",0.);
//--------------------------------------------------
//   Create the Mesh.
//--------------------------------------------------
    AMP::Mesh::MeshManagerParameters::shared_ptr mgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
    AMP::Mesh::MeshManager::shared_ptr manager ( new AMP::Mesh::MeshManager ( mgrParams ) );
    AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = manager->getMesh ( "fuel" );

//--------------------------------------------------
//  Construct PowerShape for a radial only term.
//--------------------------------------------------
    boost::shared_ptr<AMP::Database> shape_db = input_db->putDatabase("shape_db");
    shape_db->putString("coordinateSystem","cylindrical");
    shape_db->putString("type","zernikeRadial");
    shape_db->putInteger("print_info_level",10);
    
    #ifdef USE_SILO
      AMP::Operator::PowerShape::SP_HexGaussPointVariable shapeVar( new AMP::Operator::PowerShape::HexGaussPointVariable("PowerShape") );
      AMP::LinearAlgebra::Vector::shared_ptr shapeVec = meshAdapter->createVector( shapeVar );
      meshAdapter->registerVectorAsData ( shapeVec );
    #endif
      
    for ( int nMoments = 0; nMoments < 3; nMoments++ ) {
      shape_db->putInteger("numMoments", nMoments); 
      if( nMoments > 0) {
        std::vector<double> moments(nMoments,0.);
        moments[nMoments-1] = -1.;
        shape_db->putDoubleArray("Moments", moments);
      } 
      boost::shared_ptr<AMP::Operator::PowerShapeParameters> shape_params(new AMP::Operator::PowerShapeParameters( shape_db ));
      shape_params->d_MeshAdapter = meshAdapter;
      boost::shared_ptr<AMP::Operator::PowerShape> shape(new AMP::Operator::PowerShape( shape_params ));
      AMP::Operator::PowerShape::SP_HexGaussPointVariable SpecificPowerShapeVar = 
                                                         shape->createOutputVariable("SpecificPowerShape");
    
      // Create a vector associated with the Variable.
      AMP::LinearAlgebra::Vector::shared_ptr  SpecificPowerShapeVec     = meshAdapter->createVector( SpecificPowerShapeVar );
      AMP::LinearAlgebra::Vector::shared_ptr  SpecificPowerMagnitudeVec = SpecificPowerShapeVec->cloneVector( SpecificPowerShapeVar );
      SpecificPowerMagnitudeVec->setToScalar(1.0);
    
      // Set the initial value for all nodes of SpecificPowerVec to zero.
      SpecificPowerShapeVec->setToScalar(0.0);
      AMP::LinearAlgebra::Vector::shared_ptr   nullVec;
      try   { shape->apply(nullVec, SpecificPowerMagnitudeVec, SpecificPowerShapeVec, 1., 0.); }
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
      
      #ifdef USE_SILO
        shapeVec->copyVector(SpecificPowerShapeVec);
        // SpecificPowerShapeVec->copyVector(shapeVec);
        manager->writeFile<AMP::Mesh::SiloIO> ( "PowerShape-Zr", nMoments );
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
      boost::shared_ptr<AMP::Operator::PowerShapeParameters> shape_params(new AMP::Operator::PowerShapeParameters( shape_db ));
      shape_params->d_MeshAdapter = meshAdapter;
      boost::shared_ptr<AMP::Operator::PowerShape> shape(new AMP::Operator::PowerShape( shape_params ));
      AMP::Operator::PowerShape::SP_HexGaussPointVariable SpecificPowerShapeVar = 
                                                         shape->createOutputVariable("SpecificPowerShape");
    
      // Create a vector associated with the Variable.
      AMP::LinearAlgebra::Vector::shared_ptr  SpecificPowerShapeVec     = meshAdapter->createVector( SpecificPowerShapeVar );
      AMP::LinearAlgebra::Vector::shared_ptr  SpecificPowerMagnitudeVec = SpecificPowerShapeVec->cloneVector( SpecificPowerShapeVar );
      SpecificPowerMagnitudeVec->setToScalar(1.0);
    
      // Set the initial value for all nodes of SpecificPowerVec to zero.
      SpecificPowerShapeVec->setToScalar(0.0);
      AMP::LinearAlgebra::Vector::shared_ptr   nullVec;
      try   { shape->apply(nullVec, SpecificPowerMagnitudeVec, SpecificPowerShapeVec, 1., 0.); }
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
      
      #ifdef USE_SILO
        shapeVec->copyVector(SpecificPowerShapeVec);
        // SpecificPowerShapeVec->copyVector(shapeVec);
        manager->writeFile<AMP::Mesh::SiloIO> ( "PowerShape-Zr", i+3 );
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



