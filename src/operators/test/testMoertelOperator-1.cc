#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "operators/MoertelOperatorBuilder.h"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"


#include "ampmesh/MeshManager.h"
#include "ampmesh/MeshVariable.h"
#include "ampmesh/SiloIO.h"
#include "ampmesh/MeshUtils.h"

#include "../../ampmesh/test/test_Base.h"

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;


  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile("input_testMoertelOperator-1", input_db);
  input_db->printClassData(AMP::plog);

  AMP::MeshManagerParameters::shared_ptr  meshmgrParams ( new AMP::MeshManagerParameters ( input_db ) );
  AMP::MeshManager::shared_ptr  manager ( new AMP::MeshManager ( meshmgrParams ) );

  boost::shared_ptr<AMP::MoertelOperatorBuilderParameters> params ( new AMP::MoertelOperatorBuilderParameters () );
  params->d_MasterVolume = manager->getMesh ( "master" );
  params->d_SlaveVolume = manager->getMesh ( "slave" );
  params->d_DB = input_db;

  AMP::LinearAlgebra::Variable::shared_ptr tmp_slave ( new AMP::NodalScalarVariable ( "temp" , manager->getMesh( "slave" ) ) );
  AMP::LinearAlgebra::Variable::shared_ptr tmp_master ( new AMP::NodalScalarVariable ( "temp" , manager->getMesh( "master" ) ) );
  AMP::LinearAlgebra::MultiVariable  *newVar = new AMP::LinearAlgebra::MultiVariable ( "temp" );
  newVar->add ( tmp_slave );
  newVar->add ( tmp_master );
  AMP::LinearAlgebra::Variable::shared_ptr temperature_var ( newVar );

  AMP::LinearAlgebra::Vector::shared_ptr temperature = manager->createVector ( temperature_var );
  AMP::LinearAlgebra::Vector::shared_ptr master_temp = temperature->subsetVectorForVariable ( tmp_master );
  AMP::LinearAlgebra::Vector::shared_ptr slave_temp = temperature->subsetVectorForVariable ( tmp_slave );

  params->d_MasterDOFMap = manager->getMesh ( "master" )->getDOFMap ( tmp_master );
  params->d_SlaveDOFMap = manager->getMesh ( "slave" )->getDOFMap ( tmp_slave );

  params->d_MasterSurface = 4;
  params->d_SlaveSurface = 8;
  params->d_Comm = globalComm;

  manager->registerVectorAsData ( temperature , "contact" );

  if (AMP::AMP_MPI::getNodes() > 1) { 
    nut.passes("It must have at least one pass to not be an unexpected failure."); 
    nut.expected_failure("Moertel will fail in multi-core with an STL assertion."); 
  } else{
    AMP::MoertelOperatorBuilder  opBuilder ( params );

    temperature->setToScalar ( 0.0 );
    AMP::LinearAlgebra::Vector::shared_ptr lambda1 = opBuilder.createLambdaVector();
    AMP::LinearAlgebra::Vector::shared_ptr lambda2 = opBuilder.createLambdaVector();
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    boost::shared_ptr<AMP::Operator::Operator>  MOper = opBuilder.createMOperator( master_temp->getVariable() , lambda1->getVariable() );
    boost::shared_ptr<AMP::Operator::Operator>  DOper = opBuilder.createDOperator( slave_temp->getVariable() , lambda2->getVariable() );
    boost::shared_ptr<AMP::Operator::Operator>  MTOper = opBuilder.createMTOperator( lambda1->getVariable() , master_temp->getVariable() );
    boost::shared_ptr<AMP::Operator::Operator>  DTOper = opBuilder.createDTOperator( lambda2->getVariable() , slave_temp->getVariable() );

    master_temp->setToScalar ( 300.0 );
    MOper->apply ( nullVec , master_temp , lambda1 , 1.0 , 0 );

    slave_temp->setToScalar ( 300.0 );
    DOper->apply ( nullVec , slave_temp , lambda2 , 1.0 , 0 );

    manager->writeFile<AMP::SiloIO> ( "contact" , 1 );
    lambda1->add ( lambda1 , lambda2 );
    lambda1->divide ( lambda1 , lambda2 );

  if ( lambda1->maxNorm () < 0.05 )
    nut.passes ( "Moertel integrates correctly compute lambda = 0 to within 1 part in 20 for equal temps." );
  else
    nut.failure ( "Moertel integrates incorrectly compute lambda = 0 to within 1 part in 20 for equal temps." );
 }

  nut.report();
}



