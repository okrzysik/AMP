#include <string>
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"

#include "boost/shared_ptr.hpp"

#include "utils/InputDatabase.h"
#include "utils/Utilities.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/Database.h"

#include "ampmesh/SiloIO.h"
#include "ampmesh/Mesh.h"
#include "ampmesh/MultiMesh.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"

#include "vectors/Variable.h"
#include "vectors/VectorBuilder.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/VectorSelector.h"

#include "operators/boundary/ColumnBoundaryOperator.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/boundary/RobinVectorCorrection.h"
#include "operators/boundary/NeumannVectorCorrectionParameters.h"
#include "operators/boundary/DirichletMatrixCorrection.h"
#include "operators/boundary/RobinMatrixCorrection.h"
#include "operators/boundary/NeumannVectorCorrection.h"

#include "operators/VectorCopyOperator.h"

#include "operators/VolumeIntegralOperator.h"
#include "operators/CoupledOperator.h"
#include "operators/CoupledOperatorParameters.h"
#include "operators/NodeToGaussPointOperator.h"
#include "operators/map/AsyncMapColumnOperator.h"
#include "operators/map/SubchannelToCladMap.h"
#include "operators/map/CladToSubchannelMap.h"
#include "operators/map/NodeToNodeMap.h"
#include "operators/map/ScalarZAxisMap.h"
#include "operators/map/ScalarN2GZAxisMap.h"
#include "operators/map/GaussPointToGaussPointMap.h"
#include "operators/map/Map3Dto1D.h" 
#include "operators/map/Map1Dto3D.h" 
#include "operators/map/MapSurface.h" 

#include "operators/ColumnOperator.h"
#include "operators/CoupledOperator.h"
#include "operators/ElementPhysicsModelFactory.h"
#include "operators/ElementOperationFactory.h"
#include "operators/OperatorBuilder.h"
#include "operators/subchannel/SubchannelTwoEqNonlinearOperator.h"
#include "operators/subchannel/SubchannelTwoEqLinearOperator.h"
#include "operators/subchannel/SubchannelPhysicsModel.h"
#include "operators/subchannel/CoupledChannelToCladMapOperator.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"

#include "solvers/PetscKrylovSolver.h"
#include "solvers/PetscSNESSolver.h"
#include "solvers/TrilinosMLSolver.h"
#include "solvers/ColumnSolver.h"

#include "ampmesh/StructuredMeshHelper.h"

double getPower(const std::vector<double> &x) {
    return 5e6 *( 1.0 + x[0]*50 + x[1]*100 + x[2] *500 );

}

void SubchannelSolve(AMP::UnitTest *ut, std::string exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file = "output_" + exeName;
    std::string silo_name = exeName;
    AMP::PIO::logAllNodes(log_file);
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
    globalComm.barrier();

    // Read the input file
    boost::shared_ptr<AMP::InputDatabase>  global_input_db ( new AMP::InputDatabase ( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile ( input_file , global_input_db );
    global_input_db->printClassData(AMP::plog);

    // Get the Mesh database and create the mesh parameters
    boost::shared_ptr<AMP::Database> database = global_input_db->getDatabase( "Mesh" );
    boost::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(database));
    meshParams->setComm(globalComm);

    AMP::LinearAlgebra::Vector::shared_ptr  globalSolMultiVector; 
    AMP::LinearAlgebra::Vector::shared_ptr  globalRhsMultiVector; 
    AMP::LinearAlgebra::Vector::shared_ptr  globalResMultiVector; 

    AMP::LinearAlgebra::Vector::shared_ptr globalThermalSolutionVec;
    AMP::LinearAlgebra::Vector::shared_ptr globalThermalRhsVec;
    AMP::LinearAlgebra::Vector::shared_ptr globalThermalResidualVec;

    AMP::LinearAlgebra::Vector::shared_ptr  specificPowerGpVec;

    AMP::LinearAlgebra::Vector::shared_ptr flowSolVec;
    AMP::LinearAlgebra::Vector::shared_ptr flowRhsVec;
    AMP::LinearAlgebra::Vector::shared_ptr flowResVec;

    boost::shared_ptr<AMP::Solver::PetscSNESSolver>  nonlinearCoupledSolver;
    boost::shared_ptr<AMP::Solver::PetscKrylovSolver>  linearColumnSolver;
    boost::shared_ptr<AMP::Solver::ColumnSolver>  columnPreconditioner;

    AMP::Mesh::Mesh::shared_ptr manager = AMP::Mesh::Mesh::buildMesh(meshParams);
    AMP::Mesh::Mesh::shared_ptr pinMesh = manager->Subset("MultiPin");
    AMP::Mesh::Mesh::shared_ptr cladMesh;
    if ( pinMesh.get()!=NULL ) {
      pinMesh->setName("MultiPin");
      cladMesh = pinMesh->Subset("clad");
    }
    AMP::Mesh::Mesh::shared_ptr subchannelMesh = manager->Subset("subchannel");

    AMP::Mesh::Mesh::shared_ptr xyFaceMesh;
    xyFaceMesh = subchannelMesh->Subset( AMP::Mesh::StructuredMeshHelper::getXYFaceIterator( subchannelMesh , 0 ) );

    AMP::LinearAlgebra::Variable::shared_ptr thermalVariable (new AMP::LinearAlgebra::Variable("Temperature"));    // temperature on pellets and cladding
    AMP::LinearAlgebra::Variable::shared_ptr flowVariable (new AMP::LinearAlgebra::Variable("Flow"));              // enthalpy and pressure on z-faces
    AMP::LinearAlgebra::Variable::shared_ptr powerVariable (new AMP::LinearAlgebra::Variable("SpecificPowerInWattsPerGram")); // specific power on gauss points

    // Dof manager for subchannel with flow variable
    int DofsPerFace =  2;
    AMP::Discretization::DOFManager::shared_ptr faceDOFManager = AMP::Discretization::simpleDOFManager::create( subchannelMesh, 
        AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,1), AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,0), DofsPerFace );

    // dof manager for the scalar quanties on the z faces - for mapped temperature clad temp onto subchannel discretization
    AMP::Discretization::DOFManager::shared_ptr scalarFaceDOFManager = AMP::Discretization::simpleDOFManager::create( subchannelMesh, 
        AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,1), AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,0), 1);

    AMP::LinearAlgebra::Vector::shared_ptr subchannelFuelTemp = AMP::LinearAlgebra::createVector( scalarFaceDOFManager , thermalVariable );
    AMP::LinearAlgebra::Vector::shared_ptr subchannelFlowTemp = AMP::LinearAlgebra::createVector( scalarFaceDOFManager , thermalVariable );
   
    AMP::Discretization::DOFManager::shared_ptr nodalScalarDOF = AMP::Discretization::simpleDOFManager::create( pinMesh ,AMP::Mesh::Vertex,1,1,true);

    // flow temperature on clad outer surfaces and pellet temperature on clad innner surfaces, and clad inner surface temp on pellet outer surfaces
    AMP::LinearAlgebra::Vector::shared_ptr thermalMapVec = AMP::LinearAlgebra::createVector(nodalScalarDOF, thermalVariable , true);

    // create solution, rhs, and residual vectors
    flowSolVec = AMP::LinearAlgebra::createVector( faceDOFManager , flowVariable , true );
    flowRhsVec = AMP::LinearAlgebra::createVector( faceDOFManager , flowVariable , true );
    flowResVec = AMP::LinearAlgebra::createVector( faceDOFManager , flowVariable , true );

    globalThermalSolutionVec =  AMP::LinearAlgebra::createVector( nodalScalarDOF , thermalVariable )  ;
    globalThermalRhsVec      =  AMP::LinearAlgebra::createVector( nodalScalarDOF , thermalVariable )  ;
    globalThermalResidualVec =  AMP::LinearAlgebra::createVector( nodalScalarDOF , thermalVariable )  ;

    globalSolMultiVector = AMP::LinearAlgebra::MultiVector::create( "multivector" , globalComm ) ;
    globalSolMultiVector->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( globalThermalSolutionVec );
    globalSolMultiVector->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( flowSolVec );

    globalRhsMultiVector = AMP::LinearAlgebra::MultiVector::create( "multivector" , globalComm ) ;
    globalRhsMultiVector->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( globalThermalRhsVec );
    globalRhsMultiVector->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( flowRhsVec );

    globalResMultiVector = AMP::LinearAlgebra::MultiVector::create( "multivector" , globalComm ) ;
    globalResMultiVector->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( globalThermalResidualVec );
    globalResMultiVector->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( flowResVec );


    boost::shared_ptr<AMP::Operator::OperatorParameters> emptyParams;
    boost::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator (new AMP::Operator::ColumnOperator(emptyParams));
    boost::shared_ptr<AMP::Operator::ColumnOperator> linearColumnOperator (new AMP::Operator::ColumnOperator(emptyParams));
    boost::shared_ptr<AMP::Operator::ColumnOperator> volumeIntegralColumnOperator(new AMP::Operator::ColumnOperator(emptyParams));

    std::vector<AMP::Mesh::MeshID> meshIDs = manager->getBaseMeshIDs();

    // CREATE OPERATORS 
    
    for( size_t meshIndex=0; meshIndex<meshIDs.size(); meshIndex++ )
    {
      AMP::Mesh::Mesh::shared_ptr adapter =  manager->Subset( meshIDs[meshIndex] );
      if( adapter.get() == NULL ) continue;

      std::string meshName = adapter->getName();
      std::string prefix, prefixPower;

      if( meshName.compare("clad")==0 ) {
        prefix="Clad";
        prefixPower="Clad";
      } else if ( meshName.compare("pellet_1")==0 ) {
        prefix="BottomPellet";
        prefixPower="Pellet";
      } else if ( meshName.compare("pellet_3")==0 ) {
        prefix="TopPellet";
        prefixPower="Pellet";
      } else if ( meshName.compare(0,7,"pellet_")==0 ) {
        prefix="MiddlePellet";
        prefixPower="Pellet";
      } else if ( meshName.compare("subchannel")==0 ){
        continue;
      } else {
        AMP_ERROR("Unknown Mesh");
      }

      //-----------------------------------------------
      //   CREATE THE NONLINEAR THERMAL OPERATOR 1 ----
      //-----------------------------------------------
      AMP_INSIST( global_input_db->keyExists(prefix+"NonlinearThermalOperator"), "key missing!" );
      boost::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;
      boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> thermalNonlinearOperator = boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator( adapter,
            prefix+"NonlinearThermalOperator",
            global_input_db,
            thermalTransportModel));
      nonlinearColumnOperator->append(thermalNonlinearOperator);

      //-------------------------------------
      //   CREATE THE LINEAR THERMAL OPERATOR 1 ----
      //-------------------------------------
      AMP_INSIST( global_input_db->keyExists(prefix+"LinearThermalOperator"), "key missing!" );
      boost::shared_ptr<AMP::Operator::LinearBVPOperator> thermalLinearOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(adapter,
            prefix+"LinearThermalOperator",
            global_input_db,
            thermalTransportModel));
      linearColumnOperator->append(thermalLinearOperator);

      AMP_INSIST( global_input_db->keyExists(prefixPower+"VolumeIntegralOperator"), "key missing!" );
      boost::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
      boost::shared_ptr<AMP::Operator::VolumeIntegralOperator> specificPowerGpVecToPowerDensityNodalVecOperator =
        boost::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(AMP::Operator::OperatorBuilder::createOperator( adapter,
              prefixPower+"VolumeIntegralOperator",
              global_input_db,
              stransportModel));
      volumeIntegralColumnOperator->append(specificPowerGpVecToPowerDensityNodalVecOperator);


    }

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    // get subchannel physics model
    // for post processing - need water library to convert enthalpy to temperature...
    boost::shared_ptr<AMP::Database> subchannelPhysics_db = global_input_db->getDatabase("SubchannelPhysicsModel");
    boost::shared_ptr<AMP::Operator::ElementPhysicsModelParameters> params( new AMP::Operator::ElementPhysicsModelParameters(subchannelPhysics_db));
    boost::shared_ptr<AMP::Operator::SubchannelPhysicsModel>  subchannelPhysicsModel (new AMP::Operator::SubchannelPhysicsModel(params));
    boost::shared_ptr<AMP::Operator::SubchannelTwoEqNonlinearOperator> subchannelNonlinearOperator;
    boost::shared_ptr<AMP::Operator::SubchannelTwoEqLinearOperator> subchannelLinearOperator; 

    for( size_t meshIndex=0; meshIndex<meshIDs.size(); meshIndex++ )
    {
      AMP::Mesh::Mesh::shared_ptr adapter =  manager->Subset( meshIDs[meshIndex] );
      if( adapter.get() == NULL ) continue;

      std::string meshName = adapter->getName();
      if ( meshName.compare("subchannel")==0 ){

      //-----------------------------------------------
      //   CREATE THE NONLINEAR THERMAL OPERATOR 1 ----
      //-----------------------------------------------
        boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel;
        subchannelNonlinearOperator = boost::dynamic_pointer_cast<AMP::Operator::SubchannelTwoEqNonlinearOperator>(AMP::Operator::OperatorBuilder::createOperator(
                                                                                  adapter ,"SubchannelTwoEqNonlinearOperator",global_input_db,elementModel ));


        // create linear operator
        subchannelLinearOperator = boost::dynamic_pointer_cast<AMP::Operator::SubchannelTwoEqLinearOperator>(AMP::Operator::OperatorBuilder::createOperator(
                                                                                  adapter ,"SubchannelTwoEqLinearOperator",global_input_db,elementModel ));


        subchannelNonlinearOperator->setVector(subchannelFuelTemp); 

        // pass creation test
        ut->passes(exeName+": creation");
        std::cout.flush();

        nonlinearColumnOperator->append(subchannelNonlinearOperator);
        linearColumnOperator->append(subchannelLinearOperator);
      }
    }


    // CREATE MAPS
    boost::shared_ptr<AMP::Operator::ColumnOperator> mapsColumn( new AMP::Operator::ColumnOperator ( emptyParams) );
    boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator> n2nColumn( new AMP::Operator::AsyncMapColumnOperator ( emptyParams) );
    boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator> szaColumn( new AMP::Operator::AsyncMapColumnOperator ( emptyParams) );

    std::vector<AMP::Mesh::Mesh::shared_ptr> pins = boost::dynamic_pointer_cast<AMP::Mesh::MultiMesh>(manager)->getMeshes();

    for (size_t i=0; i<pins.size(); i++) {
      boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator> map = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::NodeToNodeMap>( pins[i], global_input_db->getDatabase("ThermalNodeToNodeMaps") );
      for (size_t j=0; j<map->getNumberOfOperators(); j++)
        n2nColumn->append( map->getOperator(j) );
      boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator> sza = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::ScalarZAxisMap> ( pins[i],  global_input_db->getDatabase("ThermalScalarZAxisMaps") );
      for (size_t j=0; j<sza->getNumberOfOperators(); j++)
        szaColumn->append( sza->getOperator(j) );
    }
    if ( n2nColumn->getNumberOfOperators() > 0 )
      mapsColumn->append( n2nColumn );
    if ( szaColumn->getNumberOfOperators() > 0 )
      mapsColumn->append( szaColumn );

    n2nColumn->setVector ( thermalMapVec );
    szaColumn->setVector ( thermalMapVec );

    int  curOperator = 0;
    for( size_t meshIndex=0; meshIndex<meshIDs.size(); meshIndex++ )
    {
      AMP::Mesh::Mesh::shared_ptr adapter =  manager->Subset( meshIDs[meshIndex] );
      if( adapter.get() == NULL ) continue;

      std::string meshName = adapter->getName();
      std::string prefix;

      if( meshName.compare("clad")==0 ) {
        prefix="Clad";
      } else if ( meshName.compare("pellet_1")==0 ) {
        prefix="BottomPellet";
      } else if ( meshName.compare("pellet_3")==0 ) {
        prefix="TopPellet";
      } else if ( meshName.compare(0,7,"pellet_")==0 ) {
        prefix="MiddlePellet";
      } else if ( meshName.compare("subchannel")==0 ){
        continue;
      } else {
        AMP_ERROR("Unknown Mesh");
      }

      boost::shared_ptr<AMP::Operator::NonlinearBVPOperator>  curBVPop = boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator> ( nonlinearColumnOperator->getOperator ( curOperator ) );
      boost::shared_ptr<AMP::Operator::ColumnBoundaryOperator>  curBCcol = boost::dynamic_pointer_cast<AMP::Operator::ColumnBoundaryOperator> ( curBVPop->getBoundaryOperator () );
      boost::shared_ptr<AMP::Database> operator_db = global_input_db->getDatabase ( prefix+"NonlinearThermalOperator" );
      boost::shared_ptr<AMP::Database>  curBCdb = global_input_db->getDatabase(operator_db->getString ( "BoundaryOperator" ));
      std::vector<std::string>  opNames = curBCdb->getStringArray ( "boundaryOperators" );
      for ( int curBCentry = 0 ; curBCentry != curBCdb->getInteger ( "numberOfBoundaryOperators" ) ; curBCentry++ )
      {
        if ( opNames[curBCentry] == "P2CRobinVectorCorrection" )
        {
          boost::shared_ptr<AMP::Operator::RobinVectorCorrection>  gapBC = boost::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection> ( curBCcol->getBoundaryOperator ( curBCentry ) );
          AMP_ASSERT(thermalMapVec!=NULL);
          gapBC->setVariableFlux ( thermalMapVec );
          gapBC->reset ( gapBC->getParameters() );
        }
        if ( ( opNames[curBCentry] == "BottomP2PNonlinearRobinVectorCorrection" )|| ( opNames[curBCentry] == "MiddleP2PNonlinearRobinBoundaryCondition" ) || ( opNames[curBCentry] == "TopP2PNonlinearRobinBoundaryCondition" ) )
        {
          boost::shared_ptr<AMP::Operator::RobinVectorCorrection>  p2pBC = boost::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection> ( curBCcol->getBoundaryOperator ( curBCentry ) );
          AMP_ASSERT(thermalMapVec!=NULL);
          p2pBC->setVariableFlux ( thermalMapVec );
          p2pBC->reset ( p2pBC->getParameters() );
        }
        if ( opNames[curBCentry] == "C2WBoundaryVectorCorrection" ) 
        {
          boost::shared_ptr<AMP::Database>  thisDb = global_input_db->getDatabase( opNames[curBCentry] );
          bool isCoupled = thisDb->getBoolWithDefault( "IsCoupledBoundary_0", false);
          if( isCoupled ) 
          {
            boost::shared_ptr<AMP::Operator::RobinVectorCorrection>  p2pBC = boost::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection> ( curBCcol->getBoundaryOperator ( curBCentry ) );
            AMP_ASSERT(thermalMapVec!=NULL);
            p2pBC->setVariableFlux ( thermalMapVec );
            p2pBC->reset ( p2pBC->getParameters() );
          }
        }
        if ( opNames[curBCentry] == "C2PRobinVectorCorrection" )
        {
          boost::shared_ptr<AMP::Operator::RobinVectorCorrection>  gapBC = boost::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection> ( curBCcol->getBoundaryOperator ( curBCentry ) );
          AMP_ASSERT(thermalMapVec!=NULL);
          gapBC->setVariableFlux ( thermalMapVec );
          gapBC->reset ( gapBC->getParameters() );
        }
      }

      curOperator++;
    }

    boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator>  cladToSubchannelMap, subchannelToCladMap; 
    
    boost::shared_ptr<AMP::Database> cladToSubchannelDb = global_input_db->getDatabase( "CladToSubchannelMaps" );
    cladToSubchannelMap = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::CladToSubchannelMap>( manager, cladToSubchannelDb );
    cladToSubchannelMap->setVector( subchannelFuelTemp );
    mapsColumn->append( cladToSubchannelMap );

    // unusual map - take flow variable (enthalpy and pressure) on subchannel mesh and convert to temperature then map to clad surface
    boost::shared_ptr<AMP::Database> subchannelToCladDb = global_input_db->getDatabase( "SubchannelToCladMaps" );
    subchannelToCladMap = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::SubchannelToCladMap>( manager, subchannelToCladDb );
    subchannelToCladMap->setVector( thermalMapVec );
    //
    boost::shared_ptr<AMP::InputDatabase> emptyDb (new AMP::InputDatabase("empty"));
    emptyDb->putInteger("print_info_level",0); 
    boost::shared_ptr<AMP::Operator::CoupledChannelToCladMapOperatorParameters> coupledChannelMapOperatorParams(new AMP::Operator::CoupledChannelToCladMapOperatorParameters( emptyDb ));
    coupledChannelMapOperatorParams->d_variable       = flowVariable;
    coupledChannelMapOperatorParams->d_vector         = subchannelFlowTemp;
    coupledChannelMapOperatorParams->d_mapOperator    = subchannelToCladMap;
    coupledChannelMapOperatorParams->d_subchannelMesh = subchannelMesh;
    coupledChannelMapOperatorParams->d_subchannelPhysicsModel = subchannelPhysicsModel ;
    boost::shared_ptr<AMP::Operator::Operator> coupledChannelMapOperator (new AMP::Operator::CoupledChannelToCladMapOperator(coupledChannelMapOperatorParams));
        
    mapsColumn->append( coupledChannelMapOperator );

    boost::shared_ptr<AMP::InputDatabase> copyOp_db = boost::dynamic_pointer_cast<AMP::InputDatabase>(global_input_db->getDatabase("CopyOperator"));
    boost::shared_ptr<AMP::Operator::VectorCopyOperatorParameters> vecCopyOperatorParams(new AMP::Operator::VectorCopyOperatorParameters( copyOp_db ));
    vecCopyOperatorParams->d_copyVariable = thermalVariable;
    vecCopyOperatorParams->d_copyVector = thermalMapVec;
    vecCopyOperatorParams->d_Mesh = pinMesh ;
    boost::shared_ptr<AMP::Operator::Operator> thermalCopyOperator(new AMP::Operator::VectorCopyOperator(vecCopyOperatorParams));

    boost::shared_ptr<AMP::Operator::CoupledOperatorParameters> CoupledOpParams(new AMP::Operator::CoupledOperatorParameters(emptyDb));
    CoupledOpParams->d_CopyOperator = thermalCopyOperator;
    CoupledOpParams->d_MapOperator = mapsColumn;
    CoupledOpParams->d_BVPOperator = nonlinearColumnOperator;
    boost::shared_ptr<AMP::Operator::Operator> nonlinearCoupledOperator(new AMP::Operator::CoupledOperator(CoupledOpParams));

    // get exit pressure
    double Pout = global_input_db->getDatabase( "SubchannelTwoEqNonlinearOperator" )->getDouble("Exit_Pressure");
    // calculate initial guess
    // get inlet temperature
    double Tin = global_input_db->getDatabase( "SubchannelTwoEqNonlinearOperator" )->getDouble("Inlet_Temperature");
    // compute inlet enthalpy
    std::map<std::string, boost::shared_ptr<std::vector<double> > > enthalpyArgMap;
    enthalpyArgMap.insert(std::make_pair("temperature",new std::vector<double>(1,Tin)));
    enthalpyArgMap.insert(std::make_pair("pressure",   new std::vector<double>(1,Pout)));
    std::vector<double> enthalpyResult(1);
    subchannelPhysicsModel->getProperty("Enthalpy",enthalpyResult,enthalpyArgMap); 
    double hin = enthalpyResult[0];
    std::cout<< "Enthalpy Solution:"<< hin <<std::endl;

    AMP::LinearAlgebra::Vector::shared_ptr subchannelEnthalpy = flowSolVec->select( AMP::LinearAlgebra::VS_Stride("H", 0, 2) , "H" );
    AMP::LinearAlgebra::Vector::shared_ptr subchannelPressure = flowSolVec->select( AMP::LinearAlgebra::VS_Stride("P", 1, 2) , "P" );

    subchannelEnthalpy->setToScalar(hin); 
    subchannelPressure->setToScalar(Pout); 

    // FIRST APPLY CALL
    subchannelLinearOperator->reset(subchannelNonlinearOperator->getJacobianParameters(flowSolVec));
    subchannelLinearOperator->apply( flowRhsVec, flowSolVec, flowResVec, 1.0, -1.0);

    thermalMapVec->setToScalar(400.);
    globalThermalSolutionVec->setToScalar(400.);
    nonlinearCoupledOperator->apply(globalRhsMultiVector, globalSolMultiVector, globalResMultiVector, 1.0, -1.0);

    // get nonlinear solver database
    boost::shared_ptr<AMP::Database> nonlinearSolver_db = global_input_db->getDatabase("NonlinearSolver"); 
    boost::shared_ptr<AMP::Database> linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver"); 

    // create nonlinear solver parameters
    boost::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(new AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db));
    nonlinearSolverParams->d_comm = globalComm;
    nonlinearSolverParams->d_pOperator = nonlinearCoupledOperator;
    nonlinearSolverParams->d_pInitialGuess = globalSolMultiVector;
    boost::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));

    // create preconditioner
    boost::shared_ptr<AMP::Database> columnPreconditioner_db = linearSolver_db->getDatabase("Preconditioner");
    boost::shared_ptr<AMP::Solver::SolverStrategyParameters> columnPreconditionerParams(new AMP::Solver::SolverStrategyParameters(columnPreconditioner_db));
    columnPreconditionerParams->d_pOperator = linearColumnOperator;
    columnPreconditioner.reset(new AMP::Solver::ColumnSolver(columnPreconditionerParams));

    boost::shared_ptr<AMP::Database> trilinosPreconditioner_db = columnPreconditioner_db->getDatabase("TrilinosPreconditioner");
    for(unsigned int id = 0; id != linearColumnOperator->getNumberOfOperators(); id++)
    {
      boost::shared_ptr<AMP::Solver::SolverStrategyParameters> trilinosPreconditionerParams(new AMP::Solver::SolverStrategyParameters(trilinosPreconditioner_db));
      trilinosPreconditionerParams->d_pOperator = linearColumnOperator->getOperator(id);
      boost::shared_ptr<AMP::Solver::TrilinosMLSolver> trilinosPreconditioner(new AMP::Solver::TrilinosMLSolver(trilinosPreconditionerParams));
      columnPreconditioner->append(trilinosPreconditioner);
    }

    // create linear solver
    boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver = nonlinearSolver->getKrylovSolver();
    // set preconditioner
    linearSolver->setPreconditioner(columnPreconditioner);

    // don't use zero initial guess
    nonlinearSolver->setZeroInitialGuess(false);

    AMP::Discretization::DOFManager::shared_ptr  gaussPtDOFManager;
    if ( pinMesh.get()!=NULL ) {
        gaussPtDOFManager = AMP::Discretization::simpleDOFManager::create(pinMesh,AMP::Mesh::Volume,1,8);
        specificPowerGpVec= AMP::LinearAlgebra::createVector( gaussPtDOFManager , powerVariable );
        specificPowerGpVec->setToScalar(0.0);
    }

    // Initialize the pin temperatures
    if ( pinMesh.get()!=NULL ) {
        AMP::Mesh::MeshIterator it = pinMesh->getIterator(AMP::Mesh::Volume,0);
        std::vector<size_t> dofs;
        for (size_t i=0; i<it.size(); i++) {
            gaussPtDOFManager->getDOFs(it->globalID(),dofs);
            for (size_t j=0; j<dofs.size(); j++) {
              specificPowerGpVec->setValueByGlobalID(dofs[j],getPower(it->centroid()));
            }
            ++it;
        }
    }

    specificPowerGpVec->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);

    volumeIntegralColumnOperator->apply(nullVec, specificPowerGpVec, globalThermalRhsVec , 1., 0.);
    size_t totalOp;
    if(subchannelMesh != NULL ){
      totalOp = nonlinearColumnOperator->getNumberOfOperators()-1;
    }else{
      totalOp = nonlinearColumnOperator->getNumberOfOperators();
    }
    for( size_t id = 0; id != totalOp; id++)
    {
      boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearThermalOperator = boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>( nonlinearColumnOperator->getOperator(id));
      nonlinearThermalOperator->modifyInitialSolutionVector(globalThermalSolutionVec);
      nonlinearThermalOperator->modifyRHSvector(globalThermalRhsVec);
    }
    // solve
    nonlinearSolver->solve(globalRhsMultiVector, globalSolMultiVector);

    AMP::LinearAlgebra::Vector::shared_ptr flowDenVec  = subchannelFuelTemp->cloneVector(); 
    AMP::LinearAlgebra::Vector::shared_ptr flowTempVec = subchannelFuelTemp->cloneVector(); 

    AMP::Mesh::MeshIterator face  = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
    AMP::Mesh::MeshIterator end_face = face.end();
    std::vector<size_t> dofs;
    std::vector<size_t> scalarDofs;
    for( ; face != end_face; ++face){
      faceDOFManager->getDOFs( face->globalID(), dofs );
      scalarFaceDOFManager->getDOFs( face->globalID(), scalarDofs );
      std::map<std::string, boost::shared_ptr<std::vector<double> > > outTemperatureArgMap;
      outTemperatureArgMap.insert(std::make_pair("enthalpy",new std::vector<double>(1,flowSolVec->getValueByGlobalID(dofs[0]))));
      outTemperatureArgMap.insert(std::make_pair("pressure",new std::vector<double>(1,flowSolVec->getValueByGlobalID(dofs[1]))));
      std::vector<double> outTemperatureResult(1);
      subchannelPhysicsModel->getProperty("Temperature", outTemperatureResult, outTemperatureArgMap); 
      flowTempVec->setValueByGlobalID(scalarDofs[0] ,outTemperatureResult[0]); 
     } 

#ifdef USE_SILO
    // Register the quantities to plot
    AMP::Mesh::SiloIO::shared_ptr  siloWriter( new AMP::Mesh::SiloIO );
    siloWriter->registerVector( flowSolVec, xyFaceMesh, AMP::Mesh::Face, "SubchannelFlow" );
    siloWriter->registerVector( flowTempVec, xyFaceMesh, AMP::Mesh::Face, "FlowTemp" );
    siloWriter->registerVector( globalThermalSolutionVec ,  pinMesh , AMP::Mesh::Vertex, "Temperature" );
    siloWriter->registerVector( specificPowerGpVec,  pinMesh , AMP::Mesh::Volume, "Power" );
    siloWriter->writeFile( silo_name , 0 );
#endif

}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  try {
    SubchannelSolve(&ut, "testSubchannelSolve");
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


