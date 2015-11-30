#include <string>
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"

#include "utils/shared_ptr.h"

#include "utils/InputDatabase.h"
#include "utils/Utilities.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/Database.h"
#include "ProfilerApp.h"

#include "utils/Writer.h"
#include "ampmesh/Mesh.h"
#include "ampmesh/MultiMesh.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"

#include "vectors/Variable.h"
#include "vectors/VectorBuilder.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/SimpleVector.h"
#include "vectors/VectorSelector.h"

#include "operators/boundary/ColumnBoundaryOperator.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/boundary/libmesh/RobinVectorCorrection.h"
#include "operators/boundary/libmesh/NeumannVectorCorrectionParameters.h"
#include "operators/boundary/DirichletMatrixCorrection.h"
#include "operators/boundary/libmesh/NeumannVectorCorrection.h"

#include "operators/VectorCopyOperator.h"

#include "operators/IdentityOperator.h"
#include "operators/libmesh/VolumeIntegralOperator.h"
#include "operators/CoupledOperator.h"
#include "operators/CoupledOperatorParameters.h"
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
#include "operators/subchannel/SubchannelHelpers.h"
#include "operators/subchannel/CoupledChannelToCladMapOperator.h"
#include "operators/subchannel/SubchannelConstants.h"
#include "operators/subchannel/SubchannelToPointMap.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"

#include "solvers/petsc/PetscKrylovSolver.h"
#include "solvers/petsc/PetscSNESSolver.h"
#include "solvers/trilinos/TrilinosMLSolver.h"
#include "solvers/ColumnSolver.h"
#include "solvers/BandedSolver.h"

#include "ampmesh/StructuredMeshHelper.h"
#include "discretization/structuredFaceDOFManager.h"




// Function to get an arbitrary power profile (W/kg) assuming a density of 1 kg/m^3 for the volume integral
// P is the total power of the pin, V is the volume of the pin
double getPower( const std::vector<double> &range, double P, double V, const std::vector<double> &pos ) {
    const double pi = 3.1415926535897932;
    double x = (pos[0]-range[0])/(range[1]-range[0]);
    double y = (pos[1]-range[2])/(range[3]-range[2]);
    double z = (pos[2]-range[4])/(range[5]-range[4]);
    return P/V*(0.8+0.2*x+0.2*y)*pi/2*sin(pi*z);
}



// Function to create the solution vectors
void createVectors( AMP::Mesh::Mesh::shared_ptr pinMesh, AMP::Mesh::Mesh::shared_ptr subchannelMesh, 
    AMP::LinearAlgebra::Vector::shared_ptr &globalMultiVector, 
    AMP::LinearAlgebra::Vector::shared_ptr &specificPowerGpVec )
{
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
    globalMultiVector = AMP::LinearAlgebra::MultiVector::create( "multivector" , globalComm ) ;

    AMP::LinearAlgebra::Variable::shared_ptr thermalVariable (new AMP::LinearAlgebra::Variable("Temperature"));    // temperature on pellets and cladding
    AMP::LinearAlgebra::Variable::shared_ptr flowVariable (new AMP::LinearAlgebra::Variable("Flow"));              // enthalpy and pressure on z-faces
    AMP::LinearAlgebra::Variable::shared_ptr powerVariable (new AMP::LinearAlgebra::Variable("SpecificPowerInWattsPerGram")); // specific power on gauss points

    AMP::LinearAlgebra::Vector::shared_ptr thermalVec;
    if ( pinMesh.get()!=NULL ) {
        AMP::Discretization::DOFManager::shared_ptr nodalScalarDOF = 
            AMP::Discretization::simpleDOFManager::create( pinMesh ,AMP::Mesh::Vertex,1,1,true);
        thermalVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, thermalVariable );
    }
    globalMultiVector->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( thermalVec );

    AMP::LinearAlgebra::Vector::shared_ptr flowVec;
    if ( subchannelMesh.get()!=NULL ) {
        int DOFsPerFace[3]={0,0,2};
        AMP::Discretization::DOFManager::shared_ptr faceDOFManager = 
            AMP::Discretization::structuredFaceDOFManager::create( subchannelMesh, DOFsPerFace, 0 );
        // create solution, rhs, and residual vectors
        flowVec = AMP::LinearAlgebra::createVector( faceDOFManager, flowVariable, true );
    }
    globalMultiVector->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( flowVec );

    if ( pinMesh.get()!=NULL ) {
        AMP::Discretization::DOFManager::shared_ptr gaussPtDOFManager = 
            AMP::Discretization::simpleDOFManager::create(pinMesh,AMP::Mesh::Volume,1,8);
        specificPowerGpVec= AMP::LinearAlgebra::createVector( gaussPtDOFManager , powerVariable );
        specificPowerGpVec->setToScalar(0.0);
    }
}


void SubchannelSolve(AMP::UnitTest *ut, std::string exeName )
{
    PROFILE_START("Main");
    std::string input_file = "input_" + exeName;
    std::string log_file = "output_" + exeName;
    AMP::PIO::logAllNodes(log_file);
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
    globalComm.barrier();

    // Read the input file
    AMP::shared_ptr<AMP::InputDatabase>  global_input_db ( new AMP::InputDatabase ( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile ( input_file , global_input_db );
    global_input_db->printClassData(AMP::plog);

    // Get the Mesh database and create the mesh parameters
    AMP::shared_ptr<AMP::Database> database = global_input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(database));
    meshParams->setComm(globalComm);

    // Get the meshes
    AMP::Mesh::Mesh::shared_ptr manager = AMP::Mesh::Mesh::buildMesh(meshParams);
    AMP::Mesh::Mesh::shared_ptr pinMesh = manager->Subset("MultiPin");
    AMP::Mesh::Mesh::shared_ptr cladMesh;
    if ( pinMesh.get()!=NULL ) {
        pinMesh->setName("MultiPin");
        cladMesh = pinMesh->Subset("clad");
    }
    AMP::Mesh::Mesh::shared_ptr subchannelMesh = manager->Subset("subchannel");
    AMP::Mesh::Mesh::shared_ptr xyFaceMesh;
    if ( subchannelMesh.get()!=NULL ) {
        AMP::Mesh::MeshIterator face = AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,0);
        xyFaceMesh = subchannelMesh->Subset( face );
    }

    // Variables
    AMP::LinearAlgebra::Variable::shared_ptr thermalVariable(new AMP::LinearAlgebra::Variable("Temperature"));    // temperature on pellets and cladding
    AMP::LinearAlgebra::Variable::shared_ptr flowVariable(new AMP::LinearAlgebra::Variable("Flow"));              // enthalpy and pressure on z-faces
    AMP::LinearAlgebra::Variable::shared_ptr powerVariable(new AMP::LinearAlgebra::Variable("SpecificPowerInWattsPerGram")); // specific power on gauss points

    AMP::shared_ptr<AMP::Operator::Operator> thermalCopyOperator;

    AMP::shared_ptr<AMP::Operator::OperatorParameters> emptyParams;
    AMP::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator (new AMP::Operator::ColumnOperator(emptyParams));
    AMP::shared_ptr<AMP::Operator::ColumnOperator> linearColumnOperator (new AMP::Operator::ColumnOperator(emptyParams));
    AMP::shared_ptr<AMP::Operator::ColumnOperator> volumeIntegralColumnOperator(new AMP::Operator::ColumnOperator(emptyParams));

    AMP::shared_ptr<AMP::Operator::ColumnOperator> mapsColumn( new AMP::Operator::ColumnOperator ( emptyParams) );
    AMP::shared_ptr<AMP::Operator::AsyncMapColumnOperator> n2nColumn( new AMP::Operator::AsyncMapColumnOperator ( emptyParams) );
    AMP::shared_ptr<AMP::Operator::AsyncMapColumnOperator> szaColumn( new AMP::Operator::AsyncMapColumnOperator ( emptyParams) );

    AMP::shared_ptr<AMP::Solver::PetscSNESSolver>  nonlinearCoupledSolver;
    AMP::shared_ptr<AMP::Solver::PetscKrylovSolver>  linearColumnSolver;
    AMP::shared_ptr<AMP::Solver::ColumnSolver>  columnPreconditioner;

    if ( pinMesh.get()!=NULL ) {
        std::vector<AMP::Mesh::MeshID> pinMeshIDs = pinMesh->getBaseMeshIDs();

        // CREATE OPERATORS 
        for( size_t meshIndex=0; meshIndex < pinMeshIDs.size(); meshIndex++ ) {
            AMP::Mesh::Mesh::shared_ptr adapter =  manager->Subset( pinMeshIDs[meshIndex] );
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
            } else {
                AMP_ERROR("Unknown Mesh");
            }

            // CREATE THE NONLINEAR THERMAL OPERATOR 1
            AMP_INSIST( global_input_db->keyExists(prefix+"NonlinearThermalOperator"), "key missing!" );
            AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;
            AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> thermalNonlinearOperator = 
                AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
                AMP::Operator::OperatorBuilder::createOperator( adapter,
                prefix+"NonlinearThermalOperator",
                global_input_db,
                thermalTransportModel));
            nonlinearColumnOperator->append(thermalNonlinearOperator);

            // CREATE THE LINEAR THERMAL OPERATOR 1 
            AMP_INSIST( global_input_db->keyExists(prefix+"LinearThermalOperator"), "key missing!" );
            AMP::shared_ptr<AMP::Operator::LinearBVPOperator> thermalLinearOperator = 
                AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
                AMP::Operator::OperatorBuilder::createOperator(adapter,
                prefix+"LinearThermalOperator",
                global_input_db,
                thermalTransportModel));
            linearColumnOperator->append(thermalLinearOperator);

            AMP_INSIST( global_input_db->keyExists(prefixPower+"VolumeIntegralOperator"), "key missing!" );
            AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
            AMP::shared_ptr<AMP::Operator::VolumeIntegralOperator> specificPowerGpVecToPowerDensityNodalVecOperator =
            AMP::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
                AMP::Operator::OperatorBuilder::createOperator( adapter,
                  prefixPower+"VolumeIntegralOperator",
                  global_input_db,
                  stransportModel));
            volumeIntegralColumnOperator->append(specificPowerGpVecToPowerDensityNodalVecOperator);


        }

    }

    // Get the subchannel hydraulic diameter for the temperature boundary condition 
    AMP::LinearAlgebra::Vector::shared_ptr ChannelDiameterVec = 
        AMP::Operator::Subchannel::getCladHydraulicDiameter( cladMesh, subchannelMesh, globalComm );


    AMP::LinearAlgebra::Vector::shared_ptr subchannelFuelTemp; 
    AMP::LinearAlgebra::Vector::shared_ptr subchannelFlowTemp;
    if ( subchannelMesh.get()!=NULL ) {
        int DOFsPerFace[3]={0,0,1};
        AMP::Discretization::DOFManager::shared_ptr scalarFaceDOFManager = 
            AMP::Discretization::structuredFaceDOFManager::create( subchannelMesh, DOFsPerFace, 0 );
        subchannelFuelTemp = AMP::LinearAlgebra::createVector( scalarFaceDOFManager, thermalVariable );
        subchannelFlowTemp = AMP::LinearAlgebra::createVector( scalarFaceDOFManager, thermalVariable );
    }

    // get subchannel physics model
    // for post processing - need water library to convert enthalpy to temperature...
    AMP::shared_ptr<AMP::Database> subchannelPhysics_db = global_input_db->getDatabase("SubchannelPhysicsModel");
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModelParameters> params( new AMP::Operator::ElementPhysicsModelParameters(subchannelPhysics_db));
    AMP::shared_ptr<AMP::Operator::SubchannelPhysicsModel>  subchannelPhysicsModel (new AMP::Operator::SubchannelPhysicsModel(params));
    AMP::shared_ptr<AMP::Operator::SubchannelTwoEqNonlinearOperator> subchannelNonlinearOperator;
    AMP::shared_ptr<AMP::Operator::LinearOperator> subchannelLinearOperator; 

    // Get the subchannel operators
    std::vector<double> clad_x, clad_y, clad_d;
    AMP::Operator::Subchannel::getCladProperties( globalComm, cladMesh, clad_x, clad_y,clad_d );
    if ( subchannelMesh.get()!=NULL ) {
        std::vector<AMP::Mesh::MeshID> subChannelMeshIDs = subchannelMesh->getBaseMeshIDs();

        for( size_t meshIndex=0; meshIndex < subChannelMeshIDs.size(); meshIndex++ ) {
            AMP::Mesh::Mesh::shared_ptr adapter =  manager->Subset( subChannelMeshIDs[meshIndex] );
            if( adapter.get() == NULL ) continue;

            std::string meshName = adapter->getName();
            if ( meshName.compare("subchannel")==0 ){
                // create the non-linear operator
                subchannelNonlinearOperator = AMP::dynamic_pointer_cast<AMP::Operator::SubchannelTwoEqNonlinearOperator>(
                    AMP::Operator::OperatorBuilder::createOperator( adapter ,"SubchannelTwoEqNonlinearOperator",global_input_db ) );
                AMP::shared_ptr<AMP::Operator::SubchannelOperatorParameters> nonlinearOpParams = subchannelNonlinearOperator->getParams();
                nonlinearOpParams->clad_x = clad_x;
                nonlinearOpParams->clad_y = clad_y;
                nonlinearOpParams->clad_d = clad_d;
                subchannelNonlinearOperator->reset(nonlinearOpParams);
                // create linear operator
                subchannelLinearOperator = AMP::dynamic_pointer_cast<AMP::Operator::SubchannelTwoEqLinearOperator>(
                    AMP::Operator::OperatorBuilder::createOperator( adapter ,"SubchannelTwoEqLinearOperator",
                    global_input_db, subchannelNonlinearOperator->getSubchannelPhysicsModel() ) );
                //subchannelLinearOperator.reset( new AMP::Operator::IdentityOperator( nonlinearOpParams ) );
                int DOFsPerFace[3]={0,0,2};
                AMP::Discretization::DOFManager::shared_ptr flowDOFManager = 
                    AMP::Discretization::structuredFaceDOFManager::create( subchannelMesh, DOFsPerFace, 0 );
                AMP::LinearAlgebra::Vector::shared_ptr subchannelFlow = 
                    AMP::LinearAlgebra::createVector( flowDOFManager, flowVariable );
                subchannelNonlinearOperator->setVector(subchannelFuelTemp); 
                AMP::shared_ptr<AMP::Operator::SubchannelOperatorParameters> subchannelLinearParams = 
                    AMP::dynamic_pointer_cast<AMP::Operator::SubchannelOperatorParameters>( 
                    subchannelNonlinearOperator->getParameters("Jacobian", subchannelFlow) );
                subchannelLinearParams->d_initialize = false;
                subchannelLinearOperator->reset(subchannelLinearParams);
                // pass creation test
                ut->passes(exeName+": creation");
                std::cout.flush();
                nonlinearColumnOperator->append(subchannelNonlinearOperator);
                // Do not add the subchannel to the linear operator (we will add it later)
            }
        }
    }


    // CREATE MAPS
    AMP::LinearAlgebra::Vector::shared_ptr thermalMapVec;
    AMP::LinearAlgebra::Vector::shared_ptr density_map_vec;
    if ( cladMesh.get()!=NULL ) {
        AMP::Discretization::DOFManager::shared_ptr nodalScalarDOF = 
            AMP::Discretization::simpleDOFManager::create(cladMesh,AMP::Mesh::Vertex,1,1,true);
        AMP::LinearAlgebra::Variable::shared_ptr densityVariable( new AMP::LinearAlgebra::Variable("Density") );
        density_map_vec = AMP::LinearAlgebra::createVector( nodalScalarDOF, densityVariable );
        density_map_vec->zero();
    }
    if ( pinMesh.get()!=NULL ) {
        // flow temperature on clad outer surfaces and pellet temperature on clad innner surfaces, and clad inner surface temp on pellet outer surfaces
        AMP::Discretization::DOFManager::shared_ptr nodalScalarDOF = 
            AMP::Discretization::simpleDOFManager::create( pinMesh ,AMP::Mesh::Vertex,1,1,true);
        thermalMapVec = AMP::LinearAlgebra::createVector(nodalScalarDOF, thermalVariable , true);
 
        std::vector<AMP::Mesh::MeshID> pinMeshIDs = pinMesh->getBaseMeshIDs();

        std::vector<AMP::Mesh::Mesh::shared_ptr> pins = AMP::dynamic_pointer_cast<AMP::Mesh::MultiMesh>(pinMesh)->getMeshes();

        for (size_t i=0; i<pins.size(); i++) {
            if ( global_input_db->keyExists("ThermalNodeToNodeMaps") ) {
                AMP::shared_ptr<AMP::Operator::AsyncMapColumnOperator> map = 
                    AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::NodeToNodeMap>( 
                    pins[i], global_input_db->getDatabase("ThermalNodeToNodeMaps") );
                for (size_t j=0; j<map->getNumberOfOperators(); j++)
                    n2nColumn->append( map->getOperator(j) );
            }
            if ( global_input_db->keyExists("ThermalScalarZAxisMaps") ) {
                AMP::shared_ptr<AMP::Operator::AsyncMapColumnOperator> sza = 
                    AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::ScalarZAxisMap>( 
                    pins[i],  global_input_db->getDatabase("ThermalScalarZAxisMaps") );
                for (size_t j=0; j<sza->getNumberOfOperators(); j++)
                    szaColumn->append( sza->getOperator(j) );
            }
        }
        if ( n2nColumn->getNumberOfOperators() > 0 )
            mapsColumn->append( n2nColumn );
        if ( szaColumn->getNumberOfOperators() > 0 )
            mapsColumn->append( szaColumn );

        n2nColumn->setVector ( thermalMapVec );
        szaColumn->setVector ( thermalMapVec );

        int curOperator = 0;
        for( size_t meshIndex=0; meshIndex < pinMeshIDs.size(); meshIndex++ ) {
            AMP::Mesh::Mesh::shared_ptr adapter =  manager->Subset( pinMeshIDs[meshIndex] );
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
            } else {
                AMP_ERROR("Unknown Mesh");
            }

            AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator>  curBVPop = 
                AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>( nonlinearColumnOperator->getOperator ( curOperator ) );
            AMP::shared_ptr<AMP::Operator::ColumnBoundaryOperator>  curBCcol = 
                AMP::dynamic_pointer_cast<AMP::Operator::ColumnBoundaryOperator> ( curBVPop->getBoundaryOperator () );
            AMP::shared_ptr<AMP::Database> operator_db = global_input_db->getDatabase ( prefix+"NonlinearThermalOperator" );
            AMP::shared_ptr<AMP::Database>  curBCdb = global_input_db->getDatabase(operator_db->getString ( "BoundaryOperator" ));
            std::vector<std::string>  opNames = curBCdb->getStringArray( "boundaryOperators" );
            int  opNumber = curBCdb->getInteger( "numberOfBoundaryOperators" );
            for (int curBCentry=0; curBCentry!=opNumber; curBCentry++ ) {
                if ( opNames[curBCentry] == "P2CRobinVectorCorrection" ) {
                    AMP::shared_ptr<AMP::Operator::RobinVectorCorrection>  gapBC = 
                        AMP::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection> ( curBCcol->getBoundaryOperator ( curBCentry ) );
                    AMP_ASSERT(thermalMapVec!=NULL);
                    gapBC->setVariableFlux ( thermalMapVec );
                    gapBC->reset ( gapBC->getOperatorParameters() );
                } else if ( ( opNames[curBCentry]=="BottomP2PNonlinearRobinVectorCorrection"  ) || 
                            ( opNames[curBCentry]=="MiddleP2PNonlinearRobinBoundaryCondition" ) || 
                            ( opNames[curBCentry]=="TopP2PNonlinearRobinBoundaryCondition"    ) ) {
                    AMP::shared_ptr<AMP::Operator::RobinVectorCorrection>  p2pBC = 
                        AMP::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection> ( curBCcol->getBoundaryOperator ( curBCentry ) );
                    AMP_ASSERT(thermalMapVec!=NULL);
                    p2pBC->setVariableFlux ( thermalMapVec );
                    p2pBC->reset ( p2pBC->getOperatorParameters() );
                } else if ( opNames[curBCentry] == "C2WBoundaryVectorCorrection" ) {
                    AMP::shared_ptr<AMP::Database>  thisDb = global_input_db->getDatabase( opNames[curBCentry] );
                    bool isCoupled = thisDb->getBoolWithDefault( "IsCoupledBoundary_0", false);
                    if( isCoupled ) {
                        AMP::shared_ptr<AMP::Operator::RobinVectorCorrection>  c2wBC = 
                            AMP::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection>( curBCcol->getBoundaryOperator ( curBCentry ) );
                        AMP_ASSERT(thermalMapVec!=NULL);
                        c2wBC->setVariableFlux ( thermalMapVec );
                        c2wBC->setFrozenVector( density_map_vec );
                        c2wBC->setFrozenVector( ChannelDiameterVec );
                        c2wBC->reset ( c2wBC->getOperatorParameters() );
                    }
                } else if ( opNames[curBCentry] == "C2PRobinVectorCorrection" ) {
                    AMP::shared_ptr<AMP::Operator::RobinVectorCorrection>  gapBC = 
                        AMP::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection>( curBCcol->getBoundaryOperator ( curBCentry ) );
                    AMP_ASSERT(thermalMapVec!=NULL);
                    gapBC->setVariableFlux ( thermalMapVec );
                    gapBC->reset ( gapBC->getOperatorParameters() );
                } else {
                    AMP_ERROR("Unknown boundary operator");
                }
            }
            curOperator++;
        }
    }

    // Create the maps from the clad temperature to the subchannel temperature
    AMP::shared_ptr<AMP::Operator::AsyncMapColumnOperator>  cladToSubchannelMap; 
    AMP::shared_ptr<AMP::Database> cladToSubchannelDb = global_input_db->getDatabase( "CladToSubchannelMaps" );
    cladToSubchannelMap = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::CladToSubchannelMap>( manager, cladToSubchannelDb );
    cladToSubchannelMap->setVector( subchannelFuelTemp );
    mapsColumn->append( cladToSubchannelMap );

    // Create the maps from the flow variable (enthalpy and pressure) on subchannel mesh and convert to temperature and density then map to clad surface
    AMP::shared_ptr<AMP::Operator::AsyncMapColumnOperator>  thermalSubchannelToCladMap, densitySubchannelToCladMap; 
    AMP::shared_ptr<AMP::Database> thermalCladToSubchannelDb = global_input_db->getDatabase( "ThermalSubchannelToCladMaps" );
    thermalSubchannelToCladMap = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::SubchannelToCladMap>( manager, thermalCladToSubchannelDb );
    AMP::shared_ptr<AMP::Database> densityCladToSubchannelDb = global_input_db->getDatabase( "DensitySubchannelToCladMaps" );
    densitySubchannelToCladMap = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::SubchannelToCladMap>( manager, densityCladToSubchannelDb );
    if ( cladMesh.get()!=NULL ) {
        AMP::LinearAlgebra::VS_Comm commSelector(cladMesh->getComm());
        AMP::LinearAlgebra::Vector::shared_ptr subsetTheramlVec = 
            thermalMapVec->select( commSelector, thermalMapVec->getVariable()->getName() );
        thermalSubchannelToCladMap->setVector( subsetTheramlVec );
        densitySubchannelToCladMap->setVector( density_map_vec );
    }
    AMP::shared_ptr<AMP::InputDatabase> emptyDb (new AMP::InputDatabase("empty"));
    emptyDb->putInteger("print_info_level",0); 
    AMP::shared_ptr<AMP::Operator::CoupledChannelToCladMapOperatorParameters> coupledChannelMapOperatorParams(
        new AMP::Operator::CoupledChannelToCladMapOperatorParameters( emptyDb ) );
    coupledChannelMapOperatorParams->d_variable       = flowVariable;
    coupledChannelMapOperatorParams->d_vector         = subchannelFlowTemp;
    coupledChannelMapOperatorParams->d_Mesh           = subchannelMesh;
    coupledChannelMapOperatorParams->d_thermalMapOperator = thermalSubchannelToCladMap;
    coupledChannelMapOperatorParams->d_densityMapOperator = densitySubchannelToCladMap;
    coupledChannelMapOperatorParams->d_subchannelMesh = subchannelMesh;
    coupledChannelMapOperatorParams->d_subchannelPhysicsModel = subchannelPhysicsModel ;
    AMP::shared_ptr<AMP::Operator::Operator> coupledChannelMapOperator(
        new AMP::Operator::CoupledChannelToCladMapOperator(coupledChannelMapOperatorParams) );
    mapsColumn->append( coupledChannelMapOperator );

    if ( pinMesh.get()!=NULL ) {
        AMP::shared_ptr<AMP::InputDatabase> copyOp_db = AMP::dynamic_pointer_cast<AMP::InputDatabase>(
            global_input_db->getDatabase("CopyOperator") );
        AMP::shared_ptr<AMP::Operator::VectorCopyOperatorParameters> vecCopyOperatorParams(\
            new AMP::Operator::VectorCopyOperatorParameters( copyOp_db ) );
        vecCopyOperatorParams->d_copyVariable = thermalVariable;
        vecCopyOperatorParams->d_copyVector = thermalMapVec;
        vecCopyOperatorParams->d_Mesh = pinMesh ;
        thermalCopyOperator.reset(new AMP::Operator::VectorCopyOperator(vecCopyOperatorParams));
        thermalMapVec->zero();
    }

    AMP::shared_ptr<AMP::Operator::CoupledOperatorParameters> CoupledOpParams(new AMP::Operator::CoupledOperatorParameters(emptyDb));
    CoupledOpParams->d_CopyOperator = thermalCopyOperator;
    CoupledOpParams->d_MapOperator = mapsColumn;
    CoupledOpParams->d_BVPOperator = nonlinearColumnOperator;
    AMP::shared_ptr<AMP::Operator::Operator> nonlinearCoupledOperator(new AMP::Operator::CoupledOperator(CoupledOpParams));


    // Create the solution vector
    AMP::LinearAlgebra::Vector::shared_ptr globalSolMultiVector;
    AMP::LinearAlgebra::Vector::shared_ptr specificPowerGpVec;
    createVectors( pinMesh, subchannelMesh, globalSolMultiVector, specificPowerGpVec );
    
    // Create the rhs and res vectors
    AMP::LinearAlgebra::Vector::shared_ptr globalRhsMultiVector = globalSolMultiVector->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr globalResMultiVector = globalSolMultiVector->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr flowSolVec = globalSolMultiVector->subsetVectorForVariable(flowVariable);
    AMP::LinearAlgebra::Vector::shared_ptr flowRhsVec = globalRhsMultiVector->subsetVectorForVariable(flowVariable);
    AMP::LinearAlgebra::Vector::shared_ptr flowResVec = globalResMultiVector->subsetVectorForVariable(flowVariable);
    AMP::LinearAlgebra::Vector::shared_ptr globalThermalSolVec = globalSolMultiVector->subsetVectorForVariable(thermalVariable);
    AMP::LinearAlgebra::Vector::shared_ptr globalThermalRhsVec = globalRhsMultiVector->subsetVectorForVariable(thermalVariable);
    AMP::LinearAlgebra::Vector::shared_ptr globalThermalResVec = globalResMultiVector->subsetVectorForVariable(thermalVariable);

    // create nonlinear solver
/*    AMP::shared_ptr<AMP::Solver::SolverStrategy> nonlinearSolver;
    { // Limit the scope so we can add an if else statement for Petsc vs NOX

        // get the solver databases
        AMP::shared_ptr<AMP::Database> nonlinearSolver_db = global_input_db->getDatabase("NonlinearSolver"); 
        AMP::shared_ptr<AMP::Database> linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver"); 

        // create preconditioner (thermal domains)
        AMP::shared_ptr<AMP::Database> columnPreconditioner_db = linearSolver_db->getDatabase("Preconditioner");
        AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> columnPreconditionerParams(
            new AMP::Solver::SolverStrategyParameters(columnPreconditioner_db));
        columnPreconditionerParams->d_pOperator = linearColumnOperator;
        columnPreconditioner.reset(new AMP::Solver::ColumnSolver(columnPreconditionerParams));

        AMP::shared_ptr<AMP::Database> trilinosPreconditioner_db = columnPreconditioner_db->getDatabase("TrilinosPreconditioner");
        unsigned int N_preconditioners = linearColumnOperator->getNumberOfOperators();
        for(unsigned int id=0; id<N_preconditioners; id++) {
            AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> trilinosPreconditionerParams(
                new AMP::Solver::SolverStrategyParameters(trilinosPreconditioner_db) );
            trilinosPreconditionerParams->d_pOperator = linearColumnOperator->getOperator(id);
            AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> trilinosPreconditioner(
                new AMP::Solver::TrilinosMLSolver(trilinosPreconditionerParams) );
            columnPreconditioner->append(trilinosPreconditioner);
        }

        // Create the subchannel preconditioner
        if ( subchannelLinearOperator != NULL ) {
            AMP::shared_ptr<AMP::Database> subchannelPreconditioner_db = 
                columnPreconditioner_db->getDatabase("SubchannelPreconditioner");
            AMP_ASSERT(subchannelPreconditioner_db!=NULL);
            AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> subchannelPreconditionerParams(
                new AMP::Solver::SolverStrategyParameters(subchannelPreconditioner_db) );
            subchannelPreconditionerParams->d_pOperator = subchannelLinearOperator;
            std::string preconditioner = subchannelPreconditioner_db->getString("Type");
            if ( preconditioner=="ML" ) {
                AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> subchannelPreconditioner(
                    new AMP::Solver::TrilinosMLSolver(subchannelPreconditionerParams) );
                linearColumnOperator->append(subchannelLinearOperator);
                columnPreconditioner->append(subchannelPreconditioner);
            } else if ( preconditioner=="Banded" ) {
                subchannelPreconditioner_db->putInteger("KL",3);
                subchannelPreconditioner_db->putInteger("KU",3);
                AMP::shared_ptr<AMP::Solver::BandedSolver> subchannelPreconditioner(
                    new AMP::Solver::BandedSolver(subchannelPreconditionerParams) );
                linearColumnOperator->append(subchannelLinearOperator);
                columnPreconditioner->append(subchannelPreconditioner);
            } else if ( preconditioner=="None" ) {
            } else {
                AMP_ERROR("Invalid preconditioner type");
            }
        }

        // create nonlinear solver parameters
        AMP::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(
            new AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db) );
        nonlinearSolverParams->d_comm = globalComm;
        nonlinearSolverParams->d_pOperator = nonlinearCoupledOperator;
        nonlinearSolverParams->d_pInitialGuess = globalSolMultiVector;
        nonlinearSolver.reset(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));

        // create linear solver
        AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver = 
            AMP::dynamic_pointer_cast<AMP::Solver::PetscSNESSolver>(nonlinearSolver)->getKrylovSolver();
        // set preconditioner
        linearSolver->setPreconditioner(columnPreconditioner);
        
    }

    // don't use zero initial guess
    nonlinearSolver->setZeroInitialGuess(false);*/


    // Initialize the pin temperatures
    PROFILE_START("Initialize");
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    int root_subchannel = -1;
    std::vector<double> range(6);
    if ( subchannelMesh.get()!=NULL ) {
        range = subchannelMesh->getBoundingBox();
        AMP_ASSERT(range.size()==6);
        if ( subchannelMesh->getComm().getRank()==0 )
            root_subchannel = globalComm.getRank();
    }
    root_subchannel = globalComm.maxReduce(root_subchannel);
    globalComm.bcast(&range[0],6,root_subchannel);
    // Desired power of the fuel pin (W)
    double P = global_input_db->getDatabase( "SubchannelTwoEqNonlinearOperator" )->getDouble("Rod_Power");
    // Volume of fuel in a 3.81m pin
    double V = 1.939e-4;
    if ( pinMesh.get()!=NULL ) {
        globalThermalSolVec->setToScalar(600);
        AMP::Discretization::DOFManager::shared_ptr gaussPtDOFManager = 
            AMP::Discretization::simpleDOFManager::create(pinMesh,AMP::Mesh::Volume,1,8);
        AMP::Mesh::MeshIterator it = pinMesh->getIterator(AMP::Mesh::Volume,0);
        std::vector<size_t> dofs;
        for (size_t i=0; i<it.size(); i++) {
            gaussPtDOFManager->getDOFs(it->globalID(),dofs);
            for (size_t j=0; j<dofs.size(); j++) {
                specificPowerGpVec->setValueByGlobalID(dofs[j],getPower(range,P,V,it->centroid()));
            }
            ++it;
        }
        if ( cladMesh.get()!=NULL ) {
            AMP::LinearAlgebra::VS_Mesh meshSelector(cladMesh);
            AMP::LinearAlgebra::Vector::shared_ptr cladPower = specificPowerGpVec->select(meshSelector,"cladPower");
            cladPower->zero(); 
        }
        specificPowerGpVec->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
        volumeIntegralColumnOperator->apply(specificPowerGpVec, globalThermalRhsVec);
    }

    if ( subchannelMesh.get()!=NULL ) {
        // get exit pressure
        double Pout = global_input_db->getDatabase( "SubchannelTwoEqNonlinearOperator" )->getDouble("Exit_Pressure");
        // get inlet temperature
        double Tin = global_input_db->getDatabase( "SubchannelTwoEqNonlinearOperator" )->getDouble("Inlet_Temperature");
        // compute inlet enthalpy
        std::map<std::string, AMP::shared_ptr<std::vector<double> > > enthalpyArgMap;
        enthalpyArgMap.insert(std::make_pair("temperature",AMP::shared_ptr<std::vector<double> >(new std::vector<double>(1,Tin))));
        enthalpyArgMap.insert(std::make_pair("pressure",   AMP::shared_ptr<std::vector<double> >(new std::vector<double>(1,Pout))));
        std::vector<double> enthalpyResult(1);
        subchannelPhysicsModel->getProperty("Enthalpy",enthalpyResult,enthalpyArgMap); 
        double hin = enthalpyResult[0];
        std::cout<< "Enthalpy Solution:"<< hin <<std::endl;
        std::cout<< "Outlet pressure:"<< Pout <<std::endl;

        AMP::LinearAlgebra::Vector::shared_ptr subchannelEnthalpy = flowSolVec->select( AMP::LinearAlgebra::VS_Stride(0,2), "H" );
        AMP::LinearAlgebra::Vector::shared_ptr subchannelPressure = flowSolVec->select( AMP::LinearAlgebra::VS_Stride(1,2), "P" );

        subchannelEnthalpy->setToScalar(AMP::Operator::Subchannel::scaleEnthalpy*hin); 
        subchannelPressure->setToScalar(AMP::Operator::Subchannel::scalePressure*Pout); 

        // FIRST APPLY CALL
        AMP::shared_ptr<AMP::Operator::SubchannelOperatorParameters> subchannelLinearParams = 
            AMP::dynamic_pointer_cast<AMP::Operator::SubchannelOperatorParameters>( 
            subchannelNonlinearOperator->getParameters("Jacobian", flowSolVec) );
        subchannelLinearParams->d_initialize = false;
        subchannelLinearOperator->reset(subchannelLinearParams);
        subchannelLinearOperator->residual( flowRhsVec, flowSolVec, flowResVec);
    }

    std::cout << "Reached here " << std::endl;
    nonlinearCoupledOperator->residual(globalRhsMultiVector, globalSolMultiVector, globalResMultiVector);
   
    size_t totalOp;
    if(subchannelMesh != NULL ){
        totalOp = nonlinearColumnOperator->getNumberOfOperators()-1;
    }else{
        totalOp = nonlinearColumnOperator->getNumberOfOperators();
    }
    for( size_t id = 0; id != totalOp; id++) {
        AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearThermalOperator = 
            AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>( nonlinearColumnOperator->getOperator(id));
        nonlinearThermalOperator->modifyInitialSolutionVector(globalThermalSolVec);
        nonlinearThermalOperator->modifyRHSvector(globalThermalRhsVec);
    }
    globalThermalRhsVec->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
    PROFILE_STOP("Initialize");


    // Solve
    PROFILE_START("Solve");
    AMP::pout << "Rhs norm: " << std::setprecision(13)<< globalRhsMultiVector->L2Norm()<< std::endl;
    AMP::pout << "Initial solution norm: " << std::setprecision(13)<< globalSolMultiVector->L2Norm()<< std::endl;
    nonlinearCoupledOperator->residual(globalRhsMultiVector, globalSolMultiVector, globalResMultiVector);
/*    double tempResNorm = 0.0;
    double flowResNorm = 0.0;
    if ( pinMesh!=NULL )
        tempResNorm = globalThermalResVec->L2Norm();
    if ( subchannelMesh!=NULL )
        flowResNorm = flowResVec->L2Norm();
    tempResNorm = globalComm.maxReduce(tempResNorm);
    flowResNorm = globalComm.maxReduce(flowResNorm);
    AMP::pout << "Initial residual norm: " << std::setprecision(13)<< globalResMultiVector->L2Norm()<< std::endl;
    AMP::pout << "Initial temp residual norm: " << std::setprecision(13)<< tempResNorm << std::endl;
    AMP::pout << "Initial flow residual norm: " << std::setprecision(13)<< flowResNorm << std::endl;
//    nonlinearSolver->solve(globalRhsMultiVector, globalSolMultiVector);
    nonlinearCoupledOperator->apply(globalRhsMultiVector, globalSolMultiVector, globalResMultiVector);
    AMP::pout << "Final residual norm: " << std::setprecision(13)<< globalResMultiVector->L2Norm()<< std::endl;
    PROFILE_STOP("Solve");


    // Compute the flow temperature and density
    AMP::LinearAlgebra::Vector::shared_ptr flowTempVec;
    AMP::LinearAlgebra::Vector::shared_ptr deltaFlowTempVec;
    AMP::LinearAlgebra::Vector::shared_ptr flowDensityVec;
    if(subchannelMesh != NULL ){
        flowTempVec = subchannelFuelTemp->cloneVector(); 
        flowDensityVec = subchannelFuelTemp->cloneVector(); 
        int DOFsPerFace[3]={0,0,2};
        AMP::Discretization::DOFManager::shared_ptr faceDOFManager = 
            AMP::Discretization::structuredFaceDOFManager::create( subchannelMesh, DOFsPerFace, 0 );
        DOFsPerFace[2]=1;
        AMP::Discretization::DOFManager::shared_ptr scalarFaceDOFManager = 
            AMP::Discretization::structuredFaceDOFManager::create( subchannelMesh, DOFsPerFace, 0 );
        AMP::Mesh::MeshIterator face  = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
        std::vector<size_t> dofs;
        std::vector<size_t> scalarDofs;
        const double h_scale = 1.0/AMP::Operator::Subchannel::scaleEnthalpy;    // Scale to change the input vector back to correct units
        const double P_scale = 1.0/AMP::Operator::Subchannel::scalePressure;    // Scale to change the input vector back to correct units
        for(size_t i=0; i<face.size(); i++) {
            faceDOFManager->getDOFs( face->globalID(), dofs );
            scalarFaceDOFManager->getDOFs( face->globalID(), scalarDofs );
            std::map<std::string, AMP::shared_ptr<std::vector<double> > > subchannelArgMap;
            AMP::shared_ptr<std::vector<double> > vec1(new std::vector<double>(1,h_scale*flowSolVec->getValueByGlobalID(dofs[0])));
            AMP::shared_ptr<std::vector<double> > vec2(new std::vector<double>(1,P_scale*flowSolVec->getValueByGlobalID(dofs[1])));
            subchannelArgMap.insert(std::make_pair("enthalpy",vec1));
            subchannelArgMap.insert(std::make_pair("pressure",vec2));
            std::vector<double> outTemperatureResult(1);
            subchannelPhysicsModel->getProperty("Temperature", outTemperatureResult, subchannelArgMap); 
            std::vector<double> specificVolume(1);
            subchannelPhysicsModel->getProperty("SpecificVolume", specificVolume, subchannelArgMap);
            flowTempVec->setValueByGlobalID(scalarDofs[0],outTemperatureResult[0]); 
            flowDensityVec->setValueByGlobalID(scalarDofs[0],1.0/specificVolume[0]); 
            ++face;
        } 
        flowTempVec->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
        double Tin = global_input_db->getDatabase( "SubchannelTwoEqNonlinearOperator" )->getDouble("Inlet_Temperature");
        deltaFlowTempVec = flowTempVec->cloneVector();
        deltaFlowTempVec->copyVector(flowTempVec);
        deltaFlowTempVec->addScalar(deltaFlowTempVec,-Tin);
    }
    double flowTempMin = 1e100;
    double flowTempMax = -1e100;
    if ( flowTempVec!=NULL ) {
        flowTempMin = flowTempVec->min();
        flowTempMax = flowTempVec->max();
    }
    flowTempMin = globalComm.minReduce(flowTempMin);
    flowTempMax = globalComm.maxReduce(flowTempMax);
    AMP::pout << "Subchannel Flow Temp Max : " << flowTempMax << " Min : "<< flowTempMin << std::endl;


    // Test the subchannel to point map
    AMP::shared_ptr<AMP::Operator::SubchannelToPointMapParameters> subchannelToPointMapParams(
        new AMP::Operator::SubchannelToPointMapParameters());
    subchannelToPointMapParams->d_Mesh = subchannelMesh;
    subchannelToPointMapParams->d_comm = globalComm;
    subchannelToPointMapParams->d_subchannelPhysicsModel = subchannelPhysicsModel;
    subchannelToPointMapParams->d_outputVar.reset( new AMP::LinearAlgebra::Variable("Density") );
    if(subchannelMesh != NULL ) {
        AMP::Mesh::MeshIterator face  = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
        for(size_t i=0; i<face.size(); i++) {
            std::vector<double> pos = face->centroid();
            subchannelToPointMapParams->x.push_back(pos[0]);
            subchannelToPointMapParams->y.push_back(pos[1]);
            subchannelToPointMapParams->z.push_back(pos[2]);
            ++face;
        }
        AMP_ASSERT(subchannelToPointMapParams->x.size()==flowDensityVec->getLocalSize());
    }
    AMP::Operator::SubchannelToPointMap subchannelDensityToPointMap(subchannelToPointMapParams);
    subchannelToPointMapParams->d_outputVar.reset( new AMP::LinearAlgebra::Variable("Temperature") );
    AMP::Operator::SubchannelToPointMap subchannelTemperatureToPointMap(subchannelToPointMapParams);
    AMP::LinearAlgebra::Vector::shared_ptr densityMapVec = AMP::LinearAlgebra::SimpleVector::create(
        subchannelToPointMapParams->x.size(), subchannelDensityToPointMap.getOutputVariable() );
    AMP::LinearAlgebra::Vector::shared_ptr temperatureMapVec = AMP::LinearAlgebra::SimpleVector::create(
        subchannelToPointMapParams->x.size(), subchannelTemperatureToPointMap.getOutputVariable() );
    subchannelDensityToPointMap.apply( nullVec, flowSolVec, densityMapVec );
    subchannelTemperatureToPointMap.apply( nullVec, flowSolVec, temperatureMapVec );
    if(subchannelMesh != NULL ){
        AMP::Mesh::MeshIterator face  = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
        std::vector<size_t> dofs;
        bool pass_density = true;
        bool pass_temperature = true;
        for(size_t i=0; i<face.size(); i++) {
            flowDensityVec->getDOFManager()->getDOFs( face->globalID(), dofs );
            double density1 = flowDensityVec->getValueByGlobalID(dofs[0]);
            double density2 = densityMapVec->getValueByLocalID(i);
            if ( !AMP::Utilities::approx_equal(density1,density2) )
                pass_density = false;
            double temp1 = flowTempVec->getValueByGlobalID(dofs[0]);
            double temp2 = temperatureMapVec->getValueByLocalID(i);
            if ( !AMP::Utilities::approx_equal(temp1,temp2) )
                pass_temperature = false;
            ++face;
        } 
        if ( pass_density )
            ut->passes("Subchannel density to point map");
        else
            ut->failure("Subchannel density to point map");
        if ( pass_temperature )
            ut->passes("Subchannel temperature to point map");
        else
            ut->failure("Subchannel temperature to point map");
    }
    
    
#ifdef USE_EXT_SILO
    // Rescale the solution to get the correct units
    const double h_scale = 1.0/AMP::Operator::Subchannel::scaleEnthalpy;    // Scale to change the input vector back to correct units
    const double P_scale = 1.0/AMP::Operator::Subchannel::scalePressure;    // Scale to change the input vector back to correct units
    AMP::LinearAlgebra::Vector::shared_ptr enthalpy, pressure;
    enthalpy = flowSolVec->select( AMP::LinearAlgebra::VS_Stride(0,2), "H" );
    pressure = flowSolVec->select( AMP::LinearAlgebra::VS_Stride(1,2), "P" );
    if ( enthalpy.get()!=NULL ) {
        enthalpy->scale(h_scale);
        pressure->scale(P_scale);
    }
    // Register the quantities to plot
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter("Silo");
    if( xyFaceMesh != NULL ){
        siloWriter->registerVector( flowSolVec, xyFaceMesh, AMP::Mesh::Face, "SubchannelFlow" );
        siloWriter->registerVector( flowTempVec, xyFaceMesh, AMP::Mesh::Face, "FlowTemp" );
        siloWriter->registerVector( deltaFlowTempVec, xyFaceMesh, AMP::Mesh::Face, "FlowTempDelta" );
        siloWriter->registerVector( flowDensityVec, xyFaceMesh, AMP::Mesh::Face, "FlowDensity" );
    }
    if ( pinMesh.get()!=NULL ) {
        siloWriter->registerVector( globalThermalSolVec ,  pinMesh , AMP::Mesh::Vertex, "Temperature" );
        siloWriter->registerVector( specificPowerGpVec,  pinMesh , AMP::Mesh::Volume, "Power" );
    }
    siloWriter->writeFile( exeName, 0 );
#endif
    ut->passes("test runs to completion");*/
    globalComm.barrier();
    PROFILE_STOP("Main");
    PROFILE_SAVE("exeName");
}


int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  int num_failed = 0;
  {
    AMP::UnitTest ut;
    PROFILE_ENABLE(0);

    std::string exeName = "testSubchannelSolve-1";
    if(argc == 2) exeName = argv[1];

    SubchannelSolve(&ut,exeName);

    ut.report();
    PROFILE_SAVE(exeName);
    num_failed = ut.NumFailGlobal();
  }
  AMP::AMPManager::shutdown();
  return num_failed;
}   


