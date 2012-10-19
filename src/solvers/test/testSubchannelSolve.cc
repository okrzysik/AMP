#include <string>
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"

#include "boost/shared_ptr.hpp"

#include "utils/InputDatabase.h"
#include "utils/Utilities.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/Database.h"
#include "utils/ProfilerApp.h"

#include "ampmesh/SiloIO.h"
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
#include "operators/subchannel/SubchannelHelpers.h"
#include "operators/subchannel/CoupledChannelToCladMapOperator.h"
#include "operators/subchannel/SubchannelConstants.h"
#include "operators/subchannel/SubchannelDensityToPointMap.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"

#include "solvers/PetscKrylovSolver.h"
#include "solvers/PetscSNESSolver.h"
#include "solvers/TrilinosMLSolver.h"
#include "solvers/ColumnSolver.h"

#include "ampmesh/StructuredMeshHelper.h"


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
        int DofsPerFace =  2;
        AMP::Discretization::DOFManager::shared_ptr faceDOFManager = 
            AMP::Discretization::simpleDOFManager::create( subchannelMesh, 
            AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,1), 
            AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,0), DofsPerFace );
        // dof manager for the scalar quanties on the z faces - for mapped temperature clad temp onto subchannel discretization
        /*
        AMP::Discretization::DOFManager::shared_ptr scalarFaceDOFManager = 
            AMP::Discretization::simpleDOFManager::create( subchannelMesh, 
            AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,1), 
            AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,0), 1);
        */
        // create solution, rhs, and residual vectors
        flowVec = AMP::LinearAlgebra::createVector( faceDOFManager , flowVariable , true );
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
        xyFaceMesh = subchannelMesh->Subset( AMP::Mesh::StructuredMeshHelper::getXYFaceIterator( subchannelMesh , 0 ) );
    }

    // Variables
    AMP::LinearAlgebra::Variable::shared_ptr thermalVariable (new AMP::LinearAlgebra::Variable("Temperature"));    // temperature on pellets and cladding
    AMP::LinearAlgebra::Variable::shared_ptr flowVariable (new AMP::LinearAlgebra::Variable("Flow"));              // enthalpy and pressure on z-faces
    AMP::LinearAlgebra::Variable::shared_ptr powerVariable (new AMP::LinearAlgebra::Variable("SpecificPowerInWattsPerGram")); // specific power on gauss points

    boost::shared_ptr<AMP::Operator::Operator> thermalCopyOperator;

    boost::shared_ptr<AMP::Operator::OperatorParameters> emptyParams;
    boost::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator (new AMP::Operator::ColumnOperator(emptyParams));
    boost::shared_ptr<AMP::Operator::ColumnOperator> linearColumnOperator (new AMP::Operator::ColumnOperator(emptyParams));
    boost::shared_ptr<AMP::Operator::ColumnOperator> volumeIntegralColumnOperator(new AMP::Operator::ColumnOperator(emptyParams));

    boost::shared_ptr<AMP::Operator::ColumnOperator> mapsColumn( new AMP::Operator::ColumnOperator ( emptyParams) );
    boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator> n2nColumn( new AMP::Operator::AsyncMapColumnOperator ( emptyParams) );
    boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator> szaColumn( new AMP::Operator::AsyncMapColumnOperator ( emptyParams) );

    boost::shared_ptr<AMP::Solver::PetscSNESSolver>  nonlinearCoupledSolver;
    boost::shared_ptr<AMP::Solver::PetscKrylovSolver>  linearColumnSolver;
    boost::shared_ptr<AMP::Solver::ColumnSolver>  columnPreconditioner;

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
            boost::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;
            boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> thermalNonlinearOperator = boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator( adapter,
                prefix+"NonlinearThermalOperator",
                global_input_db,
                thermalTransportModel));
            nonlinearColumnOperator->append(thermalNonlinearOperator);

            // CREATE THE LINEAR THERMAL OPERATOR 1 
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

    }

    AMP::LinearAlgebra::Vector::shared_ptr subchannelFuelTemp; 
    AMP::LinearAlgebra::Vector::shared_ptr subchannelFlowTemp;
    if ( subchannelMesh.get()!=NULL ) {
        AMP::Discretization::DOFManager::shared_ptr scalarFaceDOFManager = 
            AMP::Discretization::simpleDOFManager::create( subchannelMesh, 
            AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,1), 
            AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,0), 1);
        subchannelFuelTemp = AMP::LinearAlgebra::createVector( scalarFaceDOFManager , thermalVariable );
        subchannelFlowTemp = AMP::LinearAlgebra::createVector( scalarFaceDOFManager , thermalVariable );
    }

    // get subchannel physics model
    // for post processing - need water library to convert enthalpy to temperature...
    boost::shared_ptr<AMP::Database> subchannelPhysics_db = global_input_db->getDatabase("SubchannelPhysicsModel");
    boost::shared_ptr<AMP::Operator::ElementPhysicsModelParameters> params( new AMP::Operator::ElementPhysicsModelParameters(subchannelPhysics_db));
    boost::shared_ptr<AMP::Operator::SubchannelPhysicsModel>  subchannelPhysicsModel (new AMP::Operator::SubchannelPhysicsModel(params));
    boost::shared_ptr<AMP::Operator::SubchannelTwoEqNonlinearOperator> subchannelNonlinearOperator;
    boost::shared_ptr<AMP::Operator::SubchannelTwoEqLinearOperator> subchannelLinearOperator; 

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
                subchannelNonlinearOperator = boost::dynamic_pointer_cast<AMP::Operator::SubchannelTwoEqNonlinearOperator>(
                    AMP::Operator::OperatorBuilder::createOperator( adapter ,"SubchannelTwoEqNonlinearOperator",global_input_db ) );
                boost::shared_ptr<AMP::Operator::SubchannelOperatorParameters> nonlinearOpParams = subchannelNonlinearOperator->getParams();
                nonlinearOpParams->clad_x = clad_x;
                nonlinearOpParams->clad_y = clad_y;
                nonlinearOpParams->clad_d = clad_d;
                subchannelNonlinearOperator->reset(nonlinearOpParams);
                // create linear operator
                subchannelLinearOperator = boost::dynamic_pointer_cast<AMP::Operator::SubchannelTwoEqLinearOperator>(
                    AMP::Operator::OperatorBuilder::createOperator( adapter ,"SubchannelTwoEqLinearOperator",
                    global_input_db, subchannelNonlinearOperator->getSubchannelPhysicsModel() ) );
                subchannelNonlinearOperator->setVector(subchannelFuelTemp); 
                subchannelLinearOperator->reset(subchannelNonlinearOperator->getJacobianParameters(subchannelFuelTemp));
                // pass creation test
                ut->passes(exeName+": creation");
                std::cout.flush();
                nonlinearColumnOperator->append(subchannelNonlinearOperator);
                linearColumnOperator->append(subchannelLinearOperator);
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

        std::vector<AMP::Mesh::Mesh::shared_ptr> pins = boost::dynamic_pointer_cast<AMP::Mesh::MultiMesh>(pinMesh)->getMeshes();

        for (size_t i=0; i<pins.size(); i++) {
            if ( global_input_db->keyExists("ThermalNodeToNodeMaps") ) {
                boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator> map = 
                    AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::NodeToNodeMap>( 
                    pins[i], global_input_db->getDatabase("ThermalNodeToNodeMaps") );
                for (size_t j=0; j<map->getNumberOfOperators(); j++)
                    n2nColumn->append( map->getOperator(j) );
            }
            if ( global_input_db->keyExists("ThermalScalarZAxisMaps") ) {
                boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator> sza = 
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

            boost::shared_ptr<AMP::Operator::NonlinearBVPOperator>  curBVPop = 
                boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>( nonlinearColumnOperator->getOperator ( curOperator ) );
            boost::shared_ptr<AMP::Operator::ColumnBoundaryOperator>  curBCcol = 
                boost::dynamic_pointer_cast<AMP::Operator::ColumnBoundaryOperator> ( curBVPop->getBoundaryOperator () );
            boost::shared_ptr<AMP::Database> operator_db = global_input_db->getDatabase ( prefix+"NonlinearThermalOperator" );
            boost::shared_ptr<AMP::Database>  curBCdb = global_input_db->getDatabase(operator_db->getString ( "BoundaryOperator" ));
            std::vector<std::string>  opNames = curBCdb->getStringArray( "boundaryOperators" );
            int  opNumber = curBCdb->getInteger( "numberOfBoundaryOperators" );
            for (int curBCentry=0; curBCentry!=opNumber; curBCentry++ ) {
                if ( opNames[curBCentry] == "P2CRobinVectorCorrection" ) {
                    boost::shared_ptr<AMP::Operator::RobinVectorCorrection>  gapBC = 
                        boost::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection> ( curBCcol->getBoundaryOperator ( curBCentry ) );
                    AMP_ASSERT(thermalMapVec!=NULL);
                    gapBC->setVariableFlux ( thermalMapVec );
                    gapBC->reset ( gapBC->getParameters() );
                } else if ( ( opNames[curBCentry]=="BottomP2PNonlinearRobinVectorCorrection"  ) || 
                            ( opNames[curBCentry]=="MiddleP2PNonlinearRobinBoundaryCondition" ) || 
                            ( opNames[curBCentry]=="TopP2PNonlinearRobinBoundaryCondition"    ) ) {
                    boost::shared_ptr<AMP::Operator::RobinVectorCorrection>  p2pBC = 
                        boost::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection> ( curBCcol->getBoundaryOperator ( curBCentry ) );
                    AMP_ASSERT(thermalMapVec!=NULL);
                    p2pBC->setVariableFlux ( thermalMapVec );
                    p2pBC->reset ( p2pBC->getParameters() );
                } else if ( opNames[curBCentry] == "C2WBoundaryVectorCorrection" ) {
                    boost::shared_ptr<AMP::Database>  thisDb = global_input_db->getDatabase( opNames[curBCentry] );
                    bool isCoupled = thisDb->getBoolWithDefault( "IsCoupledBoundary_0", false);
                    if( isCoupled ) {
                        boost::shared_ptr<AMP::Operator::RobinVectorCorrection>  c2wBC = 
                            boost::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection>( curBCcol->getBoundaryOperator ( curBCentry ) );
                        AMP_ASSERT(thermalMapVec!=NULL);
                        c2wBC->setVariableFlux ( thermalMapVec );
                        c2wBC->setFrozenVector( density_map_vec );
                        c2wBC->reset ( c2wBC->getParameters() );
                    }
                } else if ( opNames[curBCentry] == "C2PRobinVectorCorrection" ) {
                    boost::shared_ptr<AMP::Operator::RobinVectorCorrection>  gapBC = 
                        boost::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection>( curBCcol->getBoundaryOperator ( curBCentry ) );
                    AMP_ASSERT(thermalMapVec!=NULL);
                    gapBC->setVariableFlux ( thermalMapVec );
                    gapBC->reset ( gapBC->getParameters() );
                } else {
                    AMP_ERROR("Unknown boundary operator");
                }
            }
            curOperator++;
        }
    }

    // Create the maps from the clad temperature to the subchannel temperature
    boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator>  cladToSubchannelMap; 
    boost::shared_ptr<AMP::Database> cladToSubchannelDb = global_input_db->getDatabase( "CladToSubchannelMaps" );
    cladToSubchannelMap = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::CladToSubchannelMap>( manager, cladToSubchannelDb );
    cladToSubchannelMap->setVector( subchannelFuelTemp );
    mapsColumn->append( cladToSubchannelMap );

    // Create the maps from the flow variable (enthalpy and pressure) on subchannel mesh and convert to temperature and density then map to clad surface
    boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator>  thermalSubchannelToCladMap, densitySubchannelToCladMap; 
    boost::shared_ptr<AMP::Database> thermalCladToSubchannelDb = global_input_db->getDatabase( "ThermalSubchannelToCladMaps" );
    thermalSubchannelToCladMap = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::SubchannelToCladMap>( manager, thermalCladToSubchannelDb );
    boost::shared_ptr<AMP::Database> densityCladToSubchannelDb = global_input_db->getDatabase( "DensitySubchannelToCladMaps" );
    densitySubchannelToCladMap = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::SubchannelToCladMap>( manager, densityCladToSubchannelDb );
    if ( cladMesh.get()!=NULL ) {
        AMP::LinearAlgebra::VS_Comm commSelector(cladMesh->getComm());
        AMP::LinearAlgebra::Vector::shared_ptr subsetTheramlVec = thermalMapVec->select(commSelector, thermalMapVec->getVariable()->getName());
        thermalSubchannelToCladMap->setVector( subsetTheramlVec );
        densitySubchannelToCladMap->setVector( density_map_vec );
    }
    boost::shared_ptr<AMP::InputDatabase> emptyDb (new AMP::InputDatabase("empty"));
    emptyDb->putInteger("print_info_level",0); 
    boost::shared_ptr<AMP::Operator::CoupledChannelToCladMapOperatorParameters> coupledChannelMapOperatorParams(new AMP::Operator::CoupledChannelToCladMapOperatorParameters( emptyDb ));
    coupledChannelMapOperatorParams->d_variable       = flowVariable;
    coupledChannelMapOperatorParams->d_vector         = subchannelFlowTemp;
    coupledChannelMapOperatorParams->d_Mesh           = subchannelMesh;
    coupledChannelMapOperatorParams->d_thermalMapOperator = thermalSubchannelToCladMap;
    coupledChannelMapOperatorParams->d_densityMapOperator = densitySubchannelToCladMap;
    coupledChannelMapOperatorParams->d_subchannelMesh = subchannelMesh;
    coupledChannelMapOperatorParams->d_subchannelPhysicsModel = subchannelPhysicsModel ;
    boost::shared_ptr<AMP::Operator::Operator> coupledChannelMapOperator (new AMP::Operator::CoupledChannelToCladMapOperator(coupledChannelMapOperatorParams));
    mapsColumn->append( coupledChannelMapOperator );

    if ( pinMesh.get()!=NULL ) {
        boost::shared_ptr<AMP::InputDatabase> copyOp_db = boost::dynamic_pointer_cast<AMP::InputDatabase>(global_input_db->getDatabase("CopyOperator"));
        boost::shared_ptr<AMP::Operator::VectorCopyOperatorParameters> vecCopyOperatorParams(new AMP::Operator::VectorCopyOperatorParameters( copyOp_db ));
        vecCopyOperatorParams->d_copyVariable = thermalVariable;
        vecCopyOperatorParams->d_copyVector = thermalMapVec;
        vecCopyOperatorParams->d_Mesh = pinMesh ;
        thermalCopyOperator.reset(new AMP::Operator::VectorCopyOperator(vecCopyOperatorParams));
        thermalMapVec->zero();
    }

    boost::shared_ptr<AMP::Operator::CoupledOperatorParameters> CoupledOpParams(new AMP::Operator::CoupledOperatorParameters(emptyDb));
    CoupledOpParams->d_CopyOperator = thermalCopyOperator;
    CoupledOpParams->d_MapOperator = mapsColumn;
    CoupledOpParams->d_BVPOperator = nonlinearColumnOperator;
    boost::shared_ptr<AMP::Operator::Operator> nonlinearCoupledOperator(new AMP::Operator::CoupledOperator(CoupledOpParams));


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
    for(unsigned int id = 0; id != linearColumnOperator->getNumberOfOperators(); id++) {
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
        volumeIntegralColumnOperator->apply(nullVec, specificPowerGpVec, globalThermalRhsVec , 1., 0.);
    }

    if ( subchannelMesh.get()!=NULL ) {
        // get exit pressure
        double Pout = global_input_db->getDatabase( "SubchannelTwoEqNonlinearOperator" )->getDouble("Exit_Pressure");
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
        std::cout<< "Outlet pressure:"<< Pout <<std::endl;

        AMP::LinearAlgebra::Vector::shared_ptr subchannelEnthalpy = flowSolVec->select( AMP::LinearAlgebra::VS_Stride(0,2), "H" );
        AMP::LinearAlgebra::Vector::shared_ptr subchannelPressure = flowSolVec->select( AMP::LinearAlgebra::VS_Stride(1,2), "P" );

        subchannelEnthalpy->setToScalar(AMP::Operator::Subchannel::scaleEnthalpy*hin); 
        subchannelPressure->setToScalar(AMP::Operator::Subchannel::scalePressure*Pout); 

        // FIRST APPLY CALL
        subchannelLinearOperator->reset(subchannelNonlinearOperator->getJacobianParameters(flowSolVec));
        subchannelLinearOperator->apply( flowRhsVec, flowSolVec, flowResVec, 1.0, -1.0);
    }

    nonlinearCoupledOperator->apply(globalRhsMultiVector, globalSolMultiVector, globalResMultiVector, 1.0, -1.0);
   
    size_t totalOp;
    if(subchannelMesh != NULL ){
        totalOp = nonlinearColumnOperator->getNumberOfOperators()-1;
    }else{
        totalOp = nonlinearColumnOperator->getNumberOfOperators();
    }
    for( size_t id = 0; id != totalOp; id++) {
        boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearThermalOperator = 
            boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>( nonlinearColumnOperator->getOperator(id));
        nonlinearThermalOperator->modifyInitialSolutionVector(globalThermalSolVec);
        nonlinearThermalOperator->modifyRHSvector(globalThermalRhsVec);
    }
    globalThermalRhsVec->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
    PROFILE_STOP("Initialize");


    // Solve
    PROFILE_START("Solve");
    AMP::pout << "Rhs norm: " << std::setprecision(13)<< globalRhsMultiVector->L2Norm()<< std::endl;
    AMP::pout << "Initial solution norm: " << std::setprecision(13)<< globalSolMultiVector->L2Norm()<< std::endl;
    nonlinearCoupledOperator->apply(globalRhsMultiVector, globalSolMultiVector, globalResMultiVector);
    double tempResNorm = 0.0;
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
    nonlinearSolver->solve(globalRhsMultiVector, globalSolMultiVector);
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
        int DofsPerFace =  2;
        AMP::Discretization::DOFManager::shared_ptr faceDOFManager = 
            AMP::Discretization::simpleDOFManager::create( subchannelMesh, 
            AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,1), 
            AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,0), DofsPerFace );
        AMP::Discretization::DOFManager::shared_ptr scalarFaceDOFManager = 
            AMP::Discretization::simpleDOFManager::create( subchannelMesh, 
            AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,1), 
            AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,0), 1);
        AMP::Mesh::MeshIterator face  = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
        AMP::Mesh::MeshIterator end_face = face.end();
        std::vector<size_t> dofs;
        std::vector<size_t> scalarDofs;
        const double h_scale = 1.0/AMP::Operator::Subchannel::scaleEnthalpy;    // Scale to change the input vector back to correct units
        const double P_scale = 1.0/AMP::Operator::Subchannel::scalePressure;    // Scale to change the input vector back to correct units
        for(size_t i=0; i<face.size(); i++) {
            faceDOFManager->getDOFs( face->globalID(), dofs );
            scalarFaceDOFManager->getDOFs( face->globalID(), scalarDofs );
            std::map<std::string, boost::shared_ptr<std::vector<double> > > subchannelArgMap;
            subchannelArgMap.insert(std::make_pair("enthalpy",new std::vector<double>(1,h_scale*flowSolVec->getValueByGlobalID(dofs[0]))));
            subchannelArgMap.insert(std::make_pair("pressure",new std::vector<double>(1,P_scale*flowSolVec->getValueByGlobalID(dofs[1]))));
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


    // Test the density to point map
    boost::shared_ptr<AMP::Operator::SubchannelDensityToPointMapParameters> subchannelToPointMapParams(
        new AMP::Operator::SubchannelDensityToPointMapParameters());
    subchannelToPointMapParams->d_Mesh = subchannelMesh;
    subchannelToPointMapParams->d_comm = globalComm;
    subchannelToPointMapParams->d_subchannelPhysicsModel = subchannelPhysicsModel;
    if(subchannelMesh != NULL ){
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
    AMP::Operator::SubchannelDensityToPointMap subchannelDensityToPointMap(subchannelToPointMapParams);
    AMP::LinearAlgebra::Vector::shared_ptr densityMapVec = AMP::LinearAlgebra::SimpleVector::create(
        subchannelToPointMapParams->x.size(), subchannelDensityToPointMap.getOutputVariable() );
    subchannelDensityToPointMap.apply( nullVec, flowSolVec, densityMapVec );
    if(subchannelMesh != NULL ){
        AMP::Mesh::MeshIterator face  = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
        std::vector<size_t> dofs;
        bool pass = true;
        for(size_t i=0; i<face.size(); i++) {
            flowDensityVec->getDOFManager()->getDOFs( face->globalID(), dofs );
            double density1 = flowDensityVec->getValueByGlobalID(dofs[0]);
            double density2 = densityMapVec->getValueByLocalID(i);
            if ( !AMP::Utilities::approx_equal(density1,density2) )
                pass = false;
            ++face;
        } 
        if ( pass )
            ut->passes("Subchannel density to point map");
        else
            ut->failure("Subchannel density to point map");
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
    AMP::Mesh::SiloIO::shared_ptr  siloWriter( new AMP::Mesh::SiloIO );
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
    siloWriter->writeFile( silo_name , 0 );
#endif
    ut->passes("test runs to completion");
    PROFILE_STOP("Main");
    PROFILE_SAVE("exeName");
}


int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;
  PROFILE_ENABLE(0);

  std::string exeName = "testSubchannelSolve-1";
  if(argc == 2) exeName = argv[1];

  SubchannelSolve(&ut,exeName);

  ut.report();
  PROFILE_SAVE(exeName);

  int num_failed = ut.NumFailGlobal();
  AMP::AMPManager::shutdown();
  return num_failed;
}   


