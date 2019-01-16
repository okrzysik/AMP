#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MultiMesh.h"
#include "AMP/ampmesh/StructuredMeshHelper.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/discretization/structuredFaceDOFManager.h"
#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/CoupledOperator.h"
#include "AMP/operators/CoupledOperatorParameters.h"
#include "AMP/operators/ElementOperationFactory.h"
#include "AMP/operators/ElementPhysicsModelFactory.h"
#include "AMP/operators/IdentityOperator.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/VectorCopyOperator.h"
#include "AMP/operators/boundary/ColumnBoundaryOperator.h"
#include "AMP/operators/boundary/DirichletMatrixCorrection.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/boundary/libmesh/NeumannVectorCorrection.h"
#include "AMP/operators/boundary/libmesh/NeumannVectorCorrectionParameters.h"
#include "AMP/operators/boundary/libmesh/RobinVectorCorrection.h"
#include "AMP/operators/libmesh/VolumeIntegralOperator.h"
#include "AMP/operators/map/AsyncMapColumnOperator.h"
#include "AMP/operators/map/CladToSubchannelMap.h"
#include "AMP/operators/map/GaussPointToGaussPointMap.h"
#include "AMP/operators/map/Map1Dto3D.h"
#include "AMP/operators/map/Map3Dto1D.h"
#include "AMP/operators/map/MapSurface.h"
#include "AMP/operators/map/NodeToNodeMap.h"
#include "AMP/operators/map/ScalarN2GZAxisMap.h"
#include "AMP/operators/map/ScalarZAxisMap.h"
#include "AMP/operators/map/SubchannelToCladMap.h"
#include "AMP/operators/subchannel/CoupledChannelToCladMapOperator.h"
#include "AMP/operators/subchannel/SubchannelConstants.h"
#include "AMP/operators/subchannel/SubchannelHelpers.h"
#include "AMP/operators/subchannel/SubchannelPhysicsModel.h"
#include "AMP/operators/subchannel/SubchannelToPointMap.h"
#include "AMP/operators/subchannel/SubchannelTwoEqLinearOperator.h"
#include "AMP/operators/subchannel/SubchannelTwoEqNonlinearOperator.h"
#include "AMP/solvers/BandedSolver.h"
#include "AMP/solvers/ColumnSolver.h"
#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#include "AMP/solvers/trilinos/nox/TrilinosNOXSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/InputDatabase.h"
#include "AMP/utils/InputManager.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/utils/shared_ptr.h"
#include "AMP/vectors/SimpleVector.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/VectorSelector.h"

#include "ProfilerApp.h"

#include <string>


// Function to get an arbitrary power profile (W/kg) assuming a density of 1 kg/m^3 for the volume
// integral
// P is the total power of the pin, V is the volume of the pin
static double
getPower( const std::vector<double> &range, double P, double V, const AMP::Mesh::Point &pos )
{
    const double pi = 3.1415926535897932;
    double x        = ( pos[0] - range[0] ) / ( range[1] - range[0] );
    double y        = ( pos[1] - range[2] ) / ( range[3] - range[2] );
    double z        = ( pos[2] - range[4] ) / ( range[5] - range[4] );
    return P / V * ( 0.8 + 0.2 * x + 0.2 * y ) * pi / 2 * sin( pi * z );
}


// Function to create the solution vectors
static void createVectors( AMP::Mesh::Mesh::shared_ptr pinMesh,
                           AMP::Mesh::Mesh::shared_ptr subchannelMesh,
                           AMP::LinearAlgebra::Vector::shared_ptr &globalMultiVector,
                           AMP::LinearAlgebra::Vector::shared_ptr &specificPowerGpVec )
{
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    auto multivec     = AMP::LinearAlgebra::MultiVector::create( "multivector", globalComm );
    globalMultiVector = multivec;

    auto thermalVariable = AMP::make_shared<AMP::LinearAlgebra::Variable>( "Temperature" );
    auto flowVariable    = AMP::make_shared<AMP::LinearAlgebra::Variable>( "Flow" );
    auto powerVariable =
        AMP::make_shared<AMP::LinearAlgebra::Variable>( "SpecificPowerInWattsPerGram" );

    AMP::LinearAlgebra::Vector::shared_ptr thermalVec;
    if ( pinMesh.get() != nullptr ) {
        auto nodalScalarDOF = AMP::Discretization::simpleDOFManager::create(
            pinMesh, AMP::Mesh::GeomType::Vertex, 1, 1, true );
        thermalVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, thermalVariable );
    }
    multivec->addVector( thermalVec );

    AMP::LinearAlgebra::Vector::shared_ptr flowVec;
    if ( subchannelMesh.get() != nullptr ) {
        int DOFsPerFace[3] = { 0, 0, 2 };
        auto faceDOFManager =
            AMP::Discretization::structuredFaceDOFManager::create( subchannelMesh, DOFsPerFace, 0 );
        // create solution, rhs, and residual vectors
        flowVec = AMP::LinearAlgebra::createVector( faceDOFManager, flowVariable, true );
    }
    multivec->addVector( flowVec );

    if ( pinMesh.get() != nullptr ) {
        auto gaussPtDOFManager = AMP::Discretization::simpleDOFManager::create(
            pinMesh, AMP::Mesh::GeomType::Volume, 1, 8 );
        specificPowerGpVec = AMP::LinearAlgebra::createVector( gaussPtDOFManager, powerVariable );
        specificPowerGpVec->setToScalar( 0.0 );
    }
}


static void SubchannelSolve( AMP::UnitTest *ut, const std::string &exeName )
{
    PROFILE_START( "Main" );
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;
    AMP::PIO::logAllNodes( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    globalComm.barrier();

    // Read the input file
    auto global_input_db = AMP::make_shared<AMP::InputDatabase>( "input_db" );
    AMP::InputManager::getManager()->parseInputFile( input_file, global_input_db );
    global_input_db->printClassData( AMP::plog );

    // Get the Mesh database and create the mesh parameters
    auto database   = global_input_db->getDatabase( "Mesh" );
    auto meshParams = AMP::make_shared<AMP::Mesh::MeshParameters>( database );
    meshParams->setComm( globalComm );

    // Get the meshes
    auto manager = AMP::Mesh::Mesh::buildMesh( meshParams );
    auto pinMesh = manager->Subset( "MultiPin" );
    AMP::Mesh::Mesh::shared_ptr cladMesh;
    if ( pinMesh.get() != nullptr ) {
        pinMesh->setName( "MultiPin" );
        cladMesh = pinMesh->Subset( "clad" );
    }
    auto subchannelMesh = manager->Subset( "subchannel" );
    AMP::Mesh::Mesh::shared_ptr xyFaceMesh;
    if ( subchannelMesh.get() != nullptr ) {
        auto face  = AMP::Mesh::StructuredMeshHelper::getXYFaceIterator( subchannelMesh, 0 );
        xyFaceMesh = subchannelMesh->Subset( face );
    }

    // Variables
    auto thermalVariable = AMP::make_shared<AMP::LinearAlgebra::Variable>( "Temperature" );
    auto flowVariable    = AMP::make_shared<AMP::LinearAlgebra::Variable>( "Flow" );
    auto powerVariable =
        AMP::make_shared<AMP::LinearAlgebra::Variable>( "SpecificPowerInWattsPerGram" );

    AMP::shared_ptr<AMP::Operator::Operator> thermalCopyOperator;

    AMP::shared_ptr<AMP::Operator::OperatorParameters> emptyParams;
    auto nonlinearColumnOperator = AMP::make_shared<AMP::Operator::ColumnOperator>( emptyParams );
    auto linearColumnOperator    = AMP::make_shared<AMP::Operator::ColumnOperator>( emptyParams );
    auto volumeIntegralColumnOperator =
        AMP::make_shared<AMP::Operator::ColumnOperator>( emptyParams );

    auto mapsColumn = AMP::make_shared<AMP::Operator::ColumnOperator>( emptyParams );
    auto n2nColumn  = AMP::make_shared<AMP::Operator::AsyncMapColumnOperator>( emptyParams );
    auto szaColumn  = AMP::make_shared<AMP::Operator::AsyncMapColumnOperator>( emptyParams );

    AMP::shared_ptr<AMP::Solver::TrilinosNOXSolver> nonlinearCoupledSolver;
    AMP::shared_ptr<AMP::Solver::ColumnSolver> columnPreconditioner;

    if ( pinMesh.get() != nullptr ) {
        auto pinMeshIDs = pinMesh->getBaseMeshIDs();

        // CREATE OPERATORS
        for ( auto pinMeshID : pinMeshIDs ) {
            auto adapter = manager->Subset( pinMeshID );
            if ( adapter.get() == nullptr )
                continue;

            std::string meshName = adapter->getName();
            std::string prefix, prefixPower;

            if ( meshName.compare( "clad" ) == 0 ) {
                prefix      = "Clad";
                prefixPower = "Clad";
            } else if ( meshName.compare( "pellet_1" ) == 0 ) {
                prefix      = "BottomPellet";
                prefixPower = "Pellet";
            } else if ( meshName.compare( "pellet_3" ) == 0 ) {
                prefix      = "TopPellet";
                prefixPower = "Pellet";
            } else if ( meshName.compare( 0, 7, "pellet_" ) == 0 ) {
                prefix      = "MiddlePellet";
                prefixPower = "Pellet";
            } else {
                AMP_ERROR( "Unknown Mesh" );
            }

            // CREATE THE NONLINEAR THERMAL OPERATOR 1
            AMP_INSIST( global_input_db->keyExists( prefix + "NonlinearThermalOperator" ),
                        "key missing!" );
            AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;
            auto thermalNonlinearOperator =
                AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
                    AMP::Operator::OperatorBuilder::createOperator( adapter,
                                                                    prefix +
                                                                        "NonlinearThermalOperator",
                                                                    global_input_db,
                                                                    thermalTransportModel ) );
            nonlinearColumnOperator->append( thermalNonlinearOperator );

            // CREATE THE LINEAR THERMAL OPERATOR 1
            AMP_INSIST( global_input_db->keyExists( prefix + "LinearThermalOperator" ),
                        "key missing!" );
            auto thermalLinearOperator =
                AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
                    AMP::Operator::OperatorBuilder::createOperator( adapter,
                                                                    prefix +
                                                                        "LinearThermalOperator",
                                                                    global_input_db,
                                                                    thermalTransportModel ) );
            linearColumnOperator->append( thermalLinearOperator );

            AMP_INSIST( global_input_db->keyExists( prefixPower + "VolumeIntegralOperator" ),
                        "key missing!" );
            AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
            auto specificPowerGpVecToPowerDensityNodalVecOperator =
                AMP::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
                    AMP::Operator::OperatorBuilder::createOperator( adapter,
                                                                    prefixPower +
                                                                        "VolumeIntegralOperator",
                                                                    global_input_db,
                                                                    stransportModel ) );
            volumeIntegralColumnOperator->append(
                specificPowerGpVecToPowerDensityNodalVecOperator );
        }
    }

    // Get the subchannel hydraulic diameter for the temperature boundary condition
    auto ChannelDiameterVec =
        AMP::Operator::Subchannel::getCladHydraulicDiameter( cladMesh, subchannelMesh, globalComm );


    AMP::LinearAlgebra::Vector::shared_ptr subchannelFuelTemp;
    AMP::LinearAlgebra::Vector::shared_ptr subchannelFlowTemp;
    if ( subchannelMesh.get() != nullptr ) {
        int DOFsPerFace[3] = { 0, 0, 1 };
        auto scalarFaceDOFManager =
            AMP::Discretization::structuredFaceDOFManager::create( subchannelMesh, DOFsPerFace, 0 );
        subchannelFuelTemp =
            AMP::LinearAlgebra::createVector( scalarFaceDOFManager, thermalVariable );
        subchannelFlowTemp =
            AMP::LinearAlgebra::createVector( scalarFaceDOFManager, thermalVariable );
    }

    // get subchannel physics model
    // for post processing - need water library to convert enthalpy to temperature...
    auto subchannelPhysics_db = global_input_db->getDatabase( "SubchannelPhysicsModel" );
    auto params =
        AMP::make_shared<AMP::Operator::ElementPhysicsModelParameters>( subchannelPhysics_db );
    auto subchannelPhysicsModel = AMP::make_shared<AMP::Operator::SubchannelPhysicsModel>( params );
    AMP::shared_ptr<AMP::Operator::SubchannelTwoEqNonlinearOperator> subchannelNonlinearOperator;
    AMP::shared_ptr<AMP::Operator::LinearOperator> subchannelLinearOperator;

    // Get the subchannel operators
    std::vector<double> clad_x, clad_y, clad_d;
    AMP::Operator::Subchannel::getCladProperties( globalComm, cladMesh, clad_x, clad_y, clad_d );
    if ( subchannelMesh.get() != nullptr ) {
        auto subChannelMeshIDs = subchannelMesh->getBaseMeshIDs();

        for ( auto subChannelMeshID : subChannelMeshIDs ) {
            auto adapter = manager->Subset( subChannelMeshID );
            if ( adapter.get() == nullptr )
                continue;

            auto meshName = adapter->getName();
            if ( meshName.compare( "subchannel" ) == 0 ) {
                // create the non-linear operator
                subchannelNonlinearOperator =
                    AMP::dynamic_pointer_cast<AMP::Operator::SubchannelTwoEqNonlinearOperator>(
                        AMP::Operator::OperatorBuilder::createOperator(
                            adapter, "SubchannelTwoEqNonlinearOperator", global_input_db ) );
                auto nonlinearOpParams    = subchannelNonlinearOperator->getParams();
                nonlinearOpParams->clad_x = clad_x;
                nonlinearOpParams->clad_y = clad_y;
                nonlinearOpParams->clad_d = clad_d;
                subchannelNonlinearOperator->reset( nonlinearOpParams );
                // create linear operator
                subchannelLinearOperator =
                    AMP::dynamic_pointer_cast<AMP::Operator::SubchannelTwoEqLinearOperator>(
                        AMP::Operator::OperatorBuilder::createOperator(
                            adapter,
                            "SubchannelTwoEqLinearOperator",
                            global_input_db,
                            subchannelNonlinearOperator->getSubchannelPhysicsModel() ) );
                // subchannelLinearOperator.reset( new AMP::Operator::IdentityOperator(
                // nonlinearOpParams ) );
                int DOFsPerFace[3]  = { 0, 0, 2 };
                auto flowDOFManager = AMP::Discretization::structuredFaceDOFManager::create(
                    subchannelMesh, DOFsPerFace, 0 );
                auto subchannelFlow =
                    AMP::LinearAlgebra::createVector( flowDOFManager, flowVariable );
                subchannelNonlinearOperator->setVector( subchannelFuelTemp );
                auto subchannelLinearParams =
                    AMP::dynamic_pointer_cast<AMP::Operator::SubchannelOperatorParameters>(
                        subchannelNonlinearOperator->getParameters( "Jacobian", subchannelFlow ) );
                subchannelLinearParams->d_initialize = false;
                subchannelLinearOperator->reset( subchannelLinearParams );
                // pass creation test
                ut->passes( exeName + ": creation" );
                std::cout.flush();
                nonlinearColumnOperator->append( subchannelNonlinearOperator );
                // Do not add the subchannel to the linear operator (we will add it later)
            }
        }
    }


    // CREATE MAPS
    AMP::LinearAlgebra::Vector::shared_ptr thermalMapVec;
    AMP::LinearAlgebra::Vector::shared_ptr density_map_vec;
    if ( cladMesh.get() != nullptr ) {
        auto nodalScalarDOF = AMP::Discretization::simpleDOFManager::create(
            cladMesh, AMP::Mesh::GeomType::Vertex, 1, 1, true );
        auto densityVariable = AMP::make_shared<AMP::LinearAlgebra::Variable>( "Density" );
        density_map_vec      = AMP::LinearAlgebra::createVector( nodalScalarDOF, densityVariable );
        density_map_vec->zero();
    }
    if ( pinMesh.get() != nullptr ) {
        // flow temperature on clad outer surfaces and pellet temperature on clad innner surfaces,
        // and clad inner
        // surface temp on pellet outer surfaces
        auto nodalScalarDOF = AMP::Discretization::simpleDOFManager::create(
            pinMesh, AMP::Mesh::GeomType::Vertex, 1, 1, true );
        thermalMapVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, thermalVariable, true );

        auto pinMeshIDs = pinMesh->getBaseMeshIDs();

        auto pins = AMP::dynamic_pointer_cast<AMP::Mesh::MultiMesh>( pinMesh )->getMeshes();

        for ( const auto &pin : pins ) {
            if ( global_input_db->keyExists( "ThermalNodeToNodeMaps" ) ) {
                auto map =
                    AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::NodeToNodeMap>(
                        pin, global_input_db->getDatabase( "ThermalNodeToNodeMaps" ) );
                for ( size_t j = 0; j < map->getNumberOfOperators(); j++ )
                    n2nColumn->append( map->getOperator( j ) );
            }
            if ( global_input_db->keyExists( "ThermalScalarZAxisMaps" ) ) {
                auto sza =
                    AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::ScalarZAxisMap>(
                        pin, global_input_db->getDatabase( "ThermalScalarZAxisMaps" ) );
                for ( size_t j = 0; j < sza->getNumberOfOperators(); j++ )
                    szaColumn->append( sza->getOperator( j ) );
            }
        }
        if ( n2nColumn->getNumberOfOperators() > 0 )
            mapsColumn->append( n2nColumn );
        if ( szaColumn->getNumberOfOperators() > 0 )
            mapsColumn->append( szaColumn );

        n2nColumn->setVector( thermalMapVec );
        szaColumn->setVector( thermalMapVec );

        int curOperator = 0;
        for ( auto pinMeshID : pinMeshIDs ) {
            auto adapter = manager->Subset( pinMeshID );
            if ( adapter.get() == nullptr )
                continue;

            std::string meshName = adapter->getName();
            std::string prefix;

            if ( meshName.compare( "clad" ) == 0 ) {
                prefix = "Clad";
            } else if ( meshName.compare( "pellet_1" ) == 0 ) {
                prefix = "BottomPellet";
            } else if ( meshName.compare( "pellet_3" ) == 0 ) {
                prefix = "TopPellet";
            } else if ( meshName.compare( 0, 7, "pellet_" ) == 0 ) {
                prefix = "MiddlePellet";
            } else {
                AMP_ERROR( "Unknown Mesh" );
            }

            auto curBVPop = AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
                nonlinearColumnOperator->getOperator( curOperator ) );
            auto curBCcol = AMP::dynamic_pointer_cast<AMP::Operator::ColumnBoundaryOperator>(
                curBVPop->getBoundaryOperator() );
            auto operator_db = global_input_db->getDatabase( prefix + "NonlinearThermalOperator" );
            auto curBCdb =
                global_input_db->getDatabase( operator_db->getString( "BoundaryOperator" ) );
            auto opNames = curBCdb->getStringArray( "boundaryOperators" );
            int opNumber = curBCdb->getInteger( "numberOfBoundaryOperators" );
            for ( int curBCentry = 0; curBCentry != opNumber; curBCentry++ ) {
                if ( opNames[curBCentry] == "P2CRobinVectorCorrection" ) {
                    auto gapBC = AMP::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection>(
                        curBCcol->getBoundaryOperator( curBCentry ) );
                    AMP_ASSERT( thermalMapVec != nullptr );
                    gapBC->setVariableFlux( thermalMapVec );
                    gapBC->reset( gapBC->getOperatorParameters() );
                } else if ( ( opNames[curBCentry] == "BottomP2PNonlinearRobinVectorCorrection" ) ||
                            ( opNames[curBCentry] == "MiddleP2PNonlinearRobinBoundaryCondition" ) ||
                            ( opNames[curBCentry] == "TopP2PNonlinearRobinBoundaryCondition" ) ) {
                    auto p2pBC = AMP::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection>(
                        curBCcol->getBoundaryOperator( curBCentry ) );
                    AMP_ASSERT( thermalMapVec != nullptr );
                    p2pBC->setVariableFlux( thermalMapVec );
                    p2pBC->reset( p2pBC->getOperatorParameters() );
                } else if ( opNames[curBCentry] == "C2WBoundaryVectorCorrection" ) {
                    auto thisDb    = global_input_db->getDatabase( opNames[curBCentry] );
                    bool isCoupled = thisDb->getBoolWithDefault( "IsCoupledBoundary_0", false );
                    if ( isCoupled ) {
                        auto c2wBC =
                            AMP::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection>(
                                curBCcol->getBoundaryOperator( curBCentry ) );
                        AMP_ASSERT( thermalMapVec != nullptr );
                        c2wBC->setVariableFlux( thermalMapVec );
                        c2wBC->setFrozenVector( density_map_vec );
                        c2wBC->setFrozenVector( ChannelDiameterVec );
                        c2wBC->reset( c2wBC->getOperatorParameters() );
                    }
                } else if ( opNames[curBCentry] == "C2PRobinVectorCorrection" ) {
                    auto gapBC = AMP::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection>(
                        curBCcol->getBoundaryOperator( curBCentry ) );
                    AMP_ASSERT( thermalMapVec != nullptr );
                    gapBC->setVariableFlux( thermalMapVec );
                    gapBC->reset( gapBC->getOperatorParameters() );
                } else {
                    AMP_ERROR( "Unknown boundary operator" );
                }
            }
            curOperator++;
        }
    }

    // Create the maps from the clad temperature to the subchannel temperature
    auto cladToSubchannelDb = global_input_db->getDatabase( "CladToSubchannelMaps" );
    auto cladToSubchannelMap =
        AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::CladToSubchannelMap>(
            manager, cladToSubchannelDb );
    cladToSubchannelMap->setVector( subchannelFuelTemp );
    mapsColumn->append( cladToSubchannelMap );

    // Create the maps from the flow variable (enthalpy and pressure) on subchannel mesh and convert
    // to temperature and
    // density then map to clad surface
    auto thermalCladToSubchannelDb = global_input_db->getDatabase( "ThermalSubchannelToCladMaps" );
    auto thermalSubchannelToCladMap =
        AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::SubchannelToCladMap>(
            manager, thermalCladToSubchannelDb );
    auto densityCladToSubchannelDb = global_input_db->getDatabase( "DensitySubchannelToCladMaps" );
    auto densitySubchannelToCladMap =
        AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::SubchannelToCladMap>(
            manager, densityCladToSubchannelDb );
    if ( cladMesh.get() != nullptr ) {
        AMP::LinearAlgebra::VS_Comm commSelector( cladMesh->getComm() );
        auto subsetTheramlVec =
            thermalMapVec->select( commSelector, thermalMapVec->getVariable()->getName() );
        thermalSubchannelToCladMap->setVector( subsetTheramlVec );
        densitySubchannelToCladMap->setVector( density_map_vec );
    }
    auto emptyDb = AMP::make_shared<AMP::InputDatabase>( "empty" );
    emptyDb->putInteger( "print_info_level", 0 );
    auto coupledChannelMapOperatorParams =
        AMP::make_shared<AMP::Operator::CoupledChannelToCladMapOperatorParameters>( emptyDb );
    coupledChannelMapOperatorParams->d_variable               = flowVariable;
    coupledChannelMapOperatorParams->d_vector                 = subchannelFlowTemp;
    coupledChannelMapOperatorParams->d_Mesh                   = subchannelMesh;
    coupledChannelMapOperatorParams->d_thermalMapOperator     = thermalSubchannelToCladMap;
    coupledChannelMapOperatorParams->d_densityMapOperator     = densitySubchannelToCladMap;
    coupledChannelMapOperatorParams->d_subchannelMesh         = subchannelMesh;
    coupledChannelMapOperatorParams->d_subchannelPhysicsModel = subchannelPhysicsModel;
    auto coupledChannelMapOperator =
        AMP::make_shared<AMP::Operator::CoupledChannelToCladMapOperator>(
            coupledChannelMapOperatorParams );
    mapsColumn->append( coupledChannelMapOperator );

    if ( pinMesh.get() != nullptr ) {
        auto copyOp_db = AMP::dynamic_pointer_cast<AMP::InputDatabase>(
            global_input_db->getDatabase( "CopyOperator" ) );
        auto vecCopyOperatorParams =
            AMP::make_shared<AMP::Operator::VectorCopyOperatorParameters>( copyOp_db );
        vecCopyOperatorParams->d_copyVariable = thermalVariable;
        vecCopyOperatorParams->d_copyVector   = thermalMapVec;
        vecCopyOperatorParams->d_Mesh         = pinMesh;
        thermalCopyOperator.reset( new AMP::Operator::VectorCopyOperator( vecCopyOperatorParams ) );
        thermalMapVec->zero();
    }

    auto CoupledOpParams = AMP::make_shared<AMP::Operator::CoupledOperatorParameters>( emptyDb );
    CoupledOpParams->d_CopyOperator = thermalCopyOperator;
    CoupledOpParams->d_MapOperator  = mapsColumn;
    CoupledOpParams->d_BVPOperator  = nonlinearColumnOperator;
    auto nonlinearCoupledOperator =
        AMP::make_shared<AMP::Operator::CoupledOperator>( CoupledOpParams );


    // Create the solution vector
    AMP::LinearAlgebra::Vector::shared_ptr globalSolMultiVector;
    AMP::LinearAlgebra::Vector::shared_ptr specificPowerGpVec;
    createVectors( pinMesh, subchannelMesh, globalSolMultiVector, specificPowerGpVec );

    // Create the rhs and res vectors
    auto globalRhsMultiVector = globalSolMultiVector->cloneVector();
    auto globalResMultiVector = globalSolMultiVector->cloneVector();
    auto flowSolVec           = globalSolMultiVector->subsetVectorForVariable( flowVariable );
    auto flowRhsVec           = globalRhsMultiVector->subsetVectorForVariable( flowVariable );
    auto flowResVec           = globalResMultiVector->subsetVectorForVariable( flowVariable );
    auto globalThermalSolVec  = globalSolMultiVector->subsetVectorForVariable( thermalVariable );
    auto globalThermalRhsVec  = globalRhsMultiVector->subsetVectorForVariable( thermalVariable );
    auto globalThermalResVec  = globalResMultiVector->subsetVectorForVariable( thermalVariable );

    // create nonlinear solver
    AMP::shared_ptr<AMP::Solver::SolverStrategy> nonlinearSolver;
    { // Limit the scope so we can add an if else statement for Petsc vs NOX

        // get the solver databases
        auto nonlinearSolver_db = global_input_db->getDatabase( "NonlinearSolver" );
        auto linearSolver_db    = nonlinearSolver_db->getDatabase( "LinearSolver" );

        // create preconditioner (thermal domains)
        auto columnPreconditioner_db = linearSolver_db->getDatabase( "Preconditioner" );
        auto columnPreconditionerParams =
            AMP::make_shared<AMP::Solver::SolverStrategyParameters>( columnPreconditioner_db );
        columnPreconditionerParams->d_pOperator = linearColumnOperator;
        columnPreconditioner.reset( new AMP::Solver::ColumnSolver( columnPreconditionerParams ) );

        auto trilinosPreconditioner_db =
            columnPreconditioner_db->getDatabase( "TrilinosPreconditioner" );
        unsigned int N_preconditioners = linearColumnOperator->getNumberOfOperators();
        for ( unsigned int id = 0; id < N_preconditioners; id++ ) {
            auto trilinosPreconditionerParams =
                AMP::make_shared<AMP::Solver::SolverStrategyParameters>(
                    trilinosPreconditioner_db );
            trilinosPreconditionerParams->d_pOperator = linearColumnOperator->getOperator( id );
            auto trilinosPreconditioner =
                AMP::make_shared<AMP::Solver::TrilinosMLSolver>( trilinosPreconditionerParams );
            columnPreconditioner->append( trilinosPreconditioner );
        }

        // Create the subchannel preconditioner
        if ( subchannelLinearOperator != nullptr ) {
            auto subchannelPreconditioner_db =
                columnPreconditioner_db->getDatabase( "SubchannelPreconditioner" );
            AMP_ASSERT( subchannelPreconditioner_db != nullptr );
            auto subchannelPreconditionerParams =
                AMP::make_shared<AMP::Solver::SolverStrategyParameters>(
                    subchannelPreconditioner_db );
            subchannelPreconditionerParams->d_pOperator = subchannelLinearOperator;
            auto preconditioner = subchannelPreconditioner_db->getString( "Type" );
            if ( preconditioner == "ML" ) {
                auto subchannelPreconditioner = AMP::make_shared<AMP::Solver::TrilinosMLSolver>(
                    subchannelPreconditionerParams );
                linearColumnOperator->append( subchannelLinearOperator );
                columnPreconditioner->append( subchannelPreconditioner );
            } else if ( preconditioner == "Banded" ) {
                subchannelPreconditioner_db->putInteger( "KL", 3 );
                subchannelPreconditioner_db->putInteger( "KU", 3 );
                auto subchannelPreconditioner =
                    AMP::make_shared<AMP::Solver::BandedSolver>( subchannelPreconditionerParams );
                linearColumnOperator->append( subchannelLinearOperator );
                columnPreconditioner->append( subchannelPreconditioner );
            } else if ( preconditioner == "None" ) {
            } else {
                AMP_ERROR( "Invalid preconditioner type" );
            }
        }

        // create nonlinear solver parameters
        auto nonlinearSolverParams =
            AMP::make_shared<AMP::Solver::TrilinosNOXSolverParameters>( nonlinearSolver_db );
        nonlinearSolverParams->d_comm            = globalComm;
        nonlinearSolverParams->d_pOperator       = nonlinearCoupledOperator;
        nonlinearSolverParams->d_pInitialGuess   = globalSolMultiVector;
        nonlinearSolverParams->d_pLinearOperator = linearColumnOperator;
        nonlinearSolver = AMP::make_shared<AMP::Solver::TrilinosNOXSolver>( nonlinearSolverParams );

        /*// create linear solver
        auto linearSolver =
        AMP::dynamic_pointer_cast<AMP::Solver::PetscSNESSolver>(nonlinearSolver)->getKrylovSolver();
        // set preconditioner
        linearSolver->setPreconditioner(columnPreconditioner);*/
    }

    // don't use zero initial guess
    nonlinearSolver->setZeroInitialGuess( false );


    // Initialize the pin temperatures
    PROFILE_START( "Initialize" );
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    int root_subchannel = -1;
    std::vector<double> range( 6 );
    if ( subchannelMesh.get() != nullptr ) {
        range = subchannelMesh->getBoundingBox();
        AMP_ASSERT( range.size() == 6 );
        if ( subchannelMesh->getComm().getRank() == 0 )
            root_subchannel = globalComm.getRank();
    }
    root_subchannel = globalComm.maxReduce( root_subchannel );
    globalComm.bcast( &range[0], 6, root_subchannel );
    // Desired power of the fuel pin (W)
    double P = global_input_db->getDatabase( "SubchannelTwoEqNonlinearOperator" )
                   ->getDouble( "Rod_Power" );
    // GeomType::Volume of fuel in a 3.81m pin
    if ( pinMesh.get() != nullptr ) {
        const double V = 1.939e-4;
        globalThermalSolVec->setToScalar( 600 );
        auto gaussPtDOFManager = AMP::Discretization::simpleDOFManager::create(
            pinMesh, AMP::Mesh::GeomType::Volume, 1, 8 );
        auto it = pinMesh->getIterator( AMP::Mesh::GeomType::Volume, 0 );
        std::vector<size_t> dofs;
        for ( size_t i = 0; i < it.size(); i++ ) {
            gaussPtDOFManager->getDOFs( it->globalID(), dofs );
            for ( size_t dof : dofs ) {
                specificPowerGpVec->setValueByGlobalID( dof,
                                                        getPower( range, P, V, it->centroid() ) );
            }
            ++it;
        }
        if ( cladMesh.get() != nullptr ) {
            AMP::LinearAlgebra::VS_Mesh meshSelector( cladMesh );
            auto cladPower = specificPowerGpVec->select( meshSelector, "cladPower" );
            cladPower->zero();
        }
        specificPowerGpVec->makeConsistent(
            AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
        volumeIntegralColumnOperator->apply( specificPowerGpVec, globalThermalRhsVec );
    }

    if ( subchannelMesh.get() != nullptr ) {
        // get exit pressure
        double Pout = global_input_db->getDatabase( "SubchannelTwoEqNonlinearOperator" )
                          ->getDouble( "Exit_Pressure" );
        // get inlet temperature
        double Tin = global_input_db->getDatabase( "SubchannelTwoEqNonlinearOperator" )
                         ->getDouble( "Inlet_Temperature" );
        // compute inlet enthalpy
        std::map<std::string, AMP::shared_ptr<std::vector<double>>> enthalpyArgMap;
        enthalpyArgMap.insert(
            std::make_pair( "temperature", AMP::make_shared<std::vector<double>>( 1, Tin ) ) );
        enthalpyArgMap.insert(
            std::make_pair( "pressure", AMP::make_shared<std::vector<double>>( 1, Pout ) ) );
        std::vector<double> enthalpyResult( 1 );
        subchannelPhysicsModel->getProperty( "Enthalpy", enthalpyResult, enthalpyArgMap );
        double hin = enthalpyResult[0];
        std::cout << "Enthalpy Solution:" << hin << std::endl;
        std::cout << "Outlet pressure:" << Pout << std::endl;

        auto subchannelEnthalpy = flowSolVec->select( AMP::LinearAlgebra::VS_Stride( 0, 2 ), "H" );
        auto subchannelPressure = flowSolVec->select( AMP::LinearAlgebra::VS_Stride( 1, 2 ), "P" );

        subchannelEnthalpy->setToScalar( AMP::Operator::Subchannel::scaleEnthalpy * hin );
        subchannelPressure->setToScalar( AMP::Operator::Subchannel::scalePressure * Pout );

        // FIRST APPLY CALL
        auto subchannelLinearParams =
            AMP::dynamic_pointer_cast<AMP::Operator::SubchannelOperatorParameters>(
                subchannelNonlinearOperator->getParameters( "Jacobian", flowSolVec ) );
        subchannelLinearParams->d_initialize = false;
        subchannelLinearOperator->reset( subchannelLinearParams );
        subchannelLinearOperator->residual( flowRhsVec, flowSolVec, flowResVec );
    }

    nonlinearCoupledOperator->residual(
        globalRhsMultiVector, globalSolMultiVector, globalResMultiVector );

    size_t totalOp;
    if ( subchannelMesh != nullptr ) {
        totalOp = nonlinearColumnOperator->getNumberOfOperators() - 1;
    } else {
        totalOp = nonlinearColumnOperator->getNumberOfOperators();
    }
    for ( size_t id = 0; id != totalOp; id++ ) {
        auto nonlinearThermalOperator =
            AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
                nonlinearColumnOperator->getOperator( id ) );
        nonlinearThermalOperator->modifyInitialSolutionVector( globalThermalSolVec );
        nonlinearThermalOperator->modifyRHSvector( globalThermalRhsVec );
    }
    globalThermalRhsVec->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
    PROFILE_STOP( "Initialize" );


    // Solve
    PROFILE_START( "Solve" );
    AMP::pout << "Rhs norm: " << std::setprecision( 13 ) << globalRhsMultiVector->L2Norm()
              << std::endl;
    AMP::pout << "Initial solution norm: " << std::setprecision( 13 )
              << globalSolMultiVector->L2Norm() << std::endl;
    nonlinearCoupledOperator->residual(
        globalRhsMultiVector, globalSolMultiVector, globalResMultiVector );
    double tempResNorm = 0.0;
    double flowResNorm = 0.0;
    if ( pinMesh != nullptr )
        tempResNorm = globalThermalResVec->L2Norm();
    if ( subchannelMesh != nullptr )
        flowResNorm = flowResVec->L2Norm();
    tempResNorm = globalComm.maxReduce( tempResNorm );
    flowResNorm = globalComm.maxReduce( flowResNorm );
    AMP::pout << "Initial residual norm: " << std::setprecision( 13 )
              << globalResMultiVector->L2Norm() << std::endl;
    AMP::pout << "Initial temp residual norm: " << std::setprecision( 13 ) << tempResNorm
              << std::endl;
    AMP::pout << "Initial flow residual norm: " << std::setprecision( 13 ) << flowResNorm
              << std::endl;
    nonlinearSolver->solve( globalRhsMultiVector, globalSolMultiVector );
    nonlinearCoupledOperator->residual(
        globalRhsMultiVector, globalSolMultiVector, globalResMultiVector );
    AMP::pout << "Final residual norm: " << std::setprecision( 13 )
              << globalResMultiVector->L2Norm() << std::endl;
    PROFILE_STOP( "Solve" );


    // Compute the flow temperature and density
    AMP::LinearAlgebra::Vector::shared_ptr flowTempVec;
    AMP::LinearAlgebra::Vector::shared_ptr deltaFlowTempVec;
    AMP::LinearAlgebra::Vector::shared_ptr flowDensityVec;
    if ( subchannelMesh != nullptr ) {
        flowTempVec        = subchannelFuelTemp->cloneVector();
        flowDensityVec     = subchannelFuelTemp->cloneVector();
        int DOFsPerFace[3] = { 0, 0, 2 };
        auto faceDOFManager =
            AMP::Discretization::structuredFaceDOFManager::create( subchannelMesh, DOFsPerFace, 0 );
        DOFsPerFace[2] = 1;
        auto scalarFaceDOFManager =
            AMP::Discretization::structuredFaceDOFManager::create( subchannelMesh, DOFsPerFace, 0 );
        auto face = xyFaceMesh->getIterator( AMP::Mesh::GeomType::Face, 0 );
        std::vector<size_t> dofs;
        std::vector<size_t> scalarDofs;
        // Scale factors to get correct units
        const double h_scale = 1.0 / AMP::Operator::Subchannel::scaleEnthalpy;
        const double P_scale = 1.0 / AMP::Operator::Subchannel::scalePressure;
        for ( size_t i = 0; i < face.size(); i++ ) {
            faceDOFManager->getDOFs( face->globalID(), dofs );
            scalarFaceDOFManager->getDOFs( face->globalID(), scalarDofs );
            std::map<std::string, AMP::shared_ptr<std::vector<double>>> subchannelArgMap;
            auto vec1 = AMP::make_shared<std::vector<double>>(
                1, h_scale * flowSolVec->getValueByGlobalID( dofs[0] ) );
            auto vec2 = AMP::make_shared<std::vector<double>>(
                1, P_scale * flowSolVec->getValueByGlobalID( dofs[1] ) );
            subchannelArgMap.insert( std::make_pair( "enthalpy", vec1 ) );
            subchannelArgMap.insert( std::make_pair( "pressure", vec2 ) );
            std::vector<double> outTemperatureResult( 1 );
            subchannelPhysicsModel->getProperty(
                "Temperature", outTemperatureResult, subchannelArgMap );
            std::vector<double> specificVolume( 1 );
            subchannelPhysicsModel->getProperty(
                "SpecificVolume", specificVolume, subchannelArgMap );
            flowTempVec->setValueByGlobalID( scalarDofs[0], outTemperatureResult[0] );
            flowDensityVec->setValueByGlobalID( scalarDofs[0], 1.0 / specificVolume[0] );
            ++face;
        }
        flowTempVec->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
        double Tin = global_input_db->getDatabase( "SubchannelTwoEqNonlinearOperator" )
                         ->getDouble( "Inlet_Temperature" );
        deltaFlowTempVec = flowTempVec->cloneVector();
        deltaFlowTempVec->copyVector( flowTempVec );
        deltaFlowTempVec->addScalar( deltaFlowTempVec, -Tin );
    }
    double flowTempMin = 1e100;
    double flowTempMax = -1e100;
    if ( flowTempVec != nullptr ) {
        flowTempMin = flowTempVec->min();
        flowTempMax = flowTempVec->max();
    }
    flowTempMin = globalComm.minReduce( flowTempMin );
    flowTempMax = globalComm.maxReduce( flowTempMax );
    AMP::pout << "Subchannel Flow Temp Max : " << flowTempMax << " Min : " << flowTempMin
              << std::endl;


    // Test the subchannel to point map
    auto subchannelToPointMapParams =
        AMP::make_shared<AMP::Operator::SubchannelToPointMapParameters>();
    subchannelToPointMapParams->d_Mesh                   = subchannelMesh;
    subchannelToPointMapParams->d_comm                   = globalComm;
    subchannelToPointMapParams->d_subchannelPhysicsModel = subchannelPhysicsModel;
    subchannelToPointMapParams->d_outputVar.reset( new AMP::LinearAlgebra::Variable( "Density" ) );
    if ( subchannelMesh != nullptr ) {
        auto face = xyFaceMesh->getIterator( AMP::Mesh::GeomType::Face, 0 );
        for ( size_t i = 0; i < face.size(); i++ ) {
            auto pos = face->centroid();
            subchannelToPointMapParams->x.push_back( pos[0] );
            subchannelToPointMapParams->y.push_back( pos[1] );
            subchannelToPointMapParams->z.push_back( pos[2] );
            ++face;
        }
        AMP_ASSERT( subchannelToPointMapParams->x.size() == flowDensityVec->getLocalSize() );
    }
    AMP::Operator::SubchannelToPointMap subchannelDensityToPointMap( subchannelToPointMapParams );
    subchannelToPointMapParams->d_outputVar.reset(
        new AMP::LinearAlgebra::Variable( "Temperature" ) );
    AMP::Operator::SubchannelToPointMap subchannelTemperatureToPointMap(
        subchannelToPointMapParams );
    auto densityMapVec = AMP::LinearAlgebra::SimpleVector<double>::create(
        subchannelToPointMapParams->x.size(), subchannelDensityToPointMap.getOutputVariable() );
    auto temperatureMapVec = AMP::LinearAlgebra::SimpleVector<double>::create(
        subchannelToPointMapParams->x.size(), subchannelTemperatureToPointMap.getOutputVariable() );
    subchannelDensityToPointMap.residual( nullVec, flowSolVec, densityMapVec );
    subchannelTemperatureToPointMap.residual( nullVec, flowSolVec, temperatureMapVec );
    if ( subchannelMesh != nullptr ) {
        auto face = xyFaceMesh->getIterator( AMP::Mesh::GeomType::Face, 0 );
        std::vector<size_t> dofs;
        bool pass_density     = true;
        bool pass_temperature = true;
        for ( size_t i = 0; i < face.size(); i++ ) {
            flowDensityVec->getDOFManager()->getDOFs( face->globalID(), dofs );
            double density1 = flowDensityVec->getValueByGlobalID( dofs[0] );
            double density2 = densityMapVec->getValueByLocalID( i );
            if ( !AMP::Utilities::approx_equal( density1, density2 ) )
                pass_density = false;
            double temp1 = flowTempVec->getValueByGlobalID( dofs[0] );
            double temp2 = temperatureMapVec->getValueByLocalID( i );
            if ( !AMP::Utilities::approx_equal( temp1, temp2 ) )
                pass_temperature = false;
            ++face;
        }
        if ( pass_density )
            ut->passes( "Subchannel density to point map" );
        else
            ut->failure( "Subchannel density to point map" );
        if ( pass_temperature )
            ut->passes( "Subchannel temperature to point map" );
        else
            ut->failure( "Subchannel temperature to point map" );
    }


#ifdef USE_EXT_SILO
    // Rescale the solution to get the correct units
    const double h_scale = 1.0 / AMP::Operator::Subchannel::scaleEnthalpy;
    const double P_scale = 1.0 / AMP::Operator::Subchannel::scalePressure;
    auto enthalpy        = flowSolVec->select( AMP::LinearAlgebra::VS_Stride( 0, 2 ), "H" );
    auto pressure        = flowSolVec->select( AMP::LinearAlgebra::VS_Stride( 1, 2 ), "P" );
    if ( enthalpy.get() != nullptr ) {
        enthalpy->scale( h_scale );
        pressure->scale( P_scale );
    }
    // Register the quantities to plot
    auto siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    if ( xyFaceMesh != nullptr ) {
        siloWriter->registerVector(
            flowSolVec, xyFaceMesh, AMP::Mesh::GeomType::Face, "SubchannelFlow" );
        siloWriter->registerVector(
            flowTempVec, xyFaceMesh, AMP::Mesh::GeomType::Face, "FlowTemp" );
        siloWriter->registerVector(
            deltaFlowTempVec, xyFaceMesh, AMP::Mesh::GeomType::Face, "FlowTempDelta" );
        siloWriter->registerVector(
            flowDensityVec, xyFaceMesh, AMP::Mesh::GeomType::Face, "FlowDensity" );
    }
    if ( pinMesh.get() != nullptr ) {
        siloWriter->registerVector(
            globalThermalSolVec, pinMesh, AMP::Mesh::GeomType::Vertex, "Temperature" );
        siloWriter->registerVector(
            specificPowerGpVec, pinMesh, AMP::Mesh::GeomType::Volume, "Power" );
    }
    siloWriter->writeFile( exeName, 0 );
#endif
    ut->passes( "test runs to completion" );
    globalComm.barrier();
    PROFILE_STOP( "Main" );
    PROFILE_SAVE( "exeName" );
}


int testSubchannelSolveNOX( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    PROFILE_ENABLE( 0 );

    std::string exeName = "testSubchannelSolveNOX-1";
    if ( argc == 2 )
        exeName = argv[1];

    SubchannelSolve( &ut, exeName );

    ut.report();
    PROFILE_SAVE( exeName );

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
