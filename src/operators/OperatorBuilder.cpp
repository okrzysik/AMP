#include "AMP/operators/OperatorBuilder.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/StructuredMeshHelper.h"
#include "AMP/operators/IdentityOperator.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/OperatorFactory.h"
#include "AMP/operators/boundary/ColumnBoundaryOperator.h"
#include "AMP/operators/boundary/DirichletMatrixCorrection.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/boundary/MassMatrixCorrection.h"
#include "AMP/operators/map/MapOperatorParameters.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/VectorBuilder.h"

#ifdef AMP_USE_LIBMESH
    #include "AMP/discretization/structuredFaceDOFManager.h"
    #include "AMP/operators/ElementOperationFactory.h"
    #include "AMP/operators/NeutronicsRhs.h"
    #include "AMP/operators/ParameterFactory.h"
    #include "AMP/operators/boundary/libmesh/NeumannVectorCorrection.h"
    #include "AMP/operators/boundary/libmesh/PressureBoundaryOperator.h"
    #include "AMP/operators/boundary/libmesh/RobinMatrixCorrection.h"
    #include "AMP/operators/boundary/libmesh/RobinVectorCorrection.h"
    #include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
    #include "AMP/operators/diffusion/DiffusionNonlinearFEOperator.h"
    #include "AMP/operators/diffusion/FickSoretNonlinearFEOperator.h"
    #include "AMP/operators/libmesh/MassLinearFEOperator.h"
    #include "AMP/operators/libmesh/VolumeIntegralOperator.h"
    #include "AMP/operators/map/libmesh/MapSurface.h"
    #include "AMP/operators/mechanics/MechanicsConstants.h"
    #include "AMP/operators/mechanics/MechanicsLinearFEOperator.h"
    #include "AMP/operators/mechanics/MechanicsNonlinearFEOperator.h"
    #include "AMP/operators/subchannel/FlowFrapconJacobian.h"
    #include "AMP/operators/subchannel/FlowFrapconOperator.h"
    #include "AMP/operators/subchannel/SubchannelFourEqNonlinearOperator.h"
    #include "AMP/operators/subchannel/SubchannelTwoEqLinearOperator.h"
    #include "AMP/operators/subchannel/SubchannelTwoEqNonlinearOperator.h"
#endif


#include <string>

#include "ProfilerApp.h"


namespace AMP::Operator::OperatorBuilder {


using IdentityOperatorParameters = OperatorParameters;


std::shared_ptr<Operator>
createOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                const std::string &operatorName,
                std::shared_ptr<AMP::Database> input_db,
                std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel );


static void setNestedOperatorMemoryLocations( std::shared_ptr<AMP::Database> input_db,
                                              std::string outerOperatorName,
                                              std::vector<std::string> nestedOperatorNames )
{
    auto outer_db = input_db->getDatabase( outerOperatorName );
    AMP_INSIST( outer_db, "OperatorBuilder: outer DB is null" );

    if ( outer_db->keyExists( "MemoryLocation" ) ) {
        // if outer operator requests a memory location it takes precedent
        auto memLoc = outer_db->getScalar<std::string>( "MemoryLocation" );
        for ( auto &innerName : nestedOperatorNames ) {
            auto inner_db = input_db->getDatabase( innerName );
            AMP_INSIST( inner_db, "OperatorBuilder: inner DB is null" );
            inner_db->putScalar(
                "MemoryLocation", memLoc, Units(), Database::Check::WarnOverwrite );
        }
    } else {
        // outer db does not specify a memory location, check if any internal one does
        // if multiple are specified use most restrictive one
        std::string memLoc{ "device" };
        auto memRestrict = []( std::string m1, std::string m2 ) -> std::string {
            int c1 = 3, c2 = 3;
            if ( m1 == "device" || m1 == "Device" ) {
                c1 = 2;
            } else if ( m1 == "managed" || m1 == "Managed" ) {
                c1 = 1;
            } else if ( m1 == "host" || m1 == "host" ) {
                c1 = 0;
            }
            if ( m2 == "device" || m2 == "Device" ) {
                c2 = 2;
            } else if ( m2 == "managed" || m2 == "Managed" ) {
                c2 = 1;
            } else if ( m2 == "host" || m2 == "host" ) {
                c2 = 0;
            }
            if ( c1 == 3 && c2 == 3 ) {
                // both spaces unrecognized
                AMP_WARNING( "Unrecognized memory space, returning host" );
                return "host";
            } else if ( c1 < c2 ) {
                return m1;
            }
            return m2;
        };
        bool found = false;
        for ( auto &innerName : nestedOperatorNames ) {
            auto inner_db = input_db->getDatabase( innerName );
            AMP_INSIST( inner_db, "OperatorBuilder: inner DB is null" );
            if ( inner_db->keyExists( "MemoryLocation" ) ) {
                found        = true;
                auto memLocI = inner_db->getScalar<std::string>( "MemoryLocation" );
                memLoc       = memRestrict( memLoc, memLocI );
            }
        }
        if ( found ) {
            outer_db->putScalar(
                "MemoryLocation", memLoc, Units(), Database::Check::WarnOverwrite );
            for ( auto &innerName : nestedOperatorNames ) {
                auto inner_db = input_db->getDatabase( innerName );
                AMP_INSIST( inner_db, "OperatorBuilder: inner DB is null" );
                inner_db->putScalar(
                    "MemoryLocation", memLoc, Units(), Database::Check::WarnOverwrite );
            }
        }
    }
}


/********************************************************
 * Create specific operators implementation              *
 ********************************************************/
#ifdef AMP_USE_LIBMESH
static std::vector<std::string> getActiveVariables( std::shared_ptr<const AMP::Database> db,
                                                    const std::string &key )
{
    std::vector<std::string> vars;
    if ( db->isDatabase( key ) ) {
        auto activeDB = db->getDatabase( key );
        for ( auto key2 : activeDB->getAllKeys() )
            vars.push_back( activeDB->getString( key2 ) );
    } else {
        vars = db->getVector<std::string>( key );
    }
    return vars;
}
static std::shared_ptr<SubchannelPhysicsModel>
createSubchannelPhysicsModel( std::shared_ptr<AMP::Database> input_db,
                              std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{
    if ( elementPhysicsModel ) {
        return std::dynamic_pointer_cast<SubchannelPhysicsModel>( elementPhysicsModel );
    } else {
        auto db             = input_db->getDatabase( "SubchannelPhysicsModel" );
        elementPhysicsModel = ElementPhysicsModelFactory::createElementPhysicsModel( db );
        return std::dynamic_pointer_cast<SubchannelPhysicsModel>( elementPhysicsModel );
    }
}
static Operator::shared_ptr
createSubchannelTwoEqLinearOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                     std::shared_ptr<AMP::Database> input_db,
                                     std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{
    AMP_ASSERT( input_db->getString( "name" ) == "SubchannelTwoEqLinearOperator" );
    auto subchannelParams    = std::make_shared<SubchannelOperatorParameters>( input_db );
    subchannelParams->d_Mesh = mesh;
    subchannelParams->d_subchannelPhysicsModel =
        createSubchannelPhysicsModel( input_db, elementPhysicsModel );
    return std::make_shared<SubchannelTwoEqLinearOperator>( subchannelParams );
}
static Operator::shared_ptr
createSubchannelTwoEqNonlinearOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                        std::shared_ptr<AMP::Database> input_db,
                                        std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{
    AMP_ASSERT( input_db->getString( "name" ) == "SubchannelTwoEqNonlinearOperator" );
    auto subchannelParams    = std::make_shared<SubchannelOperatorParameters>( input_db );
    subchannelParams->d_Mesh = mesh;
    subchannelParams->d_subchannelPhysicsModel =
        createSubchannelPhysicsModel( input_db, elementPhysicsModel );
    return std::make_shared<SubchannelTwoEqNonlinearOperator>( subchannelParams );
}
static Operator::shared_ptr
createSubchannelFourEqNonlinearOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                         std::shared_ptr<AMP::Database> input_db,
                                         std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{
    AMP_ASSERT( input_db->getString( "name" ) == "SubchannelFourEqNonlinearOperator" );
    auto subchannelParams    = std::make_shared<SubchannelOperatorParameters>( input_db );
    subchannelParams->d_Mesh = mesh;
    subchannelParams->d_subchannelPhysicsModel =
        createSubchannelPhysicsModel( input_db, elementPhysicsModel );
    return std::make_shared<SubchannelFourEqNonlinearOperator>( subchannelParams );
}

static Operator::shared_ptr
createVolumeIntegralOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                              std::shared_ptr<AMP::Database> input_db,
                              std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{
    PROFILE( "createVolumeIntegralOperator" );

    std::shared_ptr<SourcePhysicsModel> sourcePhysicsModel;
    if ( elementPhysicsModel ) {
        sourcePhysicsModel = std::dynamic_pointer_cast<SourcePhysicsModel>( elementPhysicsModel );
    } else {
        if ( input_db->keyExists( "SourcePhysicsModel" ) ) {
            auto sourceModel_db = input_db->getDatabase( "SourcePhysicsModel" );
            elementPhysicsModel =
                ElementPhysicsModelFactory::createElementPhysicsModel( sourceModel_db );
            sourcePhysicsModel =
                std::dynamic_pointer_cast<SourcePhysicsModel>( elementPhysicsModel );
        }
    }

    // next create a ElementOperation object
    AMP_INSIST( input_db->keyExists( "SourceElement" ), "Key ''SourceElement'' is missing!" );
    auto sourceNonlinearElem =
        ElementOperationFactory::createElementOperation( input_db->getDatabase( "SourceElement" ) );

    // now create the nonlinear source operator
    AMP_ASSERT( input_db->getString( "name" ) == "VolumeIntegralOperator" );
    auto volumeIntegralParameters = std::make_shared<VolumeIntegralOperatorParameters>( input_db );
    volumeIntegralParameters->d_sourcePhysicsModel = sourcePhysicsModel;
    volumeIntegralParameters->d_elemOp             = sourceNonlinearElem;
    volumeIntegralParameters->d_Mesh               = mesh;
    return std::make_shared<VolumeIntegralOperator>( volumeIntegralParameters );
}

static Operator::shared_ptr
createLinearDiffusionOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                               std::shared_ptr<AMP::Database> input_db,
                               std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{
    PROFILE( "createLinearDiffusionOperator" );

    // first create a DiffusionTransportModel
    std::shared_ptr<DiffusionTransportModel> transportModel;
    if ( elementPhysicsModel ) {
        transportModel = std::dynamic_pointer_cast<DiffusionTransportModel>( elementPhysicsModel );
    } else {
        auto db             = input_db->getDatabase( "DiffusionTransportModel" );
        elementPhysicsModel = ElementPhysicsModelFactory::createElementPhysicsModel( db );
        transportModel = std::dynamic_pointer_cast<DiffusionTransportModel>( elementPhysicsModel );
    }
    AMP_INSIST( transportModel, "NULL transport model" );

    // next create a ElementOperation object
    AMP_INSIST( input_db->keyExists( "DiffusionElement" ), "Key ''DiffusionElement'' is missing!" );
    auto diffusionLinElem = ElementOperationFactory::createElementOperation(
        input_db->getDatabase( "DiffusionElement" ) );

    // now create the linear diffusion operator
    AMP_ASSERT( input_db->getString( "name" ) == "DiffusionLinearFEOperator" );
    auto diffusionOpParams = std::make_shared<DiffusionLinearFEOperatorParameters>( input_db );
    diffusionOpParams->d_transportModel = transportModel;
    diffusionOpParams->d_elemOp         = diffusionLinElem;
    diffusionOpParams->d_Mesh           = mesh;
    diffusionOpParams->d_inDofMap       = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    diffusionOpParams->d_outDofMap = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    auto diffusionOp = std::make_shared<DiffusionLinearFEOperator>( diffusionOpParams );

    auto matrix = diffusionOp->getMatrix();
    matrix->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD );

    return diffusionOp;
}

static Operator::shared_ptr
createNonlinearDiffusionOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                  std::shared_ptr<AMP::Database> input_db,
                                  std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{

    // first create a DiffusionTransportModel
    std::shared_ptr<DiffusionTransportModel> transportModel;
    if ( elementPhysicsModel ) {
        transportModel = std::dynamic_pointer_cast<DiffusionTransportModel>( elementPhysicsModel );
    } else {
        auto db             = input_db->getDatabase( "DiffusionTransportModel" );
        elementPhysicsModel = ElementPhysicsModelFactory::createElementPhysicsModel( db );
        transportModel = std::dynamic_pointer_cast<DiffusionTransportModel>( elementPhysicsModel );
    }
    AMP_INSIST( transportModel, "NULL transport model" );

    // next create an ElementOperation object
    AMP_INSIST( input_db->keyExists( "DiffusionElement" ), "Key ''DiffusionElement'' is missing!" );
    auto diffusionNonlinearElem = ElementOperationFactory::createElementOperation(
        input_db->getDatabase( "DiffusionElement" ) );

    // now create the nonlinear diffusion operator parameters
    AMP_ASSERT( input_db->getString( "name" ) == "DiffusionNonlinearFEOperator" );
    auto diffusionNLOpParams = std::make_shared<DiffusionNonlinearFEOperatorParameters>( input_db );
    diffusionNLOpParams->d_transportModel = transportModel;
    diffusionNLOpParams->d_elemOp         = diffusionNonlinearElem;
    diffusionNLOpParams->d_Mesh           = mesh;

    // populate the parameters with frozen active variable vectors

    // nullify vectors in parameters
    diffusionNLOpParams->d_FrozenVecs.clear();

    // create variables and vectors for frozen material inputs
    auto NodalScalarDOF = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    for ( auto name : getActiveVariables( input_db, "ActiveInputVariables" ) ) {
        auto var = std::make_shared<AMP::LinearAlgebra::Variable>( name );

        auto memLoc = AMP::Utilities::memoryLocationFromString(
            diffusionNLOpParams->d_db->getWithDefault<std::string>( "MemoryLocation", "host" ) );
        auto vec = AMP::LinearAlgebra::createVector( NodalScalarDOF, var, true, memLoc );
        if ( input_db->getWithDefault<bool>( "Freeze" + name, false ) )
            diffusionNLOpParams->d_FrozenVecs[name] = vec;
    }

    // create the nonlinear diffusion operator
    return std::make_shared<DiffusionNonlinearFEOperator>( diffusionNLOpParams );
}

static Operator::shared_ptr
createNonlinearFickSoretOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                  std::string operatorName,
                                  std::shared_ptr<AMP::Database> input_db )
{
    AMP_INSIST( input_db, "NULL database object passed" );

    auto operator_db = input_db->getDatabase( operatorName );
    AMP_INSIST( operator_db, "operator database is null" );

    // operator names
    auto fickOperatorName  = operator_db->getString( "FickOperator" );
    auto soretOperatorName = operator_db->getString( "SoretOperator" );

    // Ensure consistency of operator memory locations
    setNestedOperatorMemoryLocations(
        input_db, operatorName, { fickOperatorName, soretOperatorName } );

    auto fickOperator  = createOperator( mesh, fickOperatorName, input_db );
    auto soretOperator = createOperator( mesh, soretOperatorName, input_db );

    auto db        = input_db;
    auto params    = std::make_shared<FickSoretNonlinearFEOperatorParameters>( db );
    params->d_Mesh = mesh;
    params->d_FickOperator =
        std::dynamic_pointer_cast<DiffusionNonlinearFEOperator>( fickOperator );
    params->d_SoretOperator =
        std::dynamic_pointer_cast<DiffusionNonlinearFEOperator>( soretOperator );
    params->d_name = operatorName;
    return std::make_shared<FickSoretNonlinearFEOperator>( params );
}

static Operator::shared_ptr
createLinearMechanicsOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                               std::shared_ptr<AMP::Database> input_db,
                               std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{

    // first create a MechanicsMaterialModel
    if ( !elementPhysicsModel ) {
        AMP_INSIST( input_db->keyExists( "MechanicsMaterialModel" ),
                    "Key ''MechanicsMaterialModel'' is missing!" );

        auto db             = input_db->getDatabase( "MechanicsMaterialModel" );
        elementPhysicsModel = ElementPhysicsModelFactory::createElementPhysicsModel( db );
    }
    AMP_INSIST( elementPhysicsModel, "NULL material model" );

    // next create a ElementOperation object
    AMP_INSIST( input_db->keyExists( "MechanicsElement" ), "Key ''MechanicsElement'' is missing!" );
    auto mechanicsLinElem = ElementOperationFactory::createElementOperation(
        input_db->getDatabase( "MechanicsElement" ) );

    // now create the linear mechanics operator
    AMP_ASSERT( input_db->getString( "name" ) == "MechanicsLinearFEOperator" );
    auto mechanicsOpParams = std::make_shared<MechanicsLinearFEOperatorParameters>( input_db );
    mechanicsOpParams->d_materialModel =
        std::dynamic_pointer_cast<MechanicsMaterialModel>( elementPhysicsModel );
    mechanicsOpParams->d_elemOp   = mechanicsLinElem;
    mechanicsOpParams->d_Mesh     = mesh;
    mechanicsOpParams->d_inDofMap = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 3, true );
    mechanicsOpParams->d_outDofMap = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 3, true );

    return std::make_shared<MechanicsLinearFEOperator>( mechanicsOpParams );
}

static Operator::shared_ptr
createNonlinearMechanicsOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                  std::shared_ptr<AMP::Database> input_db,
                                  std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{

    // first create a MechanicsMaterialModel
    if ( !elementPhysicsModel ) {
        AMP_INSIST( input_db->keyExists( "MechanicsMaterialModel" ),
                    "Key ''MechanicsMaterialModel'' is missing!" );

        auto db             = input_db->getDatabase( "MechanicsMaterialModel" );
        elementPhysicsModel = ElementPhysicsModelFactory::createElementPhysicsModel( db );
    }
    AMP_INSIST( elementPhysicsModel, "NULL material model" );

    // next create a ElementOperation object
    AMP_INSIST( input_db->keyExists( "MechanicsElement" ), "Key ''MechanicsElement'' is missing!" );
    auto mechanicsElem = ElementOperationFactory::createElementOperation(
        input_db->getDatabase( "MechanicsElement" ) );

    // now create the nonlinear mechanics operator
    AMP_ASSERT( input_db->getString( "name" ) == "MechanicsNonlinearFEOperator" );
    auto mechanicsOpParams = std::make_shared<MechanicsNonlinearFEOperatorParameters>( input_db );
    mechanicsOpParams->d_materialModel =
        std::dynamic_pointer_cast<MechanicsMaterialModel>( elementPhysicsModel );
    mechanicsOpParams->d_elemOp = mechanicsElem;
    mechanicsOpParams->d_Mesh   = mesh;
    mechanicsOpParams->d_dofMap[Mechanics::DISPLACEMENT] =
        AMP::Discretization::simpleDOFManager::create(
            mesh, AMP::Mesh::GeomType::Vertex, 1, 3, true );
    mechanicsOpParams->d_dofMap[Mechanics::TEMPERATURE] =
        AMP::Discretization::simpleDOFManager::create(
            mesh, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    mechanicsOpParams->d_dofMap[Mechanics::BURNUP] = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    mechanicsOpParams->d_dofMap[Mechanics::OXYGEN_CONCENTRATION] =
        AMP::Discretization::simpleDOFManager::create(
            mesh, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    mechanicsOpParams->d_dofMap[Mechanics::LHGR] = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 1, true );

    return std::make_shared<MechanicsNonlinearFEOperator>( mechanicsOpParams );
}

static Operator::shared_ptr
createMassLinearFEOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                            std::shared_ptr<AMP::Database> input_db,
                            std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{

    // first create a MassDensityModel
    std::shared_ptr<MassDensityModel> densityModel;
    if ( elementPhysicsModel ) {
        densityModel = std::dynamic_pointer_cast<MassDensityModel>( elementPhysicsModel );
    } else {
        AMP_INSIST( input_db->keyExists( "MassDensityModel" ),
                    "Key ''MassDensityModel'' is missing!" );
        auto densityModel_db = input_db->getDatabase( "MassDensityModel" );
        elementPhysicsModel =
            ElementPhysicsModelFactory::createElementPhysicsModel( densityModel_db );
        densityModel = std::dynamic_pointer_cast<MassDensityModel>( elementPhysicsModel );
    }
    AMP_INSIST( densityModel, "NULL density model" );

    // next create a ElementOperation object
    AMP_INSIST( input_db->keyExists( "MassElement" ), "Key ''MassElement'' is missing!" );
    auto densityLinElem =
        ElementOperationFactory::createElementOperation( input_db->getDatabase( "MassElement" ) );

    // now create the linear density operator
    AMP_ASSERT( input_db->getString( "name" ) == "MassLinearFEOperator" );
    auto densityOpParams            = std::make_shared<MassLinearFEOperatorParameters>( input_db );
    densityOpParams->d_densityModel = densityModel;
    densityOpParams->d_elemOp       = densityLinElem;
    densityOpParams->d_Mesh         = mesh;
    densityOpParams->d_inDofMap     = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    densityOpParams->d_outDofMap = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    return std::make_shared<MassLinearFEOperator>( densityOpParams );
}

static Operator::shared_ptr
createLinearBVPOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                         std::string operatorName,
                         std::shared_ptr<AMP::Database> input_db,
                         std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{
    PROFILE( "createLinearBVPOperator" );

    AMP_INSIST( input_db, "NULL database object passed" );

    auto operator_db = input_db->getDatabase( operatorName );
    AMP_INSIST( operator_db, "NULL database object passed" );

    // names of internal operators
    auto volumeOperatorName   = operator_db->getString( "VolumeOperator" );
    auto boundaryOperatorName = operator_db->getString( "BoundaryOperator" );

    // Ensure consistency of operator memory locations
    setNestedOperatorMemoryLocations(
        input_db, operatorName, { volumeOperatorName, boundaryOperatorName } );

    // create the volume operator
    auto volumeOperator = createOperator( mesh, volumeOperatorName, input_db, elementPhysicsModel );
    auto volumeLinearOp = std::dynamic_pointer_cast<LinearOperator>( volumeOperator );
    AMP_ASSERT( volumeLinearOp );

    // create the boundary operator
    auto boundaryOperator_db = input_db->getDatabase( boundaryOperatorName );
    AMP_INSIST( boundaryOperator_db, "NULL database object passed for boundary operator" );
    boundaryOperator_db->putScalar(
        "isAttachedToVolumeOperator", true, Units(), Database::Check::Overwrite );

    auto boundaryOperator =
        createBoundaryOperator( mesh, boundaryOperatorName, input_db, volumeLinearOp );

    auto bvpOperatorParams                = std::make_shared<BVPOperatorParameters>( operator_db );
    bvpOperatorParams->d_volumeOperator   = volumeOperator;
    bvpOperatorParams->d_boundaryOperator = boundaryOperator;

    return std::make_shared<LinearBVPOperator>( bvpOperatorParams );
}

static Operator::shared_ptr
createNonlinearBVPOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                            std::string operatorName,
                            std::shared_ptr<AMP::Database> input_db,
                            std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{
    AMP_INSIST( input_db, "NULL database object passed" );

    auto operator_db = input_db->getDatabase( operatorName );
    AMP_INSIST( operator_db, "NULL database object passed" );

    // operator names
    auto volumeOperatorName   = operator_db->getString( "VolumeOperator" );
    auto boundaryOperatorName = operator_db->getString( "BoundaryOperator" );

    // Ensure consistency of operator memory locations
    setNestedOperatorMemoryLocations(
        input_db, operatorName, { volumeOperatorName, boundaryOperatorName } );

    // create the volume operator
    auto volumeOperator = createOperator( mesh, volumeOperatorName, input_db, elementPhysicsModel );

    // create the boundary operator
    auto boundaryOperator_db = input_db->getDatabase( boundaryOperatorName );
    AMP_INSIST( boundaryOperator_db, "NULL database object passed for boundary operator" );

    boundaryOperator_db->putScalar(
        "isAttachedToVolumeOperator", true, Units(), Database::Check::Overwrite );


    auto boundaryOperator =
        createBoundaryOperator( mesh, boundaryOperatorName, input_db, volumeOperator );

    auto bvpOperatorParams                = std::make_shared<BVPOperatorParameters>( input_db );
    bvpOperatorParams->d_volumeOperator   = volumeOperator;
    bvpOperatorParams->d_boundaryOperator = boundaryOperator;

    return std::make_shared<NonlinearBVPOperator>( bvpOperatorParams );
}

static std::shared_ptr<BoundaryOperator>
createDirichletMatrixCorrection( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                 std::shared_ptr<AMP::Database> input_db,
                                 Operator::shared_ptr volumeOperator )
{
    auto linearOperator = std::dynamic_pointer_cast<LinearOperator>( volumeOperator );
    auto matrixCorrectionParameters =
        std::make_shared<DirichletMatrixCorrectionParameters>( input_db );
    matrixCorrectionParameters->d_variable    = linearOperator->getOutputVariable();
    matrixCorrectionParameters->d_inputMatrix = linearOperator->getMatrix();
    matrixCorrectionParameters->d_Mesh        = mesh;
    return std::make_shared<DirichletMatrixCorrection>( matrixCorrectionParameters );
}

static std::shared_ptr<BoundaryOperator>
createMassMatrixCorrection( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                            std::shared_ptr<AMP::Database> input_db,
                            Operator::shared_ptr volumeOperator )
{
    auto linearOperator = std::dynamic_pointer_cast<LinearOperator>( volumeOperator );
    auto matrixCorrectionParameters =
        std::make_shared<DirichletMatrixCorrectionParameters>( input_db );
    matrixCorrectionParameters->d_variable    = linearOperator->getOutputVariable();
    matrixCorrectionParameters->d_inputMatrix = linearOperator->getMatrix();
    matrixCorrectionParameters->d_Mesh        = mesh;
    return std::make_shared<MassMatrixCorrection>( matrixCorrectionParameters );
}

static std::shared_ptr<BoundaryOperator>
createRobinMatrixCorrection( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                             std::shared_ptr<AMP::Database> input_db,
                             Operator::shared_ptr volumeOperator,
                             std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{
    auto linearOperator             = std::dynamic_pointer_cast<LinearOperator>( volumeOperator );
    auto matrixCorrectionParameters = std::make_shared<RobinMatrixCorrectionParameters>( input_db );
    matrixCorrectionParameters->d_variable    = linearOperator->getOutputVariable();
    matrixCorrectionParameters->d_inputMatrix = linearOperator->getMatrix();
    matrixCorrectionParameters->d_Mesh        = mesh;

    if ( elementPhysicsModel ) {
        std::shared_ptr<RobinPhysicsModel> robinPhysicsModel;
        robinPhysicsModel = std::dynamic_pointer_cast<RobinPhysicsModel>( elementPhysicsModel );
        matrixCorrectionParameters->d_robinPhysicsModel = robinPhysicsModel;
    }

    return std::make_shared<RobinMatrixCorrection>( matrixCorrectionParameters );
}

static std::shared_ptr<BoundaryOperator>
createRobinVectorCorrection( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                             std::shared_ptr<AMP::Database> input_db,
                             Operator::shared_ptr volumeOperator,
                             std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{
    auto vectorCorrectionParameters =
        std::make_shared<NeumannVectorCorrectionParameters>( input_db );
    vectorCorrectionParameters->d_variable = volumeOperator->getOutputVariable();
    vectorCorrectionParameters->d_Mesh     = mesh;

    if ( elementPhysicsModel ) {
        std::shared_ptr<RobinPhysicsModel> robinPhysicsModel;
        robinPhysicsModel = std::dynamic_pointer_cast<RobinPhysicsModel>( elementPhysicsModel );
        vectorCorrectionParameters->d_robinPhysicsModel = robinPhysicsModel;
    }

    return std::make_shared<RobinVectorCorrection>( vectorCorrectionParameters );
}

static std::shared_ptr<BoundaryOperator>
createNeumannVectorCorrection( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                               std::shared_ptr<AMP::Database> input_db,
                               Operator::shared_ptr volumeOperator,
                               std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{
    auto vectorCorrectionParameters =
        std::make_shared<NeumannVectorCorrectionParameters>( input_db );
    vectorCorrectionParameters->d_variable = volumeOperator->getOutputVariable();
    vectorCorrectionParameters->d_Mesh     = mesh;

    if ( elementPhysicsModel && input_db->isDatabase( "RobinPhysicsModel" ) ) {
        std::shared_ptr<RobinPhysicsModel> robinPhysicsModel;
        robinPhysicsModel = std::dynamic_pointer_cast<RobinPhysicsModel>( elementPhysicsModel );
        vectorCorrectionParameters->d_robinPhysicsModel = robinPhysicsModel;
    }

    return std::make_shared<NeumannVectorCorrection>( vectorCorrectionParameters );
}

#endif


static std::shared_ptr<BoundaryOperator>
createDirichletVectorCorrection( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                 std::shared_ptr<AMP::Database> input_db,
                                 Operator::shared_ptr volumeOperator )
{
    auto vectorCorrectionParameters =
        std::make_shared<DirichletVectorCorrectionParameters>( input_db );
    vectorCorrectionParameters->d_variable = volumeOperator->getOutputVariable();
    vectorCorrectionParameters->d_Mesh     = mesh;
    return std::make_shared<DirichletVectorCorrection>( vectorCorrectionParameters );
}


static std::shared_ptr<BoundaryOperator>
createDirichletVectorCorrection( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                 std::shared_ptr<AMP::Database> input_db )
{
    auto vectorCorrectionParameters =
        std::make_shared<DirichletVectorCorrectionParameters>( input_db );
    vectorCorrectionParameters->d_Mesh = mesh;
    return std::make_shared<DirichletVectorCorrection>( vectorCorrectionParameters );
}


/********************************************************
 * createColumnBoundaryOperator                          *
 ********************************************************/
static std::shared_ptr<BoundaryOperator>
createColumnBoundaryOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                              std::string boundaryOperatorName,
                              std::shared_ptr<AMP::Database> input_db,
                              Operator::shared_ptr volumeOperator )
{

    AMP_INSIST( input_db, "NULL database object passed" );

    auto operator_db = input_db->getDatabase( boundaryOperatorName );
    AMP_INSIST( operator_db,
                "Error: createBoundaryOperator(): "
                "database object with given name not in database" );

    int numberOfBoundaryOperators =
        operator_db->getWithDefault<int>( "numberOfBoundaryOperators", 1 );

    auto boundaryOps = operator_db->getVector<std::string>( "boundaryOperators" );
    AMP_ASSERT( numberOfBoundaryOperators == (int) boundaryOps.size() );

    // Ensure consistency of operator memory locations
    setNestedOperatorMemoryLocations( input_db, boundaryOperatorName, boundaryOps );

    auto params = std::make_shared<OperatorParameters>( operator_db, mesh );

    auto columnBoundaryOperator = std::make_shared<ColumnBoundaryOperator>( params );

    for ( int i = 0; i < numberOfBoundaryOperators; i++ ) {
        auto bcOperator = createBoundaryOperator( mesh, boundaryOps[i], input_db, volumeOperator );
        AMP_ASSERT( bcOperator );
        columnBoundaryOperator->append( bcOperator );
    }

    return columnBoundaryOperator;
}


/********************************************************
 * createOperator                                        *
 ********************************************************/
void createElementPhysicsModel( std::shared_ptr<AMP::Database> operator_db,
                                std::shared_ptr<AMP::Database> input_db,
                                std::shared_ptr<ElementPhysicsModel> &elementPhysicsModel )
{
    if ( !elementPhysicsModel && operator_db->keyExists( "LocalModel" ) ) {
        if ( operator_db->isString( "LocalModel" ) ) {
            auto modelName      = operator_db->getString( "LocalModel" );
            auto model_db       = input_db->getDatabase( modelName );
            elementPhysicsModel = ElementPhysicsModelFactory::createElementPhysicsModel( model_db );
            AMP_ASSERT( elementPhysicsModel );
        } else {
            auto model_db       = input_db->getDatabase( "LocalModel" );
            elementPhysicsModel = ElementPhysicsModelFactory::createElementPhysicsModel( model_db );
            AMP_ASSERT( elementPhysicsModel );
        }
    }
}
std::shared_ptr<Operator> createOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                          const std::string &operatorName,
                                          std::shared_ptr<AMP::Database> input_db,
                                          std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{
    PROFILE( "createOperator" );

    auto operator_db = input_db;
    if ( input_db->keyExists( operatorName ) )
        operator_db = input_db->getDatabase( operatorName );
    AMP_INSIST( operator_db, "operator database is null" );

    // we create the element physics model if a database entry exists
    // and the incoming element physics model pointer is NULL
    createElementPhysicsModel( operator_db, input_db, elementPhysicsModel );

    auto operatorType = operator_db->getString( "name" );
    if ( operatorType == "DirichletMatrixCorrection" ) {
    } else if ( operatorType == "DirichletVectorCorrection" ) {
        return createDirichletVectorCorrection( mesh, operator_db );
#ifdef AMP_USE_LIBMESH
    } else if ( operatorType == "MechanicsLinearFEOperator" ) {
        return createLinearMechanicsOperator( mesh, operator_db, elementPhysicsModel );
    } else if ( operatorType == "MechanicsNonlinearFEOperator" ) {
        return createNonlinearMechanicsOperator( mesh, operator_db, elementPhysicsModel );
    } else if ( operatorType == "DiffusionLinearFEOperator" ) {
        return createLinearDiffusionOperator( mesh, operator_db, elementPhysicsModel );
    } else if ( operatorType == "DiffusionNonlinearFEOperator" ) {
        return createNonlinearDiffusionOperator( mesh, operator_db, elementPhysicsModel );
    } else if ( operatorType == "FickSoretNonlinearFEOperator" ) {
        return createNonlinearFickSoretOperator( mesh, operatorName, input_db );
    } else if ( operatorType == "SubchannelTwoEqLinearOperator" ) {
        return createSubchannelTwoEqLinearOperator( mesh, operator_db, elementPhysicsModel );
    } else if ( operatorType == "SubchannelTwoEqNonlinearOperator" ) {
        return createSubchannelTwoEqNonlinearOperator( mesh, operator_db, elementPhysicsModel );
    } else if ( operatorType == "SubchannelFourEqNonlinearOperator" ) {
        return createSubchannelFourEqNonlinearOperator( mesh, operator_db, elementPhysicsModel );
    } else if ( operatorType == "MassLinearFEOperator" ) {
        return createMassLinearFEOperator( mesh, operator_db, elementPhysicsModel );
    } else if ( operatorType == "VolumeIntegralOperator" ) {
        return createVolumeIntegralOperator( mesh, operator_db, elementPhysicsModel );
    } else if ( operatorType == "LinearBVPOperator" ) {
        // note that we pass in the full database here and not the operator db
        return createLinearBVPOperator( mesh, operatorName, input_db, elementPhysicsModel );
    } else if ( operatorType == "NonlinearBVPOperator" ) {
        // note that we pass in the full database here and not the operator db
        return createNonlinearBVPOperator( mesh, operatorName, input_db, elementPhysicsModel );
#endif
    } else if ( OperatorFactory::exists( operatorType ) ) {
        // Use the OperatorFactory to create the operator
        auto params = std::make_shared<OperatorParameters>( operator_db, mesh );
        return OperatorFactory::create( params );
    }
    AMP_ERROR( "Unable to create operator " + operatorType );
}
std::shared_ptr<Operator> createOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh1,
                                          std::shared_ptr<AMP::Mesh::Mesh> mesh2,
                                          const AMP::AMP_MPI &comm,
                                          std::shared_ptr<AMP::Database> input_db )
{
    auto params       = std::make_shared<MapOperatorParameters>( input_db );
    params->d_Mesh    = mesh1;
    params->d_MapMesh = mesh2;
    params->d_MapComm = comm;
    return OperatorFactory::create( params );
}
std::shared_ptr<Operator> createOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                          const std::string &operatorName,
                                          std::shared_ptr<AMP::Database> input_db )
{
    std::shared_ptr<ElementPhysicsModel> model;
    return createOperator( mesh, operatorName, input_db, model );
}


/********************************************************
 * createBoundaryOperator                                *
 ********************************************************/
std::shared_ptr<BoundaryOperator> createBoundaryOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                                          std::string boundaryOperatorName,
                                                          std::shared_ptr<AMP::Database> input_db,
                                                          Operator::shared_ptr volumeOperator )
{
    AMP_ASSERT( input_db );
    auto operator_db = input_db->getDatabase( boundaryOperatorName );
    AMP_INSIST( operator_db, "operator database is null" );

    // we create the element physics model if a database entry exists
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
    createElementPhysicsModel( operator_db, input_db, elementPhysicsModel );

    auto boundaryType = operator_db->getString( "name" );
    if ( boundaryType == "DirichletVectorCorrection" ) {
        // in this case the volume operator has to be a nonlinear operator
        return createDirichletVectorCorrection( mesh, operator_db, volumeOperator );
#ifdef AMP_USE_LIBMESH
    } else if ( boundaryType == "DirichletMatrixCorrection" ) {
        // in this case the volume operator has to be a linear operator
        return createDirichletMatrixCorrection( mesh, operator_db, volumeOperator );
    } else if ( boundaryType == "MassMatrixCorrection" ) {
        return createMassMatrixCorrection( mesh, operator_db, volumeOperator );
    } else if ( boundaryType == "RobinMatrixCorrection" ) {
        // in this case the volume operator has to be a linear operator
        return createRobinMatrixCorrection(
            mesh, operator_db, volumeOperator, elementPhysicsModel );
    } else if ( boundaryType == "RobinVectorCorrection" ) {
        return createRobinVectorCorrection(
            mesh, operator_db, volumeOperator, elementPhysicsModel );
    } else if ( boundaryType == "NeumannVectorCorrection" ) {
        return createNeumannVectorCorrection(
            mesh, operator_db, volumeOperator, elementPhysicsModel );
#endif
    } else if ( boundaryType == "ColumnBoundaryOperator" ) {
        // note that the global input database is passed here instead of the operator
        // database
        return createColumnBoundaryOperator( mesh, boundaryOperatorName, input_db, volumeOperator );
    } else if ( OperatorFactory::exists( boundaryType ) ) {
        // Use the OperatorFactory to create the operator
        auto params                  = std::make_shared<OperatorParameters>( operator_db, mesh );
        std::shared_ptr<Operator> op = OperatorFactory::create( params );
        auto op2                     = std::dynamic_pointer_cast<BoundaryOperator>( op );
        AMP_ASSERT( op2 );
        return op2;
    }
    return nullptr;
}


} // namespace AMP::Operator::OperatorBuilder
