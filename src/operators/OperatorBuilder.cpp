#include "AMP/operators/OperatorBuilder.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/StructuredMeshHelper.h"
#include "AMP/operators/IdentityOperator.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/VectorBuilder.h"

#ifdef USE_EXT_LIBMESH
    #include "AMP/discretization/structuredFaceDOFManager.h"
    #include "AMP/operators/ElementOperationFactory.h"
    #include "AMP/operators/NeutronicsRhs.h"
    #include "AMP/operators/ParameterFactory.h"
    #include "AMP/operators/boundary/ColumnBoundaryOperator.h"
    #include "AMP/operators/boundary/DirichletMatrixCorrection.h"
    #include "AMP/operators/boundary/DirichletVectorCorrection.h"
    #include "AMP/operators/boundary/MassMatrixCorrection.h"
    #include "AMP/operators/boundary/libmesh/NeumannVectorCorrection.h"
    #include "AMP/operators/boundary/libmesh/PressureBoundaryOperator.h"
    #include "AMP/operators/boundary/libmesh/RobinMatrixCorrection.h"
    #include "AMP/operators/boundary/libmesh/RobinVectorCorrection.h"
    #include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
    #include "AMP/operators/diffusion/DiffusionNonlinearFEOperator.h"
    #include "AMP/operators/diffusion/FickSoretNonlinearFEOperator.h"
    #include "AMP/operators/flow/NavierStokesLSWFFEOperator.h"
    #include "AMP/operators/flow/NavierStokesLSWFLinearFEOperator.h"
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


#define resetOperation( NAME )                                                      \
    do {                                                                            \
        if ( name == #NAME ) {                                                      \
            auto params = std::dynamic_pointer_cast<NAME##Parameters>( in_params ); \
            AMP_ASSERT( params.get() == in_params.get() );                          \
            retOperator.reset( new NAME( params ) );                                \
        }                                                                           \
    } while ( 0 )


namespace AMP::Operator {


using IdentityOperatorParameters = OperatorParameters;


std::shared_ptr<Operator>
OperatorBuilder::createOperator( std::shared_ptr<OperatorParameters> in_params )
{
    std::shared_ptr<Operator> retOperator;

    AMP_INSIST( in_params, "ERROR: OperatorBuilder::createOperator has NULL input" );
    AMP_INSIST( in_params->d_db,
                "ERROR: OperatorBuilder::createOperator has NULL database pointer in in_params" );

    std::string name = in_params->d_db->getString( "name" );

    resetOperation( IdentityOperator );
#ifdef USE_EXT_LIBMESH
    resetOperation( DirichletMatrixCorrection );
    resetOperation( DirichletVectorCorrection );
    resetOperation( NeumannVectorCorrection );
    resetOperation( RobinMatrixCorrection );
    resetOperation( RobinVectorCorrection );
    resetOperation( MechanicsLinearFEOperator );
    resetOperation( MechanicsNonlinearFEOperator );
    resetOperation( DiffusionLinearFEOperator );
    resetOperation( DiffusionNonlinearFEOperator );
    resetOperation( FickSoretNonlinearFEOperator );
    resetOperation( FlowFrapconOperator );
    resetOperation( FlowFrapconJacobian );
    resetOperation( NeutronicsRhs );
    // resetOperation(Mesh3Dto1D);
    if ( name == "PressureBoundaryOperator" ) {
        auto params = std::dynamic_pointer_cast<TractionBoundaryOperatorParameters>( in_params );
        AMP_ASSERT( params.get() == in_params.get() );
        retOperator.reset( new PressureBoundaryOperator( params ) );
    }
#endif
    if ( ( name == "LinearBVPOperator" ) || ( name == "NonlinearBVPOperator" ) ) {
        auto bvpOperatorParameters = std::dynamic_pointer_cast<BVPOperatorParameters>( in_params );

        AMP_INSIST( bvpOperatorParameters, "ERROR: NULL BVPOperatorParameters passed" );
        AMP_INSIST( bvpOperatorParameters->d_volumeOperatorParams,
                    "ERROR: BVPOperatorParameters has NULL volumeOperatorParams pointer" );
        AMP_INSIST( bvpOperatorParameters->d_boundaryOperatorParams,
                    "ERROR: BVPOperatorParameters has NULL boundaryOperatorParams pointer" );

        bvpOperatorParameters->d_volumeOperator =
            OperatorBuilder::createOperator( bvpOperatorParameters->d_volumeOperatorParams );
        bvpOperatorParameters->d_boundaryOperator = std::dynamic_pointer_cast<BoundaryOperator>(
            OperatorBuilder::createOperator( bvpOperatorParameters->d_boundaryOperatorParams ) );

        if ( name == "LinearBVPOperator" )
            retOperator.reset( new LinearBVPOperator(
                std::dynamic_pointer_cast<const BVPOperatorParameters>( in_params ) ) );
        else
            retOperator.reset( new NonlinearBVPOperator(
                std::dynamic_pointer_cast<const BVPOperatorParameters>( in_params ) ) );
    }

    return retOperator;
}


std::shared_ptr<Operator>
OperatorBuilder::createOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                 std::string operatorName,
                                 std::shared_ptr<AMP::Database> tmp_input_db,
                                 std::shared_ptr<ElementPhysicsModel> elementPhysicsModel,
                                 std::shared_ptr<ElementPhysicsModelFactory> localModelFactory )
{

    std::shared_ptr<Operator> retOperator;

    auto input_db = tmp_input_db;

    auto operator_db = input_db->getDatabase( operatorName );

    AMP_INSIST( operator_db,
                "Error:: OperatorBuilder::createOperator(): No operator database entry with "
                "given name exists in input database" );

    // we create the element physics model if a database entry exists
    // and the incoming element physics model pointer is NULL
    if ( ( elementPhysicsModel.get() == nullptr ) && ( operator_db->keyExists( "LocalModel" ) ) ) {
        // extract the name of the local model from the operator database
        auto localModelName = operator_db->getString( "LocalModel" );
        // check whether a database exists in the global database
        // (NOTE: not the operator database) with the given name
        AMP_INSIST( input_db->keyExists( localModelName ),
                    "Error:: OperatorBuilder::createOperator(): No local model "
                    "database entry with given name exists in input database" );

        auto localModel_db = input_db->getDatabase( localModelName );
        AMP_INSIST( localModel_db,
                    "Error:: OperatorBuilder::createOperator(): No local model database "
                    "entry with given name exists in input databaseot" );

        // If a non-NULL factory is being supplied through the argument list
        // use it, else call the AMP ElementPhysicsModelFactory interface
        if ( localModelFactory ) {
            elementPhysicsModel = localModelFactory->createElementPhysicsModel( localModel_db );
        } else {
            elementPhysicsModel =
                ElementPhysicsModelFactory::createElementPhysicsModel( localModel_db );
        }

        AMP_INSIST( elementPhysicsModel,
                    "Error:: OperatorBuilder::createOperator(): local model creation failed" );
    }

    auto operatorType = operator_db->getString( "name" );

    if ( operatorType == "IdentityOperator" ) {
        retOperator = OperatorBuilder::createIdentityOperator( meshAdapter, operator_db );
    } else if ( operatorType == "MechanicsLinearFEOperator" ) {
        retOperator = OperatorBuilder::createLinearMechanicsOperator(
            meshAdapter, operator_db, elementPhysicsModel );
    } else if ( operatorType == "MechanicsNonlinearFEOperator" ) {
        retOperator = OperatorBuilder::createNonlinearMechanicsOperator(
            meshAdapter, operator_db, elementPhysicsModel );
    } else if ( operatorType == "DiffusionLinearFEOperator" ) {
        retOperator = OperatorBuilder::createLinearDiffusionOperator(
            meshAdapter, operator_db, elementPhysicsModel );
    } else if ( operatorType == "DiffusionNonlinearFEOperator" ) {
        retOperator = OperatorBuilder::createNonlinearDiffusionOperator(
            meshAdapter, operator_db, elementPhysicsModel );
    } else if ( operatorType == "NavierStokesLSWFLinearFEOperator" ) {
        retOperator = OperatorBuilder::createLinearNavierStokesLSWFOperator(
            meshAdapter, operator_db, elementPhysicsModel );
    } else if ( operatorType == "NavierStokesLSWFFEOperator" ) {
        retOperator = OperatorBuilder::createNonlinearNavierStokesLSWFOperator(
            meshAdapter, operator_db, elementPhysicsModel );
    } else if ( operatorType == "FickSoretNonlinearFEOperator" ) {
        retOperator = OperatorBuilder::createNonlinearFickSoretOperator(
            meshAdapter, operatorName, input_db, elementPhysicsModel, localModelFactory );
    } else if ( operatorType == "FlowFrapconOperator" ) {
        retOperator = OperatorBuilder::createFlowFrapconOperator( meshAdapter, operator_db );
    } else if ( operatorType == "FlowFrapconJacobian" ) {
        retOperator = OperatorBuilder::createFlowFrapconJacobian( meshAdapter, operator_db );
    } else if ( operatorType == "SubchannelTwoEqLinearOperator" ) {
        retOperator = OperatorBuilder::createSubchannelTwoEqLinearOperator(
            meshAdapter, operator_db, elementPhysicsModel );
    } else if ( operatorType == "SubchannelTwoEqNonlinearOperator" ) {
        retOperator = OperatorBuilder::createSubchannelTwoEqNonlinearOperator(
            meshAdapter, operator_db, elementPhysicsModel );
    } else if ( operatorType == "SubchannelFourEqNonlinearOperator" ) {
        retOperator = OperatorBuilder::createSubchannelFourEqNonlinearOperator(
            meshAdapter, operator_db, elementPhysicsModel );
    } else if ( operatorType == "NeutronicsRhsOperator" ) {
        retOperator = OperatorBuilder::createNeutronicsRhsOperator( meshAdapter, operator_db );
    } else if ( operatorType == "MassLinearFEOperator" ) {
        retOperator = OperatorBuilder::createMassLinearFEOperator(
            meshAdapter, operator_db, elementPhysicsModel );
    } else if ( operatorType == "VolumeIntegralOperator" ) {
        retOperator = OperatorBuilder::createVolumeIntegralOperator(
            meshAdapter, operator_db, elementPhysicsModel );
    } else if ( operatorType == "LinearBVPOperator" ) {
        // note that we pass in the full database here and not the operator db
        retOperator = OperatorBuilder::createLinearBVPOperator(
            meshAdapter, operatorName, input_db, elementPhysicsModel, localModelFactory );
    } else if ( operatorType == "NonlinearBVPOperator" ) {
        // note that we pass in the full database here and not the operator db
        retOperator = OperatorBuilder::createNonlinearBVPOperator(
            meshAdapter, operatorName, input_db, elementPhysicsModel, localModelFactory );
    } else if ( operatorType == "DirichletMatrixCorrection" ) {
    } else if ( operatorType == "DirichletVectorCorrection" ) {
        retOperator = OperatorBuilder::createDirichletVectorCorrection(
            meshAdapter, operator_db, elementPhysicsModel );
    } else if ( operatorType == "PressureBoundaryOperator" ) {
        retOperator = OperatorBuilder::createPressureBoundaryOperator(
            meshAdapter, operator_db, elementPhysicsModel );
    } else if ( operatorType == "NeumannVectorCorrection" ) {
    } else if ( operatorType == "RobinMatrixCorrection" ) {
    } else {
        AMP_ERROR( "Unknown operator: " + operatorName );
    }

    return retOperator;
}


#ifdef USE_EXT_LIBMESH


// Create the identity operator
Operator::shared_ptr
OperatorBuilder::createIdentityOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                         std::shared_ptr<AMP::Database> input_db )
{
    AMP_INSIST( input_db, "Error: The database object for SubchannelTwoEqLinearOperator is NULL" );

    auto params       = std::make_shared<OperatorParameters>( input_db );
    params->d_Mesh    = meshAdapter;
    auto subchannelOp = std::make_shared<IdentityOperator>( params );

    return subchannelOp;
}


// Create the FlowFrapconOperator
Operator::shared_ptr
OperatorBuilder::createFlowFrapconOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                            std::shared_ptr<AMP::Database> input_db )
{

    // now create the flow frapcon operator
    std::shared_ptr<AMP::Database> flowOp_db;
    if ( input_db->getString( "name" ) == "FlowFrapconOperator" ) {
        flowOp_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    AMP_INSIST( flowOp_db, "Error: The database object for FlowFrapconOperator is NULL" );

    auto flowOpParams    = std::make_shared<FlowFrapconOperatorParameters>( flowOp_db );
    flowOpParams->d_Mesh = meshAdapter;
    auto flowOp          = std::make_shared<FlowFrapconOperator>( flowOpParams );

    return flowOp;
}

Operator::shared_ptr OperatorBuilder::createSubchannelTwoEqLinearOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{
    // first create a SubchannelPhysicsModel
    std::shared_ptr<SubchannelPhysicsModel> transportModel;

    if ( elementPhysicsModel ) {
        transportModel = std::dynamic_pointer_cast<SubchannelPhysicsModel>( elementPhysicsModel );
    } else {
        std::shared_ptr<AMP::Database> transportModel_db;
        if ( input_db->keyExists( "SubchannelPhysicsModel" ) ) {
            transportModel_db = input_db->getDatabase( "SubchannelPhysicsModel" );
        } else {
            AMP_INSIST( false, "Key ''SubchannelPhysicsModel'' is missing!" );
        }
        elementPhysicsModel =
            ElementPhysicsModelFactory::createElementPhysicsModel( transportModel_db );
        transportModel = std::dynamic_pointer_cast<SubchannelPhysicsModel>( elementPhysicsModel );
    }

    AMP_INSIST( transportModel, "NULL transport model" );
    // create the operator
    std::shared_ptr<AMP::Database> subchannel_db;
    if ( input_db->getString( "name" ) == "SubchannelTwoEqLinearOperator" ) {
        subchannel_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    AMP_INSIST( subchannel_db,
                "Error: The database object for SubchannelTwoEqLinearOperator is NULL" );

    auto subchannelParams    = std::make_shared<SubchannelOperatorParameters>( subchannel_db );
    subchannelParams->d_Mesh = meshAdapter;
    subchannelParams->d_subchannelPhysicsModel = transportModel;

    auto subchannelOp = std::make_shared<SubchannelTwoEqLinearOperator>( subchannelParams );

    return subchannelOp;
}


Operator::shared_ptr OperatorBuilder::createSubchannelTwoEqNonlinearOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{

    // first create a SubchannelPhysicsModel
    std::shared_ptr<SubchannelPhysicsModel> transportModel;

    if ( elementPhysicsModel ) {
        transportModel = std::dynamic_pointer_cast<SubchannelPhysicsModel>( elementPhysicsModel );
    } else {
        std::shared_ptr<AMP::Database> transportModel_db;
        if ( input_db->keyExists( "SubchannelPhysicsModel" ) ) {
            transportModel_db = input_db->getDatabase( "SubchannelPhysicsModel" );
        } else {
            AMP_INSIST( false, "Key ''SubchannelPhysicsModel'' is missing!" );
        }
        elementPhysicsModel =
            ElementPhysicsModelFactory::createElementPhysicsModel( transportModel_db );
        transportModel = std::dynamic_pointer_cast<SubchannelPhysicsModel>( elementPhysicsModel );
    }

    AMP_INSIST( transportModel, "NULL transport model" );

    // create the operator
    std::shared_ptr<AMP::Database> subchannel_db;
    if ( input_db->getString( "name" ) == "SubchannelTwoEqNonlinearOperator" ) {
        subchannel_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    AMP_INSIST( subchannel_db,
                "Error: The database object for SubchannelTwoEqNonlinearOperator is NULL" );

    auto subchannelParams    = std::make_shared<SubchannelOperatorParameters>( subchannel_db );
    subchannelParams->d_Mesh = meshAdapter;
    subchannelParams->d_subchannelPhysicsModel = transportModel;
    auto subchannelOp = std::make_shared<SubchannelTwoEqNonlinearOperator>( subchannelParams );

    return subchannelOp;
}


Operator::shared_ptr OperatorBuilder::createSubchannelFourEqNonlinearOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{

    // first create a SubchannelPhysicsModel
    std::shared_ptr<SubchannelPhysicsModel> transportModel;

    if ( elementPhysicsModel ) {
        transportModel = std::dynamic_pointer_cast<SubchannelPhysicsModel>( elementPhysicsModel );
    } else {
        std::shared_ptr<AMP::Database> transportModel_db;
        if ( input_db->keyExists( "SubchannelPhysicsModel" ) ) {
            transportModel_db = input_db->getDatabase( "SubchannelPhysicsModel" );
        } else {
            AMP_INSIST( false, "Key ''SubchannelPhysicsModel'' is missing!" );
        }
        elementPhysicsModel =
            ElementPhysicsModelFactory::createElementPhysicsModel( transportModel_db );
        transportModel = std::dynamic_pointer_cast<SubchannelPhysicsModel>( elementPhysicsModel );
    }

    AMP_INSIST( transportModel, "NULL transport model" );

    // create the operator
    std::shared_ptr<AMP::Database> subchannel_db;
    if ( input_db->getString( "name" ) == "SubchannelFourEqNonlinearOperator" ) {
        subchannel_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    AMP_INSIST( subchannel_db,
                "Error: The database object for SubchannelFourEqNonlinearOperator is NULL" );

    auto subchannelParams    = std::make_shared<SubchannelOperatorParameters>( subchannel_db );
    subchannelParams->d_Mesh = meshAdapter;
    subchannelParams->d_subchannelPhysicsModel = transportModel;
    auto subchannelOp = std::make_shared<SubchannelFourEqNonlinearOperator>( subchannelParams );

    return subchannelOp;
}


Operator::shared_ptr
OperatorBuilder::createNeutronicsRhsOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                              std::shared_ptr<AMP::Database> input_db )
{

    // now create the Neutronics operator
    std::shared_ptr<AMP::Database> NeutronicsOp_db;
    if ( input_db->getString( "name" ) == "NeutronicsRhsOperator" ) {
        NeutronicsOp_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    AMP_INSIST( NeutronicsOp_db,
                "Error: The database object for Neutronics Source Operator is NULL" );

    auto neutronicsOpParams    = std::make_shared<NeutronicsRhsParameters>( NeutronicsOp_db );
    neutronicsOpParams->d_Mesh = meshAdapter;
    auto neutronicsOp          = std::make_shared<NeutronicsRhs>( neutronicsOpParams );

    return neutronicsOp;
}


Operator::shared_ptr OperatorBuilder::createLinearDiffusionOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{

    // first create a DiffusionTransportModel
    std::shared_ptr<DiffusionTransportModel> transportModel;

    if ( elementPhysicsModel ) {
        transportModel = std::dynamic_pointer_cast<DiffusionTransportModel>( elementPhysicsModel );
    } else {
        std::shared_ptr<AMP::Database> transportModel_db;
        if ( input_db->keyExists( "DiffusionTransportModel" ) ) {
            transportModel_db = input_db->getDatabase( "DiffusionTransportModel" );
        } else {
            AMP_INSIST( false, "Key ''DiffusionTransportModel'' is missing!" );
        }
        elementPhysicsModel =
            ElementPhysicsModelFactory::createElementPhysicsModel( transportModel_db );
        transportModel = std::dynamic_pointer_cast<DiffusionTransportModel>( elementPhysicsModel );
    }

    AMP_INSIST( transportModel, "NULL transport model" );

    // next create a ElementOperation object
    AMP_INSIST( input_db->keyExists( "DiffusionElement" ), "Key ''DiffusionElement'' is missing!" );
    std::shared_ptr<ElementOperation> diffusionLinElem =
        ElementOperationFactory::createElementOperation(
            input_db->getDatabase( "DiffusionElement" ) );

    // now create the linear diffusion operator
    std::shared_ptr<AMP::Database> diffusionLinFEOp_db;
    if ( input_db->getString( "name" ) == "DiffusionLinearFEOperator" ) {
        diffusionLinFEOp_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    AMP_INSIST( diffusionLinFEOp_db,
                "Error: The database object for DiffusionLinearFEOperator is NULL" );

    auto diffusionOpParams =
        std::make_shared<DiffusionLinearFEOperatorParameters>( diffusionLinFEOp_db );
    diffusionOpParams->d_transportModel = transportModel;
    diffusionOpParams->d_elemOp         = diffusionLinElem;
    diffusionOpParams->d_Mesh           = meshAdapter;
    diffusionOpParams->d_inDofMap       = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    diffusionOpParams->d_outDofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    auto diffusionOp = std::make_shared<DiffusionLinearFEOperator>( diffusionOpParams );

    AMP::LinearAlgebra::Matrix::shared_ptr matrix = diffusionOp->getMatrix();
    matrix->makeConsistent();

    return diffusionOp;
}


Operator::shared_ptr OperatorBuilder::createVolumeIntegralOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{
    std::shared_ptr<SourcePhysicsModel> sourcePhysicsModel;

    if ( elementPhysicsModel ) {
        sourcePhysicsModel = std::dynamic_pointer_cast<SourcePhysicsModel>( elementPhysicsModel );
    } else {
        if ( input_db->keyExists( "SourcePhysicsModel" ) ) {
            std::shared_ptr<AMP::Database> sourceModel_db =
                input_db->getDatabase( "SourcePhysicsModel" );
            elementPhysicsModel =
                ElementPhysicsModelFactory::createElementPhysicsModel( sourceModel_db );
            sourcePhysicsModel =
                std::dynamic_pointer_cast<SourcePhysicsModel>( elementPhysicsModel );
        }
    }

    // next create a ElementOperation object
    AMP_INSIST( input_db->keyExists( "SourceElement" ), "Key ''SourceElement'' is missing!" );
    std::shared_ptr<ElementOperation> sourceNonlinearElem =
        ElementOperationFactory::createElementOperation( input_db->getDatabase( "SourceElement" ) );

    // now create the nonlinear source operator
    if ( input_db->getString( "name" ) != "VolumeIntegralOperator" ) {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    auto volumeIntegralParameters = std::make_shared<VolumeIntegralOperatorParameters>( input_db );
    volumeIntegralParameters->d_sourcePhysicsModel = sourcePhysicsModel;
    volumeIntegralParameters->d_elemOp             = sourceNonlinearElem;
    volumeIntegralParameters->d_Mesh               = meshAdapter;
    auto nonlinearSourceOp = std::make_shared<VolumeIntegralOperator>( volumeIntegralParameters );

    return nonlinearSourceOp;
}


Operator::shared_ptr OperatorBuilder::createNonlinearDiffusionOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{

    // first create a DiffusionTransportModel
    std::shared_ptr<DiffusionTransportModel> transportModel;

    if ( elementPhysicsModel ) {
        transportModel = std::dynamic_pointer_cast<DiffusionTransportModel>( elementPhysicsModel );
    } else {
        std::shared_ptr<AMP::Database> transportModel_db;
        if ( input_db->keyExists( "DiffusionTransportModel" ) ) {
            transportModel_db = input_db->getDatabase( "DiffusionTransportModel" );
        } else {
            AMP_INSIST( false, "Key ''DiffusionTransportModel'' is missing!" );
        }
        elementPhysicsModel =
            ElementPhysicsModelFactory::createElementPhysicsModel( transportModel_db );
        transportModel = std::dynamic_pointer_cast<DiffusionTransportModel>( elementPhysicsModel );
    }

    AMP_INSIST( transportModel, "NULL transport model" );

    // next create an ElementOperation object
    AMP_INSIST( input_db->keyExists( "DiffusionElement" ), "Key ''DiffusionElement'' is missing!" );
    std::shared_ptr<ElementOperation> diffusionNonlinearElem =
        ElementOperationFactory::createElementOperation(
            input_db->getDatabase( "DiffusionElement" ) );

    // now create the nonlinear diffusion operator parameters
    std::shared_ptr<AMP::Database> diffusionNLinFEOp_db;
    if ( input_db->getString( "name" ) == "DiffusionNonlinearFEOperator" ) {
        diffusionNLinFEOp_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }
    AMP_INSIST( diffusionNLinFEOp_db,
                "Error: The database object for DiffusionNonlinearFEOperator is NULL" );
    auto diffusionNLOpParams =
        std::make_shared<DiffusionNonlinearFEOperatorParameters>( diffusionNLinFEOp_db );
    diffusionNLOpParams->d_transportModel = transportModel;
    diffusionNLOpParams->d_elemOp         = diffusionNonlinearElem;
    diffusionNLOpParams->d_Mesh           = meshAdapter;

    // populate the parameters with frozen active variable vectors

    // nullify vectors in parameters
    diffusionNLOpParams->d_FrozenTemperature.reset();
    diffusionNLOpParams->d_FrozenConcentration.reset();
    diffusionNLOpParams->d_FrozenBurnup.reset();

    // create variables and vectors for frozen material inputs
    auto active_db      = diffusionNLinFEOp_db->getDatabase( "ActiveInputVariables" );
    auto NodalScalarDOF = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    std::string name;
    std::shared_ptr<AMP::LinearAlgebra::Variable> tVar;
    AMP::LinearAlgebra::Vector::shared_ptr tVec;
    std::shared_ptr<AMP::LinearAlgebra::Variable> cVar;
    AMP::LinearAlgebra::Vector::shared_ptr cVec;
    std::shared_ptr<AMP::LinearAlgebra::Variable> bVar;
    AMP::LinearAlgebra::Vector::shared_ptr bVec;
    name = active_db->getWithDefault<std::string>( "Temperature", "not_specified" );
    if ( name != "not_specified" ) {
        tVar.reset( new AMP::LinearAlgebra::Variable( name ) );
        tVec = AMP::LinearAlgebra::createVector( NodalScalarDOF, tVar, true );
        if ( diffusionNLinFEOp_db->getWithDefault<bool>( "FreezeTemperature", false ) )
            diffusionNLOpParams->d_FrozenTemperature = tVec;
    }
    name = active_db->getWithDefault<std::string>( "Concentration", "not_specified" );
    if ( name != "not_specified" ) {
        cVar.reset( new AMP::LinearAlgebra::Variable( name ) );
        cVec = AMP::LinearAlgebra::createVector( NodalScalarDOF, cVar, true );
        if ( diffusionNLinFEOp_db->getWithDefault<bool>( "FreezeConcentration", false ) )
            diffusionNLOpParams->d_FrozenConcentration = cVec;
    }
    name = active_db->getWithDefault<std::string>( "Burnup", "not_specified" );
    if ( name != "not_specified" ) {
        bVar.reset( new AMP::LinearAlgebra::Variable( name ) );
        bVec = AMP::LinearAlgebra::createVector( NodalScalarDOF, bVar, true );
        if ( diffusionNLinFEOp_db->getWithDefault<bool>( "FreezeBurnup", false ) )
            diffusionNLOpParams->d_FrozenBurnup = bVec;
    }

    // create the nonlinear diffusion operator
    auto nonlinearDiffusionOp =
        std::make_shared<DiffusionNonlinearFEOperator>( diffusionNLOpParams );

    return nonlinearDiffusionOp;
}


Operator::shared_ptr OperatorBuilder::createNonlinearFickSoretOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    std::string operatorName,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel>,
    std::shared_ptr<ElementPhysicsModelFactory> localModelFactory )
{
    std::shared_ptr<Operator> retOperator;
    AMP_INSIST( input_db, "NULL database object passed" );

    auto operator_db = input_db->getDatabase( operatorName );
    AMP_INSIST( operator_db, "NULL database object passed" );

    std::string fickOperatorName  = operator_db->getString( "FickOperator" );
    std::string soretOperatorName = operator_db->getString( "SoretOperator" );

    std::shared_ptr<ElementPhysicsModel> fickPhysicsModel;
    std::shared_ptr<ElementPhysicsModel> soretPhysicsModel;


    Operator::shared_ptr fickOperator = OperatorBuilder::createOperator(
        meshAdapter, fickOperatorName, input_db, fickPhysicsModel, localModelFactory );
    AMP_INSIST(
        fickOperator,
        "Error: unable to create Fick operator in OperatorBuilder::createFickSoretOperator" );

    Operator::shared_ptr soretOperator = OperatorBuilder::createOperator(
        meshAdapter, soretOperatorName, input_db, soretPhysicsModel, localModelFactory );

    AMP_INSIST(
        soretOperator,
        "Error: unable to create Soret operator in OperatorBuilder::createFickSoretOperator" );

    auto db        = input_db;
    auto params    = std::make_shared<FickSoretNonlinearFEOperatorParameters>( db );
    params->d_Mesh = meshAdapter;
    params->d_FickOperator =
        std::dynamic_pointer_cast<DiffusionNonlinearFEOperator>( fickOperator );
    params->d_SoretOperator =
        std::dynamic_pointer_cast<DiffusionNonlinearFEOperator>( soretOperator );
    params->d_name = operatorName;
    auto fsOp      = std::make_shared<FickSoretNonlinearFEOperator>( params );

    return fsOp;
}


Operator::shared_ptr OperatorBuilder::createLinearMechanicsOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{

    // first create a MechanicsMaterialModel
    if ( elementPhysicsModel.get() == nullptr ) {
        AMP_INSIST( input_db->keyExists( "MechanicsMaterialModel" ),
                    "Key ''MechanicsMaterialModel'' is missing!" );

        auto materialModel_db = input_db->getDatabase( "MechanicsMaterialModel" );
        elementPhysicsModel =
            ElementPhysicsModelFactory::createElementPhysicsModel( materialModel_db );
    }

    AMP_INSIST( elementPhysicsModel, "NULL material model" );

    // next create a ElementOperation object
    AMP_INSIST( input_db->keyExists( "MechanicsElement" ), "Key ''MechanicsElement'' is missing!" );
    auto mechanicsLinElem = ElementOperationFactory::createElementOperation(
        input_db->getDatabase( "MechanicsElement" ) );

    // now create the linear mechanics operator
    std::shared_ptr<AMP::Database> mechanicsLinFEOp_db;
    if ( input_db->getString( "name" ) == "MechanicsLinearFEOperator" ) {
        mechanicsLinFEOp_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    AMP_INSIST( mechanicsLinFEOp_db,
                "Error: The database object for MechanicsLinearFEOperator is NULL" );

    auto mechanicsOpParams =
        std::make_shared<MechanicsLinearFEOperatorParameters>( mechanicsLinFEOp_db );
    mechanicsOpParams->d_materialModel =
        std::dynamic_pointer_cast<MechanicsMaterialModel>( elementPhysicsModel );
    mechanicsOpParams->d_elemOp   = mechanicsLinElem;
    mechanicsOpParams->d_Mesh     = meshAdapter;
    mechanicsOpParams->d_inDofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );
    mechanicsOpParams->d_outDofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );

    auto mechanicsOp = std::make_shared<MechanicsLinearFEOperator>( mechanicsOpParams );

    return mechanicsOp;
}


Operator::shared_ptr OperatorBuilder::createNonlinearMechanicsOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{

    // first create a MechanicsMaterialModel
    if ( elementPhysicsModel.get() == nullptr ) {
        AMP_INSIST( input_db->keyExists( "MechanicsMaterialModel" ),
                    "Key ''MechanicsMaterialModel'' is missing!" );

        std::shared_ptr<AMP::Database> transportModel_db =
            input_db->getDatabase( "MechanicsMaterialModel" );
        elementPhysicsModel =
            ElementPhysicsModelFactory::createElementPhysicsModel( transportModel_db );
    }

    AMP_INSIST( elementPhysicsModel, "NULL material model" );

    // next create a ElementOperation object
    AMP_INSIST( input_db->keyExists( "MechanicsElement" ), "Key ''MechanicsElement'' is missing!" );
    std::shared_ptr<ElementOperation> mechanicsElem =
        ElementOperationFactory::createElementOperation(
            input_db->getDatabase( "MechanicsElement" ) );

    // now create the nonlinear mechanics operator
    std::shared_ptr<AMP::Database> mechanicsFEOp_db;
    if ( input_db->getString( "name" ) == "MechanicsNonlinearFEOperator" ) {
        mechanicsFEOp_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    AMP_INSIST( mechanicsFEOp_db,
                "Error: The database object for MechanicsNonlinearFEOperator is NULL" );

    auto mechanicsOpParams =
        std::make_shared<MechanicsNonlinearFEOperatorParameters>( mechanicsFEOp_db );
    mechanicsOpParams->d_materialModel =
        std::dynamic_pointer_cast<MechanicsMaterialModel>( elementPhysicsModel );
    mechanicsOpParams->d_elemOp = mechanicsElem;
    mechanicsOpParams->d_Mesh   = meshAdapter;
    mechanicsOpParams->d_dofMap[Mechanics::DISPLACEMENT] =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );
    mechanicsOpParams->d_dofMap[Mechanics::TEMPERATURE] =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    mechanicsOpParams->d_dofMap[Mechanics::BURNUP] = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    mechanicsOpParams->d_dofMap[Mechanics::OXYGEN_CONCENTRATION] =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    mechanicsOpParams->d_dofMap[Mechanics::LHGR] = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );

    auto mechanicsOp = std::make_shared<MechanicsNonlinearFEOperator>( mechanicsOpParams );

    return mechanicsOp;
}

Operator::shared_ptr OperatorBuilder::createLinearNavierStokesLSWFOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{

    if ( elementPhysicsModel.get() == nullptr ) {
        AMP_INSIST( input_db->keyExists( "FlowTransportModel" ),
                    "Key ''FlowTransportModel'' is missing!" );

        auto transportModel_db = input_db->getDatabase( "FlowTransportModel" );
        elementPhysicsModel =
            ElementPhysicsModelFactory::createElementPhysicsModel( transportModel_db );
    }

    AMP_INSIST( elementPhysicsModel, "NULL transport model" );

    // next create a ElementOperation object
    AMP_INSIST( input_db->keyExists( "FlowElement" ), "Key ''FlowElement'' is missing!" );
    auto flowLinElem =
        ElementOperationFactory::createElementOperation( input_db->getDatabase( "FlowElement" ) );

    // now create the linear flow operator
    std::shared_ptr<AMP::Database> flowLinFEOp_db;
    if ( input_db->getString( "name" ) == "NavierStokesLSWFLinearFEOperator" ) {
        flowLinFEOp_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    AMP_INSIST( flowLinFEOp_db, "Error: The database object for FlowLinearFEOperator is NULL" );

    auto flowOpParams = std::make_shared<NavierStokesLinearFEOperatorParameters>( flowLinFEOp_db );
    flowOpParams->d_transportModel =
        std::dynamic_pointer_cast<FlowTransportModel>( elementPhysicsModel );
    flowOpParams->d_elemOp   = flowLinElem;
    flowOpParams->d_Mesh     = meshAdapter;
    flowOpParams->d_inDofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 10, true );
    flowOpParams->d_outDofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 10, true );

    auto flowOp = std::make_shared<NavierStokesLSWFLinearFEOperator>( flowOpParams );

    return flowOp;
}


Operator::shared_ptr OperatorBuilder::createNonlinearNavierStokesLSWFOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{

    if ( elementPhysicsModel.get() == nullptr ) {
        AMP_INSIST( input_db->keyExists( "FlowTransportModel" ),
                    "Key ''FlowTransportModel'' is missing!" );

        auto transportModel_db = input_db->getDatabase( "FlowTransportModel" );
        elementPhysicsModel =
            ElementPhysicsModelFactory::createElementPhysicsModel( transportModel_db );
    }

    AMP_INSIST( elementPhysicsModel, "NULL material model" );

    // next create a ElementOperation object
    AMP_INSIST( input_db->keyExists( "FlowElement" ), "Key ''FlowElement'' is missing!" );
    auto flowElem =
        ElementOperationFactory::createElementOperation( input_db->getDatabase( "FlowElement" ) );

    // now create the nonlinear mechanics operator
    std::shared_ptr<AMP::Database> flowFEOp_db;
    if ( input_db->getString( "name" ) == "NavierStokesLSWFFEOperator" ) {
        flowFEOp_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    AMP_INSIST( flowFEOp_db, "Error: The database object for FlowNonlinearFEOperator is NULL" );

    auto flowOpParams = std::make_shared<NavierStokesLSWFFEOperatorParameters>( flowFEOp_db );
    flowOpParams->d_transportModel =
        std::dynamic_pointer_cast<FlowTransportModel>( elementPhysicsModel );
    flowOpParams->d_elemOp = flowElem;
    flowOpParams->d_Mesh   = meshAdapter;
    flowOpParams->d_dofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 10, true );
    auto flowOp = std::make_shared<NavierStokesLSWFFEOperator>( flowOpParams );

    return flowOp;
}


Operator::shared_ptr OperatorBuilder::createMassLinearFEOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
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
    std::shared_ptr<AMP::Database> densityLinFEOp_db;
    if ( input_db->getString( "name" ) == "MassLinearFEOperator" ) {
        densityLinFEOp_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    AMP_INSIST( densityLinFEOp_db, "Error: The database object for MassLinearFEOperator is NULL" );

    auto densityOpParams = std::make_shared<MassLinearFEOperatorParameters>( densityLinFEOp_db );
    densityOpParams->d_densityModel = densityModel;
    densityOpParams->d_elemOp       = densityLinElem;
    densityOpParams->d_Mesh         = meshAdapter;
    densityOpParams->d_inDofMap     = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    densityOpParams->d_outDofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    auto densityOp = std::make_shared<MassLinearFEOperator>( densityOpParams );

    return densityOp;
}


Operator::shared_ptr OperatorBuilder::createLinearBVPOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    std::string operatorName,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> elementPhysicsModel,
    std::shared_ptr<ElementPhysicsModelFactory> localModelFactory )
{
    std::shared_ptr<Operator> retOperator;
    AMP_INSIST( input_db, "NULL database object passed" );

    auto operator_db = input_db->getDatabase( operatorName );
    AMP_INSIST( operator_db, "NULL database object passed" );

    // create the volume operator
    std::string volumeOperatorName = operator_db->getString( "VolumeOperator" );
    // if this flag is true the same local physics model will be used for both boundary and volume
    // operators
    bool useSameLocalModelForVolumeAndBoundaryOperators = operator_db->getWithDefault<bool>(
        "useSameLocalModelForVolumeAndBoundaryOperators", false );

    auto volumeOperator = OperatorBuilder::createOperator(
        meshAdapter, volumeOperatorName, input_db, elementPhysicsModel, localModelFactory );

    auto volumeLinearOp = std::dynamic_pointer_cast<LinearOperator>( volumeOperator );
    AMP_INSIST(
        volumeLinearOp,
        "Error: unable to create linear operator in OperatorBuilder::createLinearBVPOperator" );

    // create the boundary operator
    std::string boundaryOperatorName = operator_db->getString( "BoundaryOperator" );
    auto boundaryOperator_db         = input_db->getDatabase( boundaryOperatorName );
    AMP_INSIST( boundaryOperator_db, "NULL database object passed for boundary operator" );

    boundaryOperator_db->putScalar(
        "isAttachedToVolumeOperator", true, Units(), Database::Check::Overwrite );

    std::shared_ptr<ElementPhysicsModel> boundaryLocalModel;

    if ( useSameLocalModelForVolumeAndBoundaryOperators ) {
        boundaryLocalModel = elementPhysicsModel;
    }

    std::shared_ptr<BoundaryOperator> boundaryOperator =
        OperatorBuilder::createBoundaryOperator( meshAdapter,
                                                 boundaryOperatorName,
                                                 input_db,
                                                 volumeLinearOp,
                                                 boundaryLocalModel,
                                                 localModelFactory );

    auto bvpOperatorParams                = std::make_shared<BVPOperatorParameters>( input_db );
    bvpOperatorParams->d_volumeOperator   = volumeOperator;
    bvpOperatorParams->d_boundaryOperator = boundaryOperator;

    retOperator.reset( new LinearBVPOperator( bvpOperatorParams ) );

    return retOperator;
}


Operator::shared_ptr OperatorBuilder::createNonlinearBVPOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    std::string operatorName,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> elementPhysicsModel,
    std::shared_ptr<ElementPhysicsModelFactory> localModelFactory )
{
    std::shared_ptr<Operator> retOperator;
    AMP_INSIST( input_db, "NULL database object passed" );

    auto operator_db = input_db->getDatabase( operatorName );
    AMP_INSIST( operator_db, "NULL database object passed" );

    // create the volume operator
    std::string volumeOperatorName = operator_db->getString( "VolumeOperator" );
    // if this flag is true the same local physics model will be used for both boundary and volume
    // operators
    bool useSameLocalModelForVolumeAndBoundaryOperators = operator_db->getWithDefault<bool>(
        "useSameLocalModelForVolumeAndBoundaryOperators", false );

    auto volumeOperator = OperatorBuilder::createOperator(
        meshAdapter, volumeOperatorName, input_db, elementPhysicsModel, localModelFactory );
    AMP_INSIST( volumeOperator,
                "Error: unable to create nonlinear operator in "
                "OperatorBuilder::createNonlinearBVPOperator" );

    // create the boundary operator
    std::string boundaryOperatorName = operator_db->getString( "BoundaryOperator" );
    auto boundaryOperator_db         = input_db->getDatabase( boundaryOperatorName );
    AMP_INSIST( boundaryOperator_db, "NULL database object passed for boundary operator" );

    boundaryOperator_db->putScalar(
        "isAttachedToVolumeOperator", true, Units(), Database::Check::Overwrite );

    std::shared_ptr<ElementPhysicsModel> boundaryLocalModel;

    if ( useSameLocalModelForVolumeAndBoundaryOperators ) {
        boundaryLocalModel = elementPhysicsModel;
    }

    auto boundaryOperator = OperatorBuilder::createBoundaryOperator( meshAdapter,
                                                                     boundaryOperatorName,
                                                                     input_db,
                                                                     volumeOperator,
                                                                     boundaryLocalModel,
                                                                     localModelFactory );

    auto bvpOperatorParams                = std::make_shared<BVPOperatorParameters>( input_db );
    bvpOperatorParams->d_volumeOperator   = volumeOperator;
    bvpOperatorParams->d_boundaryOperator = boundaryOperator;

    retOperator.reset( new NonlinearBVPOperator( bvpOperatorParams ) );

    return retOperator;
}


Operator::shared_ptr
OperatorBuilder::createFlowFrapconJacobian( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                            std::shared_ptr<AMP::Database> input_db )
{

    // now create the flow frapcon operator
    std::shared_ptr<AMP::Database> flowOp_db;
    if ( input_db->getString( "name" ) == "FlowFrapconJacobian" ) {
        flowOp_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    AMP_INSIST( flowOp_db, "Error: The database object for FlowFrapconJacobian is NULL" );

    auto flowOpParams    = std::make_shared<FlowFrapconJacobianParameters>( flowOp_db );
    flowOpParams->d_Mesh = meshAdapter;
    std::shared_ptr<FlowFrapconJacobian> flowOp( new FlowFrapconJacobian( flowOpParams ) );

    return flowOp;
}

std::shared_ptr<BoundaryOperator> OperatorBuilder::createBoundaryOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    std::string boundaryOperatorName,
    std::shared_ptr<AMP::Database> input_db,
    Operator::shared_ptr volumeOperator,
    std::shared_ptr<ElementPhysicsModel> elementPhysicsModel,
    std::shared_ptr<ElementPhysicsModelFactory> localModelFactory )
{
    std::shared_ptr<BoundaryOperator> retOperator;
    AMP_INSIST( input_db, "NULL database object passed" );

    auto operator_db = input_db->getDatabase( boundaryOperatorName );
    AMP_INSIST( operator_db,
                "Error: OperatorBuilder::createBoundaryOperator(): "
                "database object with given name not in database" );

    // we create the element physics model if a database entry exists
    // and the incoming element physics model pointer is NULL
    //  if( (elementPhysicsModel.get()==NULL) && (operator_db->keyExists("LocalModel" ) ) )
    //  The above Condition assumes all of the Operators inside column boundary use same Physics
    //  Model - SA
    if ( ( operator_db->keyExists( "LocalModel" ) ) ) {
        // extract the name of the local model from the operator database
        std::string localModelName = operator_db->getString( "LocalModel" );
        // check whether a database exists in the global database
        // (NOTE: not the operator database) with the given name
        AMP_INSIST( input_db->keyExists( localModelName ),
                    "Error:: OperatorBuilder::createOperator(): No local model "
                    "database entry with given name exists in input database" );

        std::shared_ptr<AMP::Database> localModel_db = input_db->getDatabase( localModelName );
        AMP_INSIST( localModel_db,
                    "Error:: OperatorBuilder::createOperator(): No local model database "
                    "entry with given name exists in input database" );

        // if a non-NULL factory is being supplied through the argument list
        // use it, else call the AMP ElementPhysicsModelFactory interface
        if ( localModelFactory ) {
            elementPhysicsModel = localModelFactory->createElementPhysicsModel( localModel_db );
        } else {
            elementPhysicsModel =
                ElementPhysicsModelFactory::createElementPhysicsModel( localModel_db );
        }

        AMP_INSIST( elementPhysicsModel,
                    "Error:: OperatorBuilder::createOperator(): local model creation failed" );
    }

    auto boundaryType = operator_db->getString( "name" );

    if ( boundaryType == "DirichletMatrixCorrection" ) {
        // in this case the volume operator has to be a linear operator
        retOperator = createDirichletMatrixCorrection(
            meshAdapter, operator_db, volumeOperator, elementPhysicsModel );
    } else if ( boundaryType == "MassMatrixCorrection" ) {
        retOperator = createMassMatrixCorrection(
            meshAdapter, operator_db, volumeOperator, elementPhysicsModel );
    } else if ( boundaryType == "RobinMatrixCorrection" ) {
        // in this case the volume operator has to be a linear operator
        retOperator = createRobinMatrixCorrection(
            meshAdapter, operator_db, volumeOperator, elementPhysicsModel );
    } else if ( boundaryType == "RobinVectorCorrection" ) {
        retOperator = createRobinVectorCorrection(
            meshAdapter, operator_db, volumeOperator, elementPhysicsModel );
    } else if ( boundaryType == "NeumannVectorCorrection" ) {
        retOperator = createNeumannVectorCorrection(
            meshAdapter, operator_db, volumeOperator, elementPhysicsModel );
    } else if ( boundaryType == "DirichletVectorCorrection" ) {
        // in this case the volume operator has to be a nonlinear operator
        retOperator = createDirichletVectorCorrection(
            meshAdapter, operator_db, volumeOperator, elementPhysicsModel );
    } else if ( boundaryType == "PressureBoundaryOperator" ) {
        retOperator =
            createPressureBoundaryOperator( meshAdapter, operator_db, elementPhysicsModel );
    } else if ( boundaryType == "ColumnBoundaryOperator" ) {
        // note that the global input database is passed here instead of the operator
        // database
        retOperator = createColumnBoundaryOperator( meshAdapter,
                                                    boundaryOperatorName,
                                                    input_db,
                                                    volumeOperator,
                                                    elementPhysicsModel,
                                                    localModelFactory );
    }

    return retOperator;
}


std::shared_ptr<BoundaryOperator> OperatorBuilder::createColumnBoundaryOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    std::string boundaryOperatorName,
    std::shared_ptr<AMP::Database> input_db,
    Operator::shared_ptr volumeOperator,
    std::shared_ptr<ElementPhysicsModel> elementPhysicsModel,
    std::shared_ptr<ElementPhysicsModelFactory> localModelFactory )
{

    AMP_INSIST( input_db, "NULL database object passed" );

    auto operator_db = input_db->getDatabase( boundaryOperatorName );
    AMP_INSIST( operator_db,
                "Error: OperatorBuilder::createBoundaryOperator(): "
                "database object with given name not in database" );

    int numberOfBoundaryOperators =
        operator_db->getWithDefault<int>( "numberOfBoundaryOperators", 1 );

    auto boundaryOps = operator_db->getVector<std::string>( "boundaryOperators" );
    AMP_ASSERT( numberOfBoundaryOperators == (int) boundaryOps.size() );

    auto params    = std::make_shared<OperatorParameters>( operator_db );
    params->d_Mesh = meshAdapter;

    auto columnBoundaryOperator = std::make_shared<ColumnBoundaryOperator>( params );

    for ( int i = 0; i < numberOfBoundaryOperators; i++ ) {
        auto bcOperator = OperatorBuilder::createBoundaryOperator( meshAdapter,
                                                                   boundaryOps[i],
                                                                   input_db,
                                                                   volumeOperator,
                                                                   elementPhysicsModel,
                                                                   localModelFactory );
        AMP_ASSERT( bcOperator );
        columnBoundaryOperator->append( bcOperator );
    }

    return columnBoundaryOperator;
}


std::shared_ptr<BoundaryOperator>
OperatorBuilder::createDirichletMatrixCorrection( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                                  std::shared_ptr<AMP::Database> input_db,
                                                  Operator::shared_ptr volumeOperator,
                                                  std::shared_ptr<ElementPhysicsModel> )
{
    std::shared_ptr<BoundaryOperator> retOperator;
    auto linearOperator = std::dynamic_pointer_cast<LinearOperator>( volumeOperator );
    auto matrixCorrectionParameters =
        std::make_shared<DirichletMatrixCorrectionParameters>( input_db );
    matrixCorrectionParameters->d_variable    = linearOperator->getOutputVariable();
    matrixCorrectionParameters->d_inputMatrix = linearOperator->getMatrix();
    matrixCorrectionParameters->d_Mesh        = meshAdapter;

    retOperator.reset( new DirichletMatrixCorrection( matrixCorrectionParameters ) );

    return retOperator;
}


std::shared_ptr<BoundaryOperator>
OperatorBuilder::createMassMatrixCorrection( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                             std::shared_ptr<AMP::Database> input_db,
                                             Operator::shared_ptr volumeOperator,
                                             std::shared_ptr<ElementPhysicsModel> )
{
    std::shared_ptr<BoundaryOperator> retOperator;
    auto linearOperator = std::dynamic_pointer_cast<LinearOperator>( volumeOperator );
    auto matrixCorrectionParameters =
        std::make_shared<DirichletMatrixCorrectionParameters>( input_db );
    matrixCorrectionParameters->d_variable    = linearOperator->getOutputVariable();
    matrixCorrectionParameters->d_inputMatrix = linearOperator->getMatrix();
    matrixCorrectionParameters->d_Mesh        = meshAdapter;

    retOperator.reset( new MassMatrixCorrection( matrixCorrectionParameters ) );

    return retOperator;
}


std::shared_ptr<BoundaryOperator> OperatorBuilder::createRobinMatrixCorrection(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    Operator::shared_ptr volumeOperator,
    std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{
    std::shared_ptr<BoundaryOperator> retOperator;
    auto linearOperator             = std::dynamic_pointer_cast<LinearOperator>( volumeOperator );
    auto matrixCorrectionParameters = std::make_shared<RobinMatrixCorrectionParameters>( input_db );
    matrixCorrectionParameters->d_variable    = linearOperator->getOutputVariable();
    matrixCorrectionParameters->d_inputMatrix = linearOperator->getMatrix();
    matrixCorrectionParameters->d_Mesh        = meshAdapter;

    if ( elementPhysicsModel ) {

        std::shared_ptr<RobinPhysicsModel> robinPhysicsModel;
        robinPhysicsModel = std::dynamic_pointer_cast<RobinPhysicsModel>( elementPhysicsModel );
        matrixCorrectionParameters->d_robinPhysicsModel = robinPhysicsModel;
    }

    retOperator.reset( new RobinMatrixCorrection( matrixCorrectionParameters ) );

    return retOperator;
}


std::shared_ptr<BoundaryOperator> OperatorBuilder::createRobinVectorCorrection(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    Operator::shared_ptr volumeOperator,
    std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{
    std::shared_ptr<BoundaryOperator> retOperator;
    auto vectorCorrectionParameters =
        std::make_shared<NeumannVectorCorrectionParameters>( input_db );
    vectorCorrectionParameters->d_variable = volumeOperator->getOutputVariable();
    vectorCorrectionParameters->d_Mesh     = meshAdapter;

    if ( elementPhysicsModel ) {
        std::shared_ptr<RobinPhysicsModel> robinPhysicsModel;
        robinPhysicsModel = std::dynamic_pointer_cast<RobinPhysicsModel>( elementPhysicsModel );
        vectorCorrectionParameters->d_robinPhysicsModel = robinPhysicsModel;
    }

    retOperator.reset( new RobinVectorCorrection( vectorCorrectionParameters ) );

    return retOperator;
}


std::shared_ptr<BoundaryOperator> OperatorBuilder::createNeumannVectorCorrection(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    Operator::shared_ptr volumeOperator,
    std::shared_ptr<ElementPhysicsModel> elementPhysicsModel )
{
    std::shared_ptr<BoundaryOperator> retOperator;
    auto vectorCorrectionParameters =
        std::make_shared<NeumannVectorCorrectionParameters>( input_db );
    vectorCorrectionParameters->d_variable = volumeOperator->getOutputVariable();
    vectorCorrectionParameters->d_Mesh     = meshAdapter;

    if ( elementPhysicsModel && input_db->isDatabase( "RobinPhysicsModel" ) ) {
        std::shared_ptr<RobinPhysicsModel> robinPhysicsModel;
        robinPhysicsModel = std::dynamic_pointer_cast<RobinPhysicsModel>( elementPhysicsModel );
        vectorCorrectionParameters->d_robinPhysicsModel = robinPhysicsModel;
    }

    retOperator.reset( new NeumannVectorCorrection( vectorCorrectionParameters ) );

    return retOperator;
}


std::shared_ptr<BoundaryOperator>
OperatorBuilder::createDirichletVectorCorrection( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                                  std::shared_ptr<AMP::Database> input_db,
                                                  Operator::shared_ptr volumeOperator,
                                                  std::shared_ptr<ElementPhysicsModel> )
{
    std::shared_ptr<BoundaryOperator> retOperator;
    auto vectorCorrectionParameters =
        std::make_shared<DirichletVectorCorrectionParameters>( input_db );
    vectorCorrectionParameters->d_variable = volumeOperator->getOutputVariable();
    vectorCorrectionParameters->d_Mesh     = meshAdapter;

    retOperator.reset( new DirichletVectorCorrection( vectorCorrectionParameters ) );

    return retOperator;
}


std::shared_ptr<BoundaryOperator>
OperatorBuilder::createDirichletVectorCorrection( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                                  std::shared_ptr<AMP::Database> input_db,
                                                  std::shared_ptr<ElementPhysicsModel> )
{
    std::shared_ptr<BoundaryOperator> retOperator;
    auto vectorCorrectionParameters =
        std::make_shared<DirichletVectorCorrectionParameters>( input_db );
    vectorCorrectionParameters->d_Mesh = meshAdapter;

    retOperator.reset( new DirichletVectorCorrection( vectorCorrectionParameters ) );

    return retOperator;
}


std::shared_ptr<BoundaryOperator>
OperatorBuilder::createPressureBoundaryOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                                 std::shared_ptr<AMP::Database> input_db,
                                                 std::shared_ptr<ElementPhysicsModel> )
{
    std::shared_ptr<BoundaryOperator> retOperator;
    auto params    = std::make_shared<OperatorParameters>( input_db );
    params->d_Mesh = meshAdapter;

    retOperator.reset( new PressureBoundaryOperator( params ) );

    return retOperator;
}


std::shared_ptr<Operator> OperatorBuilder::createOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter1,
                                                           AMP::Mesh::Mesh::shared_ptr meshAdapter2,
                                                           const AMP::AMP_MPI &comm,
                                                           std::shared_ptr<AMP::Database> input_db )
{
    std::shared_ptr<Operator> retOperator;

    std::string name = input_db->getString( "name" );
    if ( name == "MapSurface" ) {
        auto mapOperatorParameters       = std::make_shared<MapOperatorParameters>( input_db );
        mapOperatorParameters->d_Mesh    = meshAdapter1;
        mapOperatorParameters->d_MapMesh = meshAdapter2;
        mapOperatorParameters->d_MapComm = comm;
        retOperator.reset( new MapSurface( mapOperatorParameters ) );
    }
    return retOperator;
}


#else

std::shared_ptr<Operator> OperatorBuilder::createIdentityOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                                                   std::shared_ptr<AMP::Database> )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator>
    OperatorBuilder::createLinearMechanicsOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                                    std::shared_ptr<AMP::Database>,
                                                    std::shared_ptr<ElementPhysicsModel> )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator>
    OperatorBuilder::createNonlinearMechanicsOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                                       std::shared_ptr<AMP::Database>,
                                                       std::shared_ptr<ElementPhysicsModel> )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator>
    OperatorBuilder::createLinearDiffusionOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                                    std::shared_ptr<AMP::Database>,
                                                    std::shared_ptr<ElementPhysicsModel> )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator>
    OperatorBuilder::createNonlinearDiffusionOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                                       std::shared_ptr<AMP::Database>,
                                                       std::shared_ptr<ElementPhysicsModel> )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator>
    OperatorBuilder::createLinearNavierStokesLSWFOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                                           std::shared_ptr<AMP::Database>,
                                                           std::shared_ptr<ElementPhysicsModel> )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator>
    OperatorBuilder::createNonlinearNavierStokesLSWFOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                                              std::shared_ptr<AMP::Database>,
                                                              std::shared_ptr<ElementPhysicsModel> )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator> OperatorBuilder::createNonlinearFickSoretOperator(
    std::shared_ptr<AMP::Mesh::Mesh>,
    std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char>>,
    std::shared_ptr<AMP::Database>,
    std::shared_ptr<ElementPhysicsModel>,
    std::shared_ptr<ElementPhysicsModelFactory> )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator>
    OperatorBuilder::createFlowFrapconOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                                std::shared_ptr<AMP::Database> )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator>
    OperatorBuilder::createFlowFrapconJacobian( std::shared_ptr<AMP::Mesh::Mesh>,
                                                std::shared_ptr<AMP::Database> )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator>
    OperatorBuilder::createSubchannelTwoEqLinearOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                                          std::shared_ptr<AMP::Database>,
                                                          std::shared_ptr<ElementPhysicsModel> )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator>
    OperatorBuilder::createSubchannelTwoEqNonlinearOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                                             std::shared_ptr<AMP::Database>,
                                                             std::shared_ptr<ElementPhysicsModel> )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator>
    OperatorBuilder::createSubchannelFourEqNonlinearOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                                              std::shared_ptr<AMP::Database>,
                                                              std::shared_ptr<ElementPhysicsModel> )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator>
    OperatorBuilder::createNeutronicsRhsOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                                  std::shared_ptr<AMP::Database> )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator>
    OperatorBuilder::createMassLinearFEOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                                 std::shared_ptr<AMP::Database>,
                                                 std::shared_ptr<ElementPhysicsModel> )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator>
    OperatorBuilder::createVolumeIntegralOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                                   std::shared_ptr<AMP::Database>,
                                                   std::shared_ptr<ElementPhysicsModel> )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator> OperatorBuilder::createLinearBVPOperator(
    std::shared_ptr<AMP::Mesh::Mesh>,
    std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char>>,
    std::shared_ptr<AMP::Database>,
    std::shared_ptr<ElementPhysicsModel>,
    std::shared_ptr<ElementPhysicsModelFactory> )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator> OperatorBuilder::createNonlinearBVPOperator(
    std::shared_ptr<AMP::Mesh::Mesh>,
    std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char>>,
    std::shared_ptr<AMP::Database>,
    std::shared_ptr<ElementPhysicsModel>,
    std::shared_ptr<ElementPhysicsModelFactory> )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<BoundaryOperator>
    OperatorBuilder::createDirichletVectorCorrection( std::shared_ptr<AMP::Mesh::Mesh>,
                                                      std::shared_ptr<AMP::Database>,
                                                      std::shared_ptr<ElementPhysicsModel> )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<BoundaryOperator>
    OperatorBuilder::createPressureBoundaryOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                                     std::shared_ptr<AMP::Database>,
                                                     std::shared_ptr<ElementPhysicsModel> )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}


#endif

} // namespace AMP::Operator
