#include "OperatorBuilder.h"
#include "utils/Utilities.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/Variable.h"
#include "vectors/VectorBuilder.h"

#include "ampmesh/StructuredMeshHelper.h"

#include <string>


#include "operators/IdentityOperator.h"

#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"

#ifdef USE_EXT_LIBMESH
#include "discretization/structuredFaceDOFManager.h"
#include "operators/ElementOperationFactory.h"
#include "operators/NeutronicsRhs.h"
#include "operators/ParameterFactory.h"
#include "operators/boundary/ColumnBoundaryOperator.h"
#include "operators/boundary/DirichletMatrixCorrection.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/boundary/MassMatrixCorrection.h"
#include "operators/boundary/libmesh/NeumannVectorCorrection.h"
#include "operators/boundary/libmesh/PressureBoundaryOperator.h"
#include "operators/boundary/libmesh/RobinMatrixCorrection.h"
#include "operators/boundary/libmesh/RobinVectorCorrection.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "operators/diffusion/FickSoretNonlinearFEOperator.h"
#include "operators/flow/NavierStokesLSWFFEOperator.h"
#include "operators/flow/NavierStokesLSWFLinearFEOperator.h"
#include "operators/libmesh/MassLinearFEOperator.h"
#include "operators/libmesh/VolumeIntegralOperator.h"
#include "operators/map/MapSurface.h"
#include "operators/mechanics/MechanicsConstants.h"
#include "operators/mechanics/MechanicsLinearFEOperator.h"
#include "operators/mechanics/MechanicsNonlinearFEOperator.h"
#include "operators/subchannel/FlowFrapconJacobian.h"
#include "operators/subchannel/FlowFrapconOperator.h"
#include "operators/subchannel/SubchannelFourEqNonlinearOperator.h"
#include "operators/subchannel/SubchannelTwoEqLinearOperator.h"
#include "operators/subchannel/SubchannelTwoEqNonlinearOperator.h"
#endif


#define resetOperation( NAME )                                            \
    do {                                                                  \
        if ( name == #NAME ) {                                            \
            AMP::shared_ptr<NAME##Parameters> params =                    \
                AMP::dynamic_pointer_cast<NAME##Parameters>( in_params ); \
            AMP_ASSERT( params.get() == in_params.get() );                \
            retOperator.reset( new NAME( params ) );                      \
        }                                                                 \
    } while ( 0 )


namespace AMP {
namespace Operator {


typedef OperatorParameters IdentityOperatorParameters;


AMP::shared_ptr<Operator>
OperatorBuilder::createOperator( AMP::shared_ptr<OperatorParameters> in_params )
{
    AMP::shared_ptr<Operator> retOperator;

    AMP_INSIST( in_params.get() != nullptr,
                "ERROR: OperatorBuilder::createOperator has NULL input" );
    AMP_INSIST( in_params->d_db.get() != nullptr,
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
        AMP::shared_ptr<TractionBoundaryOperatorParameters> params =
            AMP::dynamic_pointer_cast<TractionBoundaryOperatorParameters>( in_params );
        AMP_ASSERT( params.get() == in_params.get() );
        retOperator.reset( new PressureBoundaryOperator( params ) );
    }
#endif
    if ( ( name == "LinearBVPOperator" ) || ( name == "NonlinearBVPOperator" ) ) {
        AMP::shared_ptr<BVPOperatorParameters> bvpOperatorParameters =
            AMP::dynamic_pointer_cast<BVPOperatorParameters>( in_params );

        AMP_INSIST( bvpOperatorParameters.get() != nullptr,
                    "ERROR: NULL BVPOperatorParameters passed" );
        AMP_INSIST( bvpOperatorParameters->d_volumeOperatorParams.get() != nullptr,
                    "ERROR: BVPOperatorParameters has NULL volumeOperatorParams pointer" );
        AMP_INSIST( bvpOperatorParameters->d_boundaryOperatorParams.get() != nullptr,
                    "ERROR: BVPOperatorParameters has NULL boundaryOperatorParams pointer" );

        bvpOperatorParameters->d_volumeOperator =
            OperatorBuilder::createOperator( bvpOperatorParameters->d_volumeOperatorParams );
        bvpOperatorParameters->d_boundaryOperator = AMP::dynamic_pointer_cast<BoundaryOperator>(
            OperatorBuilder::createOperator( bvpOperatorParameters->d_boundaryOperatorParams ) );

        if ( name == "LinearBVPOperator" )
            retOperator.reset( new LinearBVPOperator(
                AMP::dynamic_pointer_cast<BVPOperatorParameters>( in_params ) ) );
        else
            retOperator.reset( new NonlinearBVPOperator(
                AMP::dynamic_pointer_cast<BVPOperatorParameters>( in_params ) ) );
    }

    return retOperator;
}


AMP::shared_ptr<Operator>
OperatorBuilder::createOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                 std::string operatorName,
                                 AMP::shared_ptr<AMP::Database>
                                     tmp_input_db,
                                 AMP::shared_ptr<AMP::Operator::ElementPhysicsModel>
                                     elementPhysicsModel,
                                 AMP::shared_ptr<AMP::Operator::ElementPhysicsModelFactory>
                                     localModelFactory )
{

    AMP::shared_ptr<Operator> retOperator;

    AMP::shared_ptr<AMP::Database> input_db =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>( tmp_input_db );

    AMP::shared_ptr<AMP::Database> operator_db =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>( input_db->getDatabase( operatorName ) );

    AMP_INSIST( operator_db.get() != nullptr,
                "Error:: OperatorBuilder::createOperator(): No operator database entry with "
                "given name exists in input database" );

    // we create the element physics model if a database entry exists
    // and the incoming element physics model pointer is NULL
    if ( ( elementPhysicsModel.get() == nullptr ) && ( operator_db->keyExists( "LocalModel" ) ) ) {
        // extract the name of the local model from the operator database
        std::string localModelName = operator_db->getString( "LocalModel" );
        // check whether a database exists in the global database
        // (NOTE: not the operator database) with the given name
        AMP_INSIST( input_db->keyExists( localModelName ),
                    "Error:: OperatorBuilder::createOperator(): No local model "
                    "database entry with given name exists in input database" );

        AMP::shared_ptr<AMP::Database> localModel_db = input_db->getDatabase( localModelName );
        AMP_INSIST( localModel_db.get() != nullptr,
                    "Error:: OperatorBuilder::createOperator(): No local model database "
                    "entry with given name exists in input databaseot" );

        // If a non-NULL factory is being supplied through the argument list
        // use it, else call the AMP ElementPhysicsModelFactory interface
        if ( localModelFactory.get() != nullptr ) {
            elementPhysicsModel = localModelFactory->createElementPhysicsModel( localModel_db );
        } else {
            elementPhysicsModel =
                ElementPhysicsModelFactory::createElementPhysicsModel( localModel_db );
        }

        AMP_INSIST( elementPhysicsModel.get() != nullptr,
                    "Error:: OperatorBuilder::createOperator(): local model creation failed" );
    }

    std::string operatorType = operator_db->getString( "name" );

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
    }

    return retOperator;
}


#ifdef USE_EXT_LIBMESH


// Create the identity operator
AMP::Operator::Operator::shared_ptr OperatorBuilder::createIdentityOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter, AMP::shared_ptr<AMP::Database> input_db )
{
    AMP_INSIST( input_db.get() != nullptr,
                "Error: The database object for SubchannelTwoEqLinearOperator is NULL" );

    AMP::shared_ptr<AMP::Operator::OperatorParameters> params(
        new AMP::Operator::OperatorParameters( input_db ) );
    params->d_Mesh = meshAdapter;
    AMP::shared_ptr<AMP::Operator::IdentityOperator> subchannelOp(
        new AMP::Operator::IdentityOperator( params ) );

    return subchannelOp;
}


// Create the FlowFrapconOperator
AMP::Operator::Operator::shared_ptr OperatorBuilder::createFlowFrapconOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter, AMP::shared_ptr<AMP::Database> input_db )
{

    // now create the flow frapcon operator
    AMP::shared_ptr<AMP::Database> flowOp_db;
    if ( input_db->getString( "name" ) == "FlowFrapconOperator" ) {
        flowOp_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    AMP_INSIST( flowOp_db.get() != nullptr,
                "Error: The database object for FlowFrapconOperator is NULL" );

    AMP::shared_ptr<AMP::Operator::FlowFrapconOperatorParameters> flowOpParams(
        new AMP::Operator::FlowFrapconOperatorParameters( flowOp_db ) );
    flowOpParams->d_Mesh = meshAdapter;
    AMP::shared_ptr<AMP::Operator::FlowFrapconOperator> flowOp(
        new AMP::Operator::FlowFrapconOperator( flowOpParams ) );

    return flowOp;
}

AMP::Operator::Operator::shared_ptr OperatorBuilder::createSubchannelTwoEqLinearOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    AMP::shared_ptr<AMP::Database>
        input_db,
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel>
        elementPhysicsModel )
{
    // first create a SubchannelPhysicsModel
    AMP::shared_ptr<AMP::Operator::SubchannelPhysicsModel> transportModel;

    if ( elementPhysicsModel.get() != nullptr ) {
        transportModel =
            AMP::dynamic_pointer_cast<AMP::Operator::SubchannelPhysicsModel>( elementPhysicsModel );
    } else {
        AMP::shared_ptr<AMP::Database> transportModel_db;
        if ( input_db->keyExists( "SubchannelPhysicsModel" ) ) {
            transportModel_db = input_db->getDatabase( "SubchannelPhysicsModel" );
        } else {
            AMP_INSIST( false, "Key ''SubchannelPhysicsModel'' is missing!" );
        }
        elementPhysicsModel =
            ElementPhysicsModelFactory::createElementPhysicsModel( transportModel_db );
        transportModel = AMP::dynamic_pointer_cast<SubchannelPhysicsModel>( elementPhysicsModel );
    }

    AMP_INSIST( transportModel.get() != nullptr, "NULL transport model" );
    // create the operator
    AMP::shared_ptr<AMP::Database> subchannel_db;
    if ( input_db->getString( "name" ) == "SubchannelTwoEqLinearOperator" ) {
        subchannel_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    AMP_INSIST( subchannel_db.get() != nullptr,
                "Error: The database object for SubchannelTwoEqLinearOperator is NULL" );

    AMP::shared_ptr<AMP::Operator::SubchannelOperatorParameters> subchannelParams(
        new AMP::Operator::SubchannelOperatorParameters( subchannel_db ) );
    subchannelParams->d_Mesh                   = meshAdapter;
    subchannelParams->d_subchannelPhysicsModel = transportModel;

    AMP::shared_ptr<AMP::Operator::SubchannelTwoEqLinearOperator> subchannelOp(
        new AMP::Operator::SubchannelTwoEqLinearOperator( subchannelParams ) );

    return subchannelOp;
}


AMP::Operator::Operator::shared_ptr OperatorBuilder::createSubchannelTwoEqNonlinearOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    AMP::shared_ptr<AMP::Database>
        input_db,
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel>
        elementPhysicsModel )
{

    // first create a SubchannelPhysicsModel
    AMP::shared_ptr<AMP::Operator::SubchannelPhysicsModel> transportModel;

    if ( elementPhysicsModel.get() != nullptr ) {
        transportModel =
            AMP::dynamic_pointer_cast<AMP::Operator::SubchannelPhysicsModel>( elementPhysicsModel );
    } else {
        AMP::shared_ptr<AMP::Database> transportModel_db;
        if ( input_db->keyExists( "SubchannelPhysicsModel" ) ) {
            transportModel_db = input_db->getDatabase( "SubchannelPhysicsModel" );
        } else {
            AMP_INSIST( false, "Key ''SubchannelPhysicsModel'' is missing!" );
        }
        elementPhysicsModel =
            ElementPhysicsModelFactory::createElementPhysicsModel( transportModel_db );
        transportModel = AMP::dynamic_pointer_cast<SubchannelPhysicsModel>( elementPhysicsModel );
    }

    AMP_INSIST( transportModel.get() != nullptr, "NULL transport model" );

    // create the operator
    AMP::shared_ptr<AMP::Database> subchannel_db;
    if ( input_db->getString( "name" ) == "SubchannelTwoEqNonlinearOperator" ) {
        subchannel_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    AMP_INSIST( subchannel_db.get() != nullptr,
                "Error: The database object for SubchannelTwoEqNonlinearOperator is NULL" );

    AMP::shared_ptr<AMP::Operator::SubchannelOperatorParameters> subchannelParams(
        new AMP::Operator::SubchannelOperatorParameters( subchannel_db ) );
    subchannelParams->d_Mesh                   = meshAdapter;
    subchannelParams->d_subchannelPhysicsModel = transportModel;
    AMP::shared_ptr<AMP::Operator::SubchannelTwoEqNonlinearOperator> subchannelOp(
        new AMP::Operator::SubchannelTwoEqNonlinearOperator( subchannelParams ) );

    return subchannelOp;
}


AMP::Operator::Operator::shared_ptr OperatorBuilder::createSubchannelFourEqNonlinearOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    AMP::shared_ptr<AMP::Database>
        input_db,
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel>
        elementPhysicsModel )
{

    // first create a SubchannelPhysicsModel
    AMP::shared_ptr<AMP::Operator::SubchannelPhysicsModel> transportModel;

    if ( elementPhysicsModel.get() != nullptr ) {
        transportModel =
            AMP::dynamic_pointer_cast<AMP::Operator::SubchannelPhysicsModel>( elementPhysicsModel );
    } else {
        AMP::shared_ptr<AMP::Database> transportModel_db;
        if ( input_db->keyExists( "SubchannelPhysicsModel" ) ) {
            transportModel_db = input_db->getDatabase( "SubchannelPhysicsModel" );
        } else {
            AMP_INSIST( false, "Key ''SubchannelPhysicsModel'' is missing!" );
        }
        elementPhysicsModel =
            ElementPhysicsModelFactory::createElementPhysicsModel( transportModel_db );
        transportModel = AMP::dynamic_pointer_cast<SubchannelPhysicsModel>( elementPhysicsModel );
    }

    AMP_INSIST( transportModel.get() != nullptr, "NULL transport model" );

    // create the operator
    AMP::shared_ptr<AMP::Database> subchannel_db;
    if ( input_db->getString( "name" ) == "SubchannelFourEqNonlinearOperator" ) {
        subchannel_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    AMP_INSIST( subchannel_db.get() != nullptr,
                "Error: The database object for SubchannelFourEqNonlinearOperator is NULL" );

    AMP::shared_ptr<AMP::Operator::SubchannelOperatorParameters> subchannelParams(
        new AMP::Operator::SubchannelOperatorParameters( subchannel_db ) );
    subchannelParams->d_Mesh                   = meshAdapter;
    subchannelParams->d_subchannelPhysicsModel = transportModel;
    AMP::shared_ptr<AMP::Operator::SubchannelFourEqNonlinearOperator> subchannelOp(
        new AMP::Operator::SubchannelFourEqNonlinearOperator( subchannelParams ) );

    return subchannelOp;
}


AMP::Operator::Operator::shared_ptr OperatorBuilder::createNeutronicsRhsOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter, AMP::shared_ptr<AMP::Database> input_db )
{

    // now create the Neutronics operator
    AMP::shared_ptr<AMP::Database> NeutronicsOp_db;
    if ( input_db->getString( "name" ) == "NeutronicsRhsOperator" ) {
        NeutronicsOp_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    AMP_INSIST( NeutronicsOp_db.get() != nullptr,
                "Error: The database object for Neutronics Source Operator is NULL" );

    AMP::shared_ptr<AMP::Operator::NeutronicsRhsParameters> neutronicsOpParams(
        new AMP::Operator::NeutronicsRhsParameters( NeutronicsOp_db ) );
    neutronicsOpParams->d_Mesh = meshAdapter;
    AMP::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOp(
        new AMP::Operator::NeutronicsRhs( neutronicsOpParams ) );

    return neutronicsOp;
}


AMP::Operator::Operator::shared_ptr
OperatorBuilder::createLinearDiffusionOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                                AMP::shared_ptr<AMP::Database>
                                                    input_db,
                                                AMP::shared_ptr<AMP::Operator::ElementPhysicsModel>
                                                    elementPhysicsModel )
{

    // first create a DiffusionTransportModel
    AMP::shared_ptr<AMP::Operator::DiffusionTransportModel> transportModel;

    if ( elementPhysicsModel.get() != nullptr ) {
        transportModel = AMP::dynamic_pointer_cast<AMP::Operator::DiffusionTransportModel>(
            elementPhysicsModel );
    } else {
        AMP::shared_ptr<AMP::Database> transportModel_db;
        if ( input_db->keyExists( "DiffusionTransportModel" ) ) {
            transportModel_db = input_db->getDatabase( "DiffusionTransportModel" );
        } else {
            AMP_INSIST( false, "Key ''DiffusionTransportModel'' is missing!" );
        }
        elementPhysicsModel =
            ElementPhysicsModelFactory::createElementPhysicsModel( transportModel_db );
        transportModel = AMP::dynamic_pointer_cast<DiffusionTransportModel>( elementPhysicsModel );
    }

    AMP_INSIST( transportModel.get() != nullptr, "NULL transport model" );

    // next create a ElementOperation object
    AMP_INSIST( input_db->keyExists( "DiffusionElement" ), "Key ''DiffusionElement'' is missing!" );
    AMP::shared_ptr<AMP::Operator::ElementOperation> diffusionLinElem =
        ElementOperationFactory::createElementOperation(
            input_db->getDatabase( "DiffusionElement" ) );

    // now create the linear diffusion operator
    AMP::shared_ptr<AMP::Database> diffusionLinFEOp_db;
    if ( input_db->getString( "name" ) == "DiffusionLinearFEOperator" ) {
        diffusionLinFEOp_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    AMP_INSIST( diffusionLinFEOp_db.get() != nullptr,
                "Error: The database object for DiffusionLinearFEOperator is NULL" );

    AMP::shared_ptr<AMP::Operator::DiffusionLinearFEOperatorParameters> diffusionOpParams(
        new AMP::Operator::DiffusionLinearFEOperatorParameters( diffusionLinFEOp_db ) );
    diffusionOpParams->d_transportModel = transportModel;
    diffusionOpParams->d_elemOp         = diffusionLinElem;
    diffusionOpParams->d_Mesh           = meshAdapter;
    diffusionOpParams->d_inDofMap =
        AMP::Discretization::simpleDOFManager::create( meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    diffusionOpParams->d_outDofMap =
        AMP::Discretization::simpleDOFManager::create( meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    AMP::shared_ptr<AMP::Operator::DiffusionLinearFEOperator> diffusionOp(
        new AMP::Operator::DiffusionLinearFEOperator( diffusionOpParams ) );

    AMP::LinearAlgebra::Matrix::shared_ptr matrix = diffusionOp->getMatrix();
    matrix->makeConsistent();

    return diffusionOp;
}


AMP::Operator::Operator::shared_ptr
OperatorBuilder::createVolumeIntegralOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                               AMP::shared_ptr<AMP::Database>
                                                   input_db,
                                               AMP::shared_ptr<AMP::Operator::ElementPhysicsModel>
                                                   elementPhysicsModel )
{
    AMP::shared_ptr<AMP::Operator::SourcePhysicsModel> sourcePhysicsModel;

    if ( elementPhysicsModel.get() != nullptr ) {
        sourcePhysicsModel =
            AMP::dynamic_pointer_cast<AMP::Operator::SourcePhysicsModel>( elementPhysicsModel );
    } else {
        if ( input_db->keyExists( "SourcePhysicsModel" ) ) {
            AMP::shared_ptr<AMP::Database> sourceModel_db =
                input_db->getDatabase( "SourcePhysicsModel" );
            elementPhysicsModel =
                ElementPhysicsModelFactory::createElementPhysicsModel( sourceModel_db );
            sourcePhysicsModel =
                AMP::dynamic_pointer_cast<SourcePhysicsModel>( elementPhysicsModel );
        }
    }

    // next create a ElementOperation object
    AMP_INSIST( input_db->keyExists( "SourceElement" ), "Key ''SourceElement'' is missing!" );
    AMP::shared_ptr<AMP::Operator::ElementOperation> sourceNonlinearElem =
        ElementOperationFactory::createElementOperation( input_db->getDatabase( "SourceElement" ) );

    // now create the nonlinear source operator
    AMP::shared_ptr<AMP::Database> sourceNLinFEOp_db;
    if ( input_db->getString( "name" ) == "VolumeIntegralOperator" ) {
        sourceNLinFEOp_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    AMP::shared_ptr<AMP::Operator::VolumeIntegralOperatorParameters> volumeIntegralParameters(
        new AMP::Operator::VolumeIntegralOperatorParameters( input_db ) );
    volumeIntegralParameters->d_sourcePhysicsModel = sourcePhysicsModel;
    volumeIntegralParameters->d_elemOp             = sourceNonlinearElem;
    volumeIntegralParameters->d_Mesh               = meshAdapter;
    AMP::shared_ptr<AMP::Operator::VolumeIntegralOperator> nonlinearSourceOp(
        new AMP::Operator::VolumeIntegralOperator( volumeIntegralParameters ) );

    return nonlinearSourceOp;
}


AMP::Operator::Operator::shared_ptr OperatorBuilder::createNonlinearDiffusionOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    AMP::shared_ptr<AMP::Database>
        input_db,
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel>
        elementPhysicsModel )
{

    // first create a DiffusionTransportModel
    AMP::shared_ptr<AMP::Operator::DiffusionTransportModel> transportModel;

    if ( elementPhysicsModel.get() != nullptr ) {
        transportModel = AMP::dynamic_pointer_cast<AMP::Operator::DiffusionTransportModel>(
            elementPhysicsModel );
    } else {
        AMP::shared_ptr<AMP::Database> transportModel_db;
        if ( input_db->keyExists( "DiffusionTransportModel" ) ) {
            transportModel_db = input_db->getDatabase( "DiffusionTransportModel" );
        } else {
            AMP_INSIST( false, "Key ''DiffusionTransportModel'' is missing!" );
        }
        elementPhysicsModel =
            ElementPhysicsModelFactory::createElementPhysicsModel( transportModel_db );
        transportModel = AMP::dynamic_pointer_cast<DiffusionTransportModel>( elementPhysicsModel );
    }

    AMP_INSIST( transportModel.get() != nullptr, "NULL transport model" );

    // next create an ElementOperation object
    AMP_INSIST( input_db->keyExists( "DiffusionElement" ), "Key ''DiffusionElement'' is missing!" );
    AMP::shared_ptr<AMP::Operator::ElementOperation> diffusionNonlinearElem =
        ElementOperationFactory::createElementOperation(
            input_db->getDatabase( "DiffusionElement" ) );

    // now create the nonlinear diffusion operator parameters
    AMP::shared_ptr<AMP::Database> diffusionNLinFEOp_db;
    if ( input_db->getString( "name" ) == "DiffusionNonlinearFEOperator" ) {
        diffusionNLinFEOp_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }
    AMP_INSIST( diffusionNLinFEOp_db.get() != nullptr,
                "Error: The database object for DiffusionNonlinearFEOperator is NULL" );
    AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperatorParameters> diffusionNLOpParams(
        new AMP::Operator::DiffusionNonlinearFEOperatorParameters( diffusionNLinFEOp_db ) );
    diffusionNLOpParams->d_transportModel = transportModel;
    diffusionNLOpParams->d_elemOp         = diffusionNonlinearElem;
    diffusionNLOpParams->d_Mesh           = meshAdapter;

    // populate the parameters with frozen active variable vectors

    // nullify vectors in parameters
    diffusionNLOpParams->d_FrozenTemperature.reset();
    diffusionNLOpParams->d_FrozenConcentration.reset();
    diffusionNLOpParams->d_FrozenBurnup.reset();

    // create variables and vectors for frozen material inputs
    AMP::shared_ptr<AMP::Database> active_db =
        diffusionNLinFEOp_db->getDatabase( "ActiveInputVariables" );
    AMP::Discretization::DOFManager::shared_ptr NodalScalarDOF =
        AMP::Discretization::simpleDOFManager::create( meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    std::string name;
    AMP::LinearAlgebra::Variable::shared_ptr tVar;
    AMP::LinearAlgebra::Vector::shared_ptr tVec;
    AMP::LinearAlgebra::Variable::shared_ptr cVar;
    AMP::LinearAlgebra::Vector::shared_ptr cVec;
    AMP::LinearAlgebra::Variable::shared_ptr bVar;
    AMP::LinearAlgebra::Vector::shared_ptr bVec;
    name = active_db->getStringWithDefault( "Temperature", "not_specified" );
    if ( name != "not_specified" ) {
        tVar.reset( new AMP::LinearAlgebra::Variable( name ) );
        tVec = AMP::LinearAlgebra::createVector( NodalScalarDOF, tVar, true );
        if ( diffusionNLinFEOp_db->getBoolWithDefault( "FreezeTemperature", false ) )
            diffusionNLOpParams->d_FrozenTemperature = tVec;
    }
    name = active_db->getStringWithDefault( "Concentration", "not_specified" );
    if ( name != "not_specified" ) {
        cVar.reset( new AMP::LinearAlgebra::Variable( name ) );
        cVec = AMP::LinearAlgebra::createVector( NodalScalarDOF, cVar, true );
        if ( diffusionNLinFEOp_db->getBoolWithDefault( "FreezeConcentration", false ) )
            diffusionNLOpParams->d_FrozenConcentration = cVec;
    }
    name = active_db->getStringWithDefault( "Burnup", "not_specified" );
    if ( name != "not_specified" ) {
        bVar.reset( new AMP::LinearAlgebra::Variable( name ) );
        bVec = AMP::LinearAlgebra::createVector( NodalScalarDOF, bVar, true );
        if ( diffusionNLinFEOp_db->getBoolWithDefault( "FreezeBurnup", false ) )
            diffusionNLOpParams->d_FrozenBurnup = bVec;
    }

    // create the nonlinear diffusion operator
    AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> nonlinearDiffusionOp(
        new AMP::Operator::DiffusionNonlinearFEOperator( diffusionNLOpParams ) );

    return nonlinearDiffusionOp;
}


AMP::Operator::Operator::shared_ptr OperatorBuilder::createNonlinearFickSoretOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    std::string operatorName,
    AMP::shared_ptr<AMP::Database>
        input_db,
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel>,
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModelFactory>
        localModelFactory )
{
    AMP::shared_ptr<Operator> retOperator;
    AMP_INSIST( input_db.get() != nullptr, "NULL database object passed" );

    AMP::shared_ptr<AMP::Database> operator_db =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>( input_db->getDatabase( operatorName ) );
    AMP_INSIST( operator_db.get() != nullptr, "NULL database object passed" );

    std::string fickOperatorName  = operator_db->getString( "FickOperator" );
    std::string soretOperatorName = operator_db->getString( "SoretOperator" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> fickPhysicsModel;
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> soretPhysicsModel;


    AMP::Operator::Operator::shared_ptr fickOperator = OperatorBuilder::createOperator(
        meshAdapter, fickOperatorName, input_db, fickPhysicsModel, localModelFactory );
    AMP_INSIST(
        fickOperator.get() != nullptr,
        "Error: unable to create Fick operator in OperatorBuilder::createFickSoretOperator" );

    AMP::Operator::Operator::shared_ptr soretOperator = OperatorBuilder::createOperator(
        meshAdapter, soretOperatorName, input_db, soretPhysicsModel, localModelFactory );

    AMP_INSIST(
        soretOperator.get() != nullptr,
        "Error: unable to create Soret operator in OperatorBuilder::createFickSoretOperator" );

    AMP::shared_ptr<AMP::Database> db = AMP::dynamic_pointer_cast<AMP::Database>( input_db );
    AMP::shared_ptr<FickSoretNonlinearFEOperatorParameters> params(
        new FickSoretNonlinearFEOperatorParameters( db ) );
    params->d_Mesh = meshAdapter;
    params->d_FickOperator =
        AMP::dynamic_pointer_cast<DiffusionNonlinearFEOperator>( fickOperator );
    params->d_SoretOperator =
        AMP::dynamic_pointer_cast<DiffusionNonlinearFEOperator>( soretOperator );
    params->d_name = operatorName;
    FickSoretNonlinearFEOperator::shared_ptr fsOp( new FickSoretNonlinearFEOperator( params ) );

    return fsOp;
}


AMP::Operator::Operator::shared_ptr
OperatorBuilder::createLinearMechanicsOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                                AMP::shared_ptr<AMP::Database>
                                                    input_db,
                                                AMP::shared_ptr<AMP::Operator::ElementPhysicsModel>
                                                    elementPhysicsModel )
{

    // first create a MechanicsMaterialModel
    if ( elementPhysicsModel.get() == nullptr ) {
        AMP_INSIST( input_db->keyExists( "MechanicsMaterialModel" ),
                    "Key ''MechanicsMaterialModel'' is missing!" );

        AMP::shared_ptr<AMP::Database> materialModel_db =
            input_db->getDatabase( "MechanicsMaterialModel" );
        elementPhysicsModel =
            ElementPhysicsModelFactory::createElementPhysicsModel( materialModel_db );
    }

    AMP_INSIST( elementPhysicsModel.get() != nullptr, "NULL material model" );

    // next create a ElementOperation object
    AMP_INSIST( input_db->keyExists( "MechanicsElement" ), "Key ''MechanicsElement'' is missing!" );
    AMP::shared_ptr<AMP::Operator::ElementOperation> mechanicsLinElem =
        ElementOperationFactory::createElementOperation(
            input_db->getDatabase( "MechanicsElement" ) );

    // now create the linear mechanics operator
    AMP::shared_ptr<AMP::Database> mechanicsLinFEOp_db;
    if ( input_db->getString( "name" ) == "MechanicsLinearFEOperator" ) {
        mechanicsLinFEOp_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    AMP_INSIST( mechanicsLinFEOp_db.get() != nullptr,
                "Error: The database object for MechanicsLinearFEOperator is NULL" );

    AMP::shared_ptr<AMP::Operator::MechanicsLinearFEOperatorParameters> mechanicsOpParams(
        new AMP::Operator::MechanicsLinearFEOperatorParameters( mechanicsLinFEOp_db ) );
    mechanicsOpParams->d_materialModel =
        AMP::dynamic_pointer_cast<MechanicsMaterialModel>( elementPhysicsModel );
    mechanicsOpParams->d_elemOp = mechanicsLinElem;
    mechanicsOpParams->d_Mesh   = meshAdapter;
    mechanicsOpParams->d_inDofMap =
        AMP::Discretization::simpleDOFManager::create( meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );
    mechanicsOpParams->d_outDofMap =
        AMP::Discretization::simpleDOFManager::create( meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );

    AMP::shared_ptr<AMP::Operator::MechanicsLinearFEOperator> mechanicsOp(
        new AMP::Operator::MechanicsLinearFEOperator( mechanicsOpParams ) );

    return mechanicsOp;
}


AMP::Operator::Operator::shared_ptr OperatorBuilder::createNonlinearMechanicsOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    AMP::shared_ptr<AMP::Database>
        input_db,
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel>
        elementPhysicsModel )
{

    // first create a MechanicsMaterialModel
    if ( elementPhysicsModel.get() == nullptr ) {
        AMP_INSIST( input_db->keyExists( "MechanicsMaterialModel" ),
                    "Key ''MechanicsMaterialModel'' is missing!" );

        AMP::shared_ptr<AMP::Database> transportModel_db =
            input_db->getDatabase( "MechanicsMaterialModel" );
        elementPhysicsModel =
            ElementPhysicsModelFactory::createElementPhysicsModel( transportModel_db );
    }

    AMP_INSIST( elementPhysicsModel.get() != nullptr, "NULL material model" );

    // next create a ElementOperation object
    AMP_INSIST( input_db->keyExists( "MechanicsElement" ), "Key ''MechanicsElement'' is missing!" );
    AMP::shared_ptr<AMP::Operator::ElementOperation> mechanicsElem =
        ElementOperationFactory::createElementOperation(
            input_db->getDatabase( "MechanicsElement" ) );

    // now create the nonlinear mechanics operator
    AMP::shared_ptr<AMP::Database> mechanicsFEOp_db;
    if ( input_db->getString( "name" ) == "MechanicsNonlinearFEOperator" ) {
        mechanicsFEOp_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    AMP_INSIST( mechanicsFEOp_db.get() != nullptr,
                "Error: The database object for MechanicsNonlinearFEOperator is NULL" );

    AMP::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperatorParameters> mechanicsOpParams(
        new AMP::Operator::MechanicsNonlinearFEOperatorParameters( mechanicsFEOp_db ) );
    mechanicsOpParams->d_materialModel =
        AMP::dynamic_pointer_cast<MechanicsMaterialModel>( elementPhysicsModel );
    mechanicsOpParams->d_elemOp = mechanicsElem;
    mechanicsOpParams->d_Mesh   = meshAdapter;
    mechanicsOpParams->d_dofMap[Mechanics::DISPLACEMENT] =
        AMP::Discretization::simpleDOFManager::create( meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );
    mechanicsOpParams->d_dofMap[Mechanics::TEMPERATURE] =
        AMP::Discretization::simpleDOFManager::create( meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    mechanicsOpParams->d_dofMap[Mechanics::BURNUP] =
        AMP::Discretization::simpleDOFManager::create( meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    mechanicsOpParams->d_dofMap[Mechanics::OXYGEN_CONCENTRATION] =
        AMP::Discretization::simpleDOFManager::create( meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    mechanicsOpParams->d_dofMap[Mechanics::LHGR] =
        AMP::Discretization::simpleDOFManager::create( meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );

    AMP::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> mechanicsOp(
        new AMP::Operator::MechanicsNonlinearFEOperator( mechanicsOpParams ) );

    return mechanicsOp;
}

AMP::Operator::Operator::shared_ptr OperatorBuilder::createLinearNavierStokesLSWFOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    AMP::shared_ptr<AMP::Database>
        input_db,
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel>
        elementPhysicsModel )
{

    if ( elementPhysicsModel.get() == nullptr ) {
        AMP_INSIST( input_db->keyExists( "FlowTransportModel" ),
                    "Key ''FlowTransportModel'' is missing!" );

        AMP::shared_ptr<AMP::Database> transportModel_db =
            input_db->getDatabase( "FlowTransportModel" );
        elementPhysicsModel =
            ElementPhysicsModelFactory::createElementPhysicsModel( transportModel_db );
    }

    AMP_INSIST( elementPhysicsModel.get() != nullptr, "NULL transport model" );

    // next create a ElementOperation object
    AMP_INSIST( input_db->keyExists( "FlowElement" ), "Key ''FlowElement'' is missing!" );
    AMP::shared_ptr<AMP::Operator::ElementOperation> flowLinElem =
        ElementOperationFactory::createElementOperation( input_db->getDatabase( "FlowElement" ) );

    // now create the linear flow operator
    AMP::shared_ptr<AMP::Database> flowLinFEOp_db;
    if ( input_db->getString( "name" ) == "NavierStokesLSWFLinearFEOperator" ) {
        flowLinFEOp_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    AMP_INSIST( flowLinFEOp_db.get() != nullptr,
                "Error: The database object for FlowLinearFEOperator is NULL" );

    AMP::shared_ptr<AMP::Operator::NavierStokesLinearFEOperatorParameters> flowOpParams(
        new AMP::Operator::NavierStokesLinearFEOperatorParameters( flowLinFEOp_db ) );
    flowOpParams->d_transportModel =
        AMP::dynamic_pointer_cast<FlowTransportModel>( elementPhysicsModel );
    flowOpParams->d_elemOp   = flowLinElem;
    flowOpParams->d_Mesh     = meshAdapter;
    flowOpParams->d_inDofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 10, true );
    flowOpParams->d_outDofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 10, true );

    AMP::shared_ptr<AMP::Operator::NavierStokesLSWFLinearFEOperator> flowOp(
        new AMP::Operator::NavierStokesLSWFLinearFEOperator( flowOpParams ) );

    return flowOp;
}


AMP::Operator::Operator::shared_ptr OperatorBuilder::createNonlinearNavierStokesLSWFOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    AMP::shared_ptr<AMP::Database>
        input_db,
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel>
        elementPhysicsModel )
{

    if ( elementPhysicsModel.get() == nullptr ) {
        AMP_INSIST( input_db->keyExists( "FlowTransportModel" ),
                    "Key ''FlowTransportModel'' is missing!" );

        AMP::shared_ptr<AMP::Database> transportModel_db =
            input_db->getDatabase( "FlowTransportModel" );
        elementPhysicsModel =
            ElementPhysicsModelFactory::createElementPhysicsModel( transportModel_db );
    }

    AMP_INSIST( elementPhysicsModel.get() != nullptr, "NULL material model" );

    // next create a ElementOperation object
    AMP_INSIST( input_db->keyExists( "FlowElement" ), "Key ''FlowElement'' is missing!" );
    AMP::shared_ptr<AMP::Operator::ElementOperation> flowElem =
        ElementOperationFactory::createElementOperation( input_db->getDatabase( "FlowElement" ) );

    // now create the nonlinear mechanics operator
    AMP::shared_ptr<AMP::Database> flowFEOp_db;
    if ( input_db->getString( "name" ) == "NavierStokesLSWFFEOperator" ) {
        flowFEOp_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    AMP_INSIST( flowFEOp_db.get() != nullptr,
                "Error: The database object for FlowNonlinearFEOperator is NULL" );

    AMP::shared_ptr<AMP::Operator::NavierStokesLSWFFEOperatorParameters> flowOpParams(
        new AMP::Operator::NavierStokesLSWFFEOperatorParameters( flowFEOp_db ) );
    flowOpParams->d_transportModel =
        AMP::dynamic_pointer_cast<FlowTransportModel>( elementPhysicsModel );
    flowOpParams->d_elemOp = flowElem;
    flowOpParams->d_Mesh   = meshAdapter;
    flowOpParams->d_dofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 10, true );
    AMP::shared_ptr<AMP::Operator::NavierStokesLSWFFEOperator> flowOp(
        new AMP::Operator::NavierStokesLSWFFEOperator( flowOpParams ) );

    return flowOp;
}


AMP::Operator::Operator::shared_ptr
OperatorBuilder::createMassLinearFEOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                             AMP::shared_ptr<AMP::Database>
                                                 input_db,
                                             AMP::shared_ptr<AMP::Operator::ElementPhysicsModel>
                                                 elementPhysicsModel )
{

    // first create a MassDensityModel
    AMP::shared_ptr<AMP::Operator::MassDensityModel> densityModel;

    if ( elementPhysicsModel.get() != nullptr ) {
        densityModel =
            AMP::dynamic_pointer_cast<AMP::Operator::MassDensityModel>( elementPhysicsModel );
    } else {
        AMP_INSIST( input_db->keyExists( "MassDensityModel" ),
                    "Key ''MassDensityModel'' is missing!" );
        AMP::shared_ptr<AMP::Database> densityModel_db =
            input_db->getDatabase( "MassDensityModel" );
        elementPhysicsModel =
            ElementPhysicsModelFactory::createElementPhysicsModel( densityModel_db );
        densityModel = AMP::dynamic_pointer_cast<MassDensityModel>( elementPhysicsModel );
    }

    AMP_INSIST( densityModel.get() != nullptr, "NULL density model" );

    // next create a ElementOperation object
    AMP_INSIST( input_db->keyExists( "MassElement" ), "Key ''MassElement'' is missing!" );
    AMP::shared_ptr<AMP::Operator::ElementOperation> densityLinElem =
        ElementOperationFactory::createElementOperation( input_db->getDatabase( "MassElement" ) );

    // now create the linear density operator
    AMP::shared_ptr<AMP::Database> densityLinFEOp_db;
    if ( input_db->getString( "name" ) == "MassLinearFEOperator" ) {
        densityLinFEOp_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    AMP_INSIST( densityLinFEOp_db.get() != nullptr,
                "Error: The database object for MassLinearFEOperator is NULL" );

    AMP::shared_ptr<AMP::Operator::MassLinearFEOperatorParameters> densityOpParams(
        new AMP::Operator::MassLinearFEOperatorParameters( densityLinFEOp_db ) );
    densityOpParams->d_densityModel = densityModel;
    densityOpParams->d_elemOp       = densityLinElem;
    densityOpParams->d_Mesh         = meshAdapter;
    densityOpParams->d_inDofMap =
        AMP::Discretization::simpleDOFManager::create( meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    densityOpParams->d_outDofMap =
        AMP::Discretization::simpleDOFManager::create( meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    AMP::shared_ptr<AMP::Operator::MassLinearFEOperator> densityOp(
        new AMP::Operator::MassLinearFEOperator( densityOpParams ) );

    return densityOp;
}


AMP::Operator::Operator::shared_ptr
OperatorBuilder::createLinearBVPOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                          std::string operatorName,
                                          AMP::shared_ptr<AMP::Database>
                                              input_db,
                                          AMP::shared_ptr<AMP::Operator::ElementPhysicsModel>
                                              elementPhysicsModel,
                                          AMP::shared_ptr<AMP::Operator::ElementPhysicsModelFactory>
                                              localModelFactory )
{
    AMP::shared_ptr<Operator> retOperator;
    AMP_INSIST( input_db.get() != nullptr, "NULL database object passed" );

    AMP::shared_ptr<AMP::Database> operator_db =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>( input_db->getDatabase( operatorName ) );
    AMP_INSIST( operator_db.get() != nullptr, "NULL database object passed" );

    // create the volume operator
    std::string volumeOperatorName = operator_db->getString( "VolumeOperator" );
    // if this flag is true the same local physics model will be used for both boundary and volume
    // operators
    bool useSameLocalModelForVolumeAndBoundaryOperators =
        operator_db->getBoolWithDefault( "useSameLocalModelForVolumeAndBoundaryOperators", false );

    AMP::Operator::Operator::shared_ptr volumeOperator = OperatorBuilder::createOperator(
        meshAdapter, volumeOperatorName, input_db, elementPhysicsModel, localModelFactory );

    AMP::shared_ptr<AMP::Operator::LinearOperator> volumeLinearOp =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearOperator>( volumeOperator );
    AMP_INSIST(
        volumeLinearOp.get() != nullptr,
        "Error: unable to create linear operator in OperatorBuilder::createLinearBVPOperator" );

    // create the boundary operator
    std::string boundaryOperatorName = operator_db->getString( "BoundaryOperator" );
    AMP::shared_ptr<AMP::Database> boundaryOperator_db =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>(
            input_db->getDatabase( boundaryOperatorName ) );
    AMP_INSIST( boundaryOperator_db.get() != nullptr,
                "NULL database object passed for boundary operator" );

    boundaryOperator_db->putBool( "isAttachedToVolumeOperator", true );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> boundaryLocalModel;

    if ( useSameLocalModelForVolumeAndBoundaryOperators ) {
        boundaryLocalModel = elementPhysicsModel;
    }

    AMP::shared_ptr<AMP::Operator::BoundaryOperator> boundaryOperator =
        OperatorBuilder::createBoundaryOperator( meshAdapter,
                                                 boundaryOperatorName,
                                                 input_db,
                                                 volumeLinearOp,
                                                 boundaryLocalModel,
                                                 localModelFactory );

    AMP::shared_ptr<AMP::Operator::BVPOperatorParameters> bvpOperatorParams(
        new AMP::Operator::BVPOperatorParameters( input_db ) );
    bvpOperatorParams->d_volumeOperator   = volumeOperator;
    bvpOperatorParams->d_boundaryOperator = boundaryOperator;

    retOperator.reset( new AMP::Operator::LinearBVPOperator( bvpOperatorParams ) );

    return retOperator;
}


AMP::Operator::Operator::shared_ptr OperatorBuilder::createNonlinearBVPOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    std::string operatorName,
    AMP::shared_ptr<AMP::Database>
        input_db,
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel>
        elementPhysicsModel,
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModelFactory>
        localModelFactory )
{
    AMP::shared_ptr<Operator> retOperator;
    AMP_INSIST( input_db.get() != nullptr, "NULL database object passed" );

    AMP::shared_ptr<AMP::Database> operator_db =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>( input_db->getDatabase( operatorName ) );
    AMP_INSIST( operator_db.get() != nullptr, "NULL database object passed" );

    // create the volume operator
    std::string volumeOperatorName = operator_db->getString( "VolumeOperator" );
    // if this flag is true the same local physics model will be used for both boundary and volume
    // operators
    bool useSameLocalModelForVolumeAndBoundaryOperators =
        operator_db->getBoolWithDefault( "useSameLocalModelForVolumeAndBoundaryOperators", false );

    AMP::Operator::Operator::shared_ptr volumeOperator = OperatorBuilder::createOperator(
        meshAdapter, volumeOperatorName, input_db, elementPhysicsModel, localModelFactory );
    AMP_INSIST( volumeOperator.get() != nullptr,
                "Error: unable to create nonlinear operator in "
                "OperatorBuilder::createNonlinearBVPOperator" );

    // create the boundary operator
    std::string boundaryOperatorName = operator_db->getString( "BoundaryOperator" );
    AMP::shared_ptr<AMP::Database> boundaryOperator_db =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>(
            input_db->getDatabase( boundaryOperatorName ) );
    AMP_INSIST( boundaryOperator_db.get() != nullptr,
                "NULL database object passed for boundary operator" );

    boundaryOperator_db->putBool( "isAttachedToVolumeOperator", true );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> boundaryLocalModel;

    if ( useSameLocalModelForVolumeAndBoundaryOperators ) {
        boundaryLocalModel = elementPhysicsModel;
    }

    AMP::shared_ptr<AMP::Operator::BoundaryOperator> boundaryOperator =
        OperatorBuilder::createBoundaryOperator( meshAdapter,
                                                 boundaryOperatorName,
                                                 input_db,
                                                 volumeOperator,
                                                 boundaryLocalModel,
                                                 localModelFactory );

    AMP::shared_ptr<AMP::Operator::BVPOperatorParameters> bvpOperatorParams(
        new AMP::Operator::BVPOperatorParameters( input_db ) );
    bvpOperatorParams->d_volumeOperator   = volumeOperator;
    bvpOperatorParams->d_boundaryOperator = boundaryOperator;

    retOperator.reset( new AMP::Operator::NonlinearBVPOperator( bvpOperatorParams ) );

    return retOperator;
}


AMP::Operator::Operator::shared_ptr OperatorBuilder::createFlowFrapconJacobian(
    AMP::Mesh::Mesh::shared_ptr meshAdapter, AMP::shared_ptr<AMP::Database> input_db )
{

    // now create the flow frapcon operator
    AMP::shared_ptr<AMP::Database> flowOp_db;
    if ( input_db->getString( "name" ) == "FlowFrapconJacobian" ) {
        flowOp_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    AMP_INSIST( flowOp_db.get() != nullptr,
                "Error: The database object for FlowFrapconJacobian is NULL" );

    AMP::shared_ptr<AMP::Operator::FlowFrapconJacobianParameters> flowOpParams(
        new AMP::Operator::FlowFrapconJacobianParameters( flowOp_db ) );
    flowOpParams->d_Mesh = meshAdapter;
    AMP::shared_ptr<AMP::Operator::FlowFrapconJacobian> flowOp(
        new AMP::Operator::FlowFrapconJacobian( flowOpParams ) );

    return flowOp;
}

AMP::shared_ptr<BoundaryOperator>
OperatorBuilder::createBoundaryOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                         std::string boundaryOperatorName,
                                         AMP::shared_ptr<AMP::Database>
                                             input_db,
                                         AMP::Operator::Operator::shared_ptr volumeOperator,
                                         AMP::shared_ptr<AMP::Operator::ElementPhysicsModel>
                                             elementPhysicsModel,
                                         AMP::shared_ptr<AMP::Operator::ElementPhysicsModelFactory>
                                             localModelFactory )
{
    AMP::shared_ptr<BoundaryOperator> retOperator;
    AMP_INSIST( input_db.get() != nullptr, "NULL database object passed" );

    AMP::shared_ptr<AMP::Database> operator_db = AMP::dynamic_pointer_cast<AMP::InputDatabase>(
        input_db->getDatabase( boundaryOperatorName ) );
    AMP_INSIST( operator_db.get() != nullptr,
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

        AMP::shared_ptr<AMP::Database> localModel_db = input_db->getDatabase( localModelName );
        AMP_INSIST( localModel_db.get() != nullptr,
                    "Error:: OperatorBuilder::createOperator(): No local model database "
                    "entry with given name exists in input database" );

        // if a non-NULL factory is being supplied through the argument list
        // use it, else call the AMP ElementPhysicsModelFactory interface
        if ( localModelFactory.get() != nullptr ) {
            elementPhysicsModel = localModelFactory->createElementPhysicsModel( localModel_db );
        } else {
            elementPhysicsModel =
                ElementPhysicsModelFactory::createElementPhysicsModel( localModel_db );
        }

        AMP_INSIST( elementPhysicsModel.get() != nullptr,
                    "Error:: OperatorBuilder::createOperator(): local model creation failed" );
    }

    std::string boundaryType = operator_db->getString( "name" );

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


AMP::shared_ptr<BoundaryOperator> OperatorBuilder::createColumnBoundaryOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    std::string boundaryOperatorName,
    AMP::shared_ptr<AMP::Database>
        input_db,
    AMP::Operator::Operator::shared_ptr volumeOperator,
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel>
        elementPhysicsModel,
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModelFactory>
        localModelFactory )
{

    AMP_INSIST( input_db.get() != nullptr, "NULL database object passed" );

    AMP::shared_ptr<AMP::Database> operator_db = AMP::dynamic_pointer_cast<AMP::InputDatabase>(
        input_db->getDatabase( boundaryOperatorName ) );
    AMP_INSIST( operator_db.get() != nullptr,
                "Error: OperatorBuilder::createBoundaryOperator(): "
                "database object with given name not in database" );

    int numberOfBoundaryOperators =
        operator_db->getIntegerWithDefault( "numberOfBoundaryOperators", 1 );

    auto boundaryOps = new std::string[numberOfBoundaryOperators];
    operator_db->getStringArray( "boundaryOperators", boundaryOps, numberOfBoundaryOperators );

    AMP::shared_ptr<OperatorParameters> params( new OperatorParameters( operator_db ) );
    params->d_Mesh = meshAdapter;

    AMP::shared_ptr<AMP::Operator::ColumnBoundaryOperator> columnBoundaryOperator(
        new AMP::Operator::ColumnBoundaryOperator( params ) );

    for ( int i = 0; i < numberOfBoundaryOperators; i++ ) {
        AMP::shared_ptr<BoundaryOperator> bcOperator =
            OperatorBuilder::createBoundaryOperator( meshAdapter,
                                                     boundaryOps[i],
                                                     input_db,
                                                     volumeOperator,
                                                     elementPhysicsModel,
                                                     localModelFactory );
        AMP_ASSERT( bcOperator );
        columnBoundaryOperator->append( bcOperator );
    }
    delete[] boundaryOps;

    return columnBoundaryOperator;
}


AMP::shared_ptr<BoundaryOperator> OperatorBuilder::createDirichletMatrixCorrection(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    AMP::shared_ptr<AMP::Database>
        input_db,
    AMP::Operator::Operator::shared_ptr volumeOperator,
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> )
{
    AMP::shared_ptr<BoundaryOperator> retOperator;
    AMP::shared_ptr<AMP::Operator::LinearOperator> linearOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearOperator>( volumeOperator );
    AMP::shared_ptr<AMP::Operator::DirichletMatrixCorrectionParameters> matrixCorrectionParameters(
        new AMP::Operator::DirichletMatrixCorrectionParameters( input_db ) );
    matrixCorrectionParameters->d_variable    = linearOperator->getOutputVariable();
    matrixCorrectionParameters->d_inputMatrix = linearOperator->getMatrix();
    matrixCorrectionParameters->d_Mesh        = meshAdapter;

    retOperator.reset( new AMP::Operator::DirichletMatrixCorrection( matrixCorrectionParameters ) );

    return retOperator;
}


AMP::shared_ptr<BoundaryOperator>
OperatorBuilder::createMassMatrixCorrection( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                             AMP::shared_ptr<AMP::Database>
                                                 input_db,
                                             AMP::Operator::Operator::shared_ptr volumeOperator,
                                             AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> )
{
    AMP::shared_ptr<BoundaryOperator> retOperator;
    AMP::shared_ptr<AMP::Operator::LinearOperator> linearOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearOperator>( volumeOperator );
    AMP::shared_ptr<AMP::Operator::DirichletMatrixCorrectionParameters> matrixCorrectionParameters(
        new AMP::Operator::DirichletMatrixCorrectionParameters( input_db ) );
    matrixCorrectionParameters->d_variable    = linearOperator->getOutputVariable();
    matrixCorrectionParameters->d_inputMatrix = linearOperator->getMatrix();
    matrixCorrectionParameters->d_Mesh        = meshAdapter;

    retOperator.reset( new AMP::Operator::MassMatrixCorrection( matrixCorrectionParameters ) );

    return retOperator;
}


AMP::shared_ptr<BoundaryOperator>
OperatorBuilder::createRobinMatrixCorrection( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                              AMP::shared_ptr<AMP::Database>
                                                  input_db,
                                              AMP::Operator::Operator::shared_ptr volumeOperator,
                                              AMP::shared_ptr<AMP::Operator::ElementPhysicsModel>
                                                  elementPhysicsModel )
{
    AMP::shared_ptr<BoundaryOperator> retOperator;
    AMP::shared_ptr<AMP::Operator::LinearOperator> linearOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearOperator>( volumeOperator );
    AMP::shared_ptr<AMP::Operator::RobinMatrixCorrectionParameters> matrixCorrectionParameters(
        new AMP::Operator::RobinMatrixCorrectionParameters( input_db ) );
    matrixCorrectionParameters->d_variable    = linearOperator->getOutputVariable();
    matrixCorrectionParameters->d_inputMatrix = linearOperator->getMatrix();
    matrixCorrectionParameters->d_Mesh        = meshAdapter;

    if ( elementPhysicsModel.get() != nullptr ) {

        AMP::shared_ptr<AMP::Operator::RobinPhysicsModel> robinPhysicsModel;
        robinPhysicsModel =
            AMP::dynamic_pointer_cast<AMP::Operator::RobinPhysicsModel>( elementPhysicsModel );
        matrixCorrectionParameters->d_robinPhysicsModel = robinPhysicsModel;
    }

    retOperator.reset( new AMP::Operator::RobinMatrixCorrection( matrixCorrectionParameters ) );

    return retOperator;
}


AMP::shared_ptr<BoundaryOperator>
OperatorBuilder::createRobinVectorCorrection( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                              AMP::shared_ptr<AMP::Database>
                                                  input_db,
                                              AMP::Operator::Operator::shared_ptr volumeOperator,
                                              AMP::shared_ptr<AMP::Operator::ElementPhysicsModel>
                                                  elementPhysicsModel )
{
    AMP::shared_ptr<BoundaryOperator> retOperator;
    AMP::shared_ptr<AMP::Operator::NeumannVectorCorrectionParameters> vectorCorrectionParameters(
        new AMP::Operator::NeumannVectorCorrectionParameters( input_db ) );
    vectorCorrectionParameters->d_variable = volumeOperator->getOutputVariable();
    vectorCorrectionParameters->d_Mesh     = meshAdapter;

    if ( elementPhysicsModel.get() != nullptr ) {
        AMP::shared_ptr<AMP::Operator::RobinPhysicsModel> robinPhysicsModel;
        robinPhysicsModel =
            AMP::dynamic_pointer_cast<AMP::Operator::RobinPhysicsModel>( elementPhysicsModel );
        vectorCorrectionParameters->d_robinPhysicsModel = robinPhysicsModel;
    }

    retOperator.reset( new AMP::Operator::RobinVectorCorrection( vectorCorrectionParameters ) );

    return retOperator;
}


AMP::shared_ptr<BoundaryOperator>
OperatorBuilder::createNeumannVectorCorrection( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                                AMP::shared_ptr<AMP::Database>
                                                    input_db,
                                                AMP::Operator::Operator::shared_ptr volumeOperator,
                                                AMP::shared_ptr<AMP::Operator::ElementPhysicsModel>
                                                    elementPhysicsModel )
{
    AMP::shared_ptr<BoundaryOperator> retOperator;
    AMP::shared_ptr<AMP::Operator::NeumannVectorCorrectionParameters> vectorCorrectionParameters(
        new AMP::Operator::NeumannVectorCorrectionParameters( input_db ) );
    vectorCorrectionParameters->d_variable = volumeOperator->getOutputVariable();
    vectorCorrectionParameters->d_Mesh     = meshAdapter;

    if ( elementPhysicsModel.get() != nullptr && input_db->isDatabase( "RobinPhysicsModel" ) ) {
        AMP::shared_ptr<AMP::Operator::RobinPhysicsModel> robinPhysicsModel;
        robinPhysicsModel =
            AMP::dynamic_pointer_cast<AMP::Operator::RobinPhysicsModel>( elementPhysicsModel );
        vectorCorrectionParameters->d_robinPhysicsModel = robinPhysicsModel;
    }

    retOperator.reset( new AMP::Operator::NeumannVectorCorrection( vectorCorrectionParameters ) );

    return retOperator;
}


AMP::shared_ptr<BoundaryOperator> OperatorBuilder::createDirichletVectorCorrection(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    AMP::shared_ptr<AMP::Database>
        input_db,
    AMP::Operator::Operator::shared_ptr volumeOperator,
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> )
{
    AMP::shared_ptr<BoundaryOperator> retOperator;
    AMP::shared_ptr<AMP::Operator::DirichletVectorCorrectionParameters> vectorCorrectionParameters(
        new AMP::Operator::DirichletVectorCorrectionParameters( input_db ) );
    vectorCorrectionParameters->d_variable = volumeOperator->getOutputVariable();
    vectorCorrectionParameters->d_Mesh     = meshAdapter;

    retOperator.reset( new AMP::Operator::DirichletVectorCorrection( vectorCorrectionParameters ) );

    return retOperator;
}


AMP::shared_ptr<BoundaryOperator> OperatorBuilder::createDirichletVectorCorrection(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    AMP::shared_ptr<AMP::Database>
        input_db,
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> )
{
    AMP::shared_ptr<BoundaryOperator> retOperator;
    AMP::shared_ptr<AMP::Operator::DirichletVectorCorrectionParameters> vectorCorrectionParameters(
        new AMP::Operator::DirichletVectorCorrectionParameters( input_db ) );
    vectorCorrectionParameters->d_Mesh = meshAdapter;

    retOperator.reset( new AMP::Operator::DirichletVectorCorrection( vectorCorrectionParameters ) );

    return retOperator;
}


AMP::shared_ptr<BoundaryOperator> OperatorBuilder::createPressureBoundaryOperator(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    AMP::shared_ptr<AMP::Database>
        input_db,
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> )
{
    AMP::shared_ptr<BoundaryOperator> retOperator;
    AMP::shared_ptr<AMP::Operator::OperatorParameters> params(
        new AMP::Operator::OperatorParameters( input_db ) );
    params->d_Mesh = meshAdapter;

    retOperator.reset( new AMP::Operator::PressureBoundaryOperator( params ) );

    return retOperator;
}


AMP::shared_ptr<Operator> OperatorBuilder::createOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter1,
                                                           AMP::Mesh::Mesh::shared_ptr meshAdapter2,
                                                           AMP::AMP_MPI comm,
                                                           AMP::shared_ptr<AMP::Database>
                                                               input_db )
{
    AMP::shared_ptr<Operator> retOperator;

    std::string name = input_db->getString( "name" );
    if ( name == "MapSurface" ) {
        AMP::shared_ptr<AMP::Operator::MapOperatorParameters> mapOperatorParameters(
            new AMP::Operator::MapOperatorParameters( input_db ) );
        mapOperatorParameters->d_Mesh    = meshAdapter1;
        mapOperatorParameters->d_MapMesh = meshAdapter2;
        mapOperatorParameters->d_MapComm = comm;
        retOperator.reset( new AMP::Operator::MapSurface( mapOperatorParameters ) );
    }
    return retOperator;
}


#endif
}
}
