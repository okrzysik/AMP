#include "AMP/operators/OperatorBuilder.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/StructuredMeshHelper.h"
#include "AMP/operators/IdentityOperator.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/VectorBuilder.h"

#ifdef AMP_USE_LIBMESH
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

#include "ProfilerApp.h"

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

void OperatorBuilder::setNestedOperatorMemoryLocations(
    std::shared_ptr<AMP::Database> input_db,
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

std::vector<std::string>
OperatorBuilder::getActiveVariables( std::shared_ptr<const AMP::Database> db,
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


std::shared_ptr<Operator>
OperatorBuilder::createOperator( std::shared_ptr<OperatorParameters> in_params )
{
    PROFILE( "OperatorBuilder::createOperator" );

    std::shared_ptr<Operator> retOperator;

    AMP_INSIST( in_params, "ERROR: OperatorBuilder::createOperator has NULL input" );
    AMP_INSIST( in_params->d_db,
                "ERROR: OperatorBuilder::createOperator has NULL database pointer in in_params" );

    std::string name = in_params->d_db->getString( "name" );

    resetOperation( IdentityOperator );
#ifdef AMP_USE_LIBMESH
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


std::shared_ptr<Operator> OperatorBuilder::createOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                                           const std::string &operatorName,
                                                           std::shared_ptr<AMP::Database> input_db )
{
    PROFILE( "OperatorBuilder::createOperator" );

    std::shared_ptr<ElementPhysicsModel> model;
    std::shared_ptr<ElementPhysicsModelFactory> factory;
    return createOperator( mesh, operatorName, input_db, model, factory );
}


std::shared_ptr<Operator>
OperatorBuilder::createOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                 const std::string &operatorName,
                                 std::shared_ptr<AMP::Database> input_db,
                                 std::shared_ptr<ElementPhysicsModel> &elementPhysicsModel,
                                 std::shared_ptr<ElementPhysicsModelFactory> localModelFactory )
{
    PROFILE( "OperatorBuilder::createOperator" );

    auto operator_db = input_db->getDatabase( operatorName );

    AMP_INSIST( operator_db,
                "Error:: OperatorBuilder::createOperator(): No operator database entry with "
                "given name exists in input database: " +
                    operatorName );

    // we create the element physics model if a database entry exists
    // and the incoming element physics model pointer is NULL
    if ( !elementPhysicsModel && operator_db->keyExists( "LocalModel" ) ) {
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
    std::shared_ptr<Operator> retOperator;
    if ( operatorType == "IdentityOperator" ) {
        retOperator = OperatorBuilder::createIdentityOperator( mesh, operator_db );
    } else if ( operatorType == "MechanicsLinearFEOperator" ) {
        retOperator = OperatorBuilder::createLinearMechanicsOperator(
            mesh, operator_db, elementPhysicsModel );
    } else if ( operatorType == "MechanicsNonlinearFEOperator" ) {
        retOperator = OperatorBuilder::createNonlinearMechanicsOperator(
            mesh, operator_db, elementPhysicsModel );
    } else if ( operatorType == "DiffusionLinearFEOperator" ) {
        retOperator = OperatorBuilder::createLinearDiffusionOperator(
            mesh, operator_db, elementPhysicsModel );
    } else if ( operatorType == "DiffusionNonlinearFEOperator" ) {
        retOperator = OperatorBuilder::createNonlinearDiffusionOperator(
            mesh, operator_db, elementPhysicsModel );
    } else if ( operatorType == "NavierStokesLSWFLinearFEOperator" ) {
        retOperator = OperatorBuilder::createLinearNavierStokesLSWFOperator(
            mesh, operator_db, elementPhysicsModel );
    } else if ( operatorType == "NavierStokesLSWFFEOperator" ) {
        retOperator = OperatorBuilder::createNonlinearNavierStokesLSWFOperator(
            mesh, operator_db, elementPhysicsModel );
    } else if ( operatorType == "FickSoretNonlinearFEOperator" ) {
        retOperator = OperatorBuilder::createNonlinearFickSoretOperator(
            mesh, operatorName, input_db, elementPhysicsModel, localModelFactory );
    } else if ( operatorType == "FlowFrapconOperator" ) {
        retOperator = OperatorBuilder::createFlowFrapconOperator( mesh, operator_db );
    } else if ( operatorType == "FlowFrapconJacobian" ) {
        retOperator = OperatorBuilder::createFlowFrapconJacobian( mesh, operator_db );
    } else if ( operatorType == "SubchannelTwoEqLinearOperator" ) {
        retOperator = OperatorBuilder::createSubchannelTwoEqLinearOperator(
            mesh, operator_db, elementPhysicsModel );
    } else if ( operatorType == "SubchannelTwoEqNonlinearOperator" ) {
        retOperator = OperatorBuilder::createSubchannelTwoEqNonlinearOperator(
            mesh, operator_db, elementPhysicsModel );
    } else if ( operatorType == "SubchannelFourEqNonlinearOperator" ) {
        retOperator = OperatorBuilder::createSubchannelFourEqNonlinearOperator(
            mesh, operator_db, elementPhysicsModel );
    } else if ( operatorType == "NeutronicsRhsOperator" ) {
        retOperator = OperatorBuilder::createNeutronicsRhsOperator( mesh, operator_db );
    } else if ( operatorType == "MassLinearFEOperator" ) {
        retOperator =
            OperatorBuilder::createMassLinearFEOperator( mesh, operator_db, elementPhysicsModel );
    } else if ( operatorType == "VolumeIntegralOperator" ) {
        retOperator =
            OperatorBuilder::createVolumeIntegralOperator( mesh, operator_db, elementPhysicsModel );
    } else if ( operatorType == "LinearBVPOperator" ) {
        // note that we pass in the full database here and not the operator db
        retOperator = OperatorBuilder::createLinearBVPOperator(
            mesh, operatorName, input_db, elementPhysicsModel, localModelFactory );
    } else if ( operatorType == "NonlinearBVPOperator" ) {
        // note that we pass in the full database here and not the operator db
        retOperator = OperatorBuilder::createNonlinearBVPOperator(
            mesh, operatorName, input_db, elementPhysicsModel, localModelFactory );
    } else if ( operatorType == "DirichletMatrixCorrection" ) {
    } else if ( operatorType == "DirichletVectorCorrection" ) {
        retOperator = OperatorBuilder::createDirichletVectorCorrection(
            mesh, operator_db, elementPhysicsModel );
    } else if ( operatorType == "PressureBoundaryOperator" ) {
        retOperator = OperatorBuilder::createPressureBoundaryOperator(
            mesh, operator_db, elementPhysicsModel );
    } else if ( operatorType == "NeumannVectorCorrection" ) {
    } else if ( operatorType == "RobinMatrixCorrection" ) {
    } else {
        AMP_ERROR( "Unknown operator: " + operatorName );
    }

    return retOperator;
}


#ifdef AMP_USE_LIBMESH


// Create the identity operator
Operator::shared_ptr
OperatorBuilder::createIdentityOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                         std::shared_ptr<AMP::Database> input_db )
{
    PROFILE( "OperatorBuilder::createIdentityOperator" );

    AMP_INSIST( input_db, "Error: The database object for SubchannelTwoEqLinearOperator is NULL" );

    auto params       = std::make_shared<OperatorParameters>( input_db );
    params->d_Mesh    = mesh;
    auto subchannelOp = std::make_shared<IdentityOperator>( params );

    return subchannelOp;
}


// Create the FlowFrapconOperator
Operator::shared_ptr
OperatorBuilder::createFlowFrapconOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                            std::shared_ptr<AMP::Database> input_db )
{
    PROFILE( "OperatorBuilder::createFlowFrapconOperator" );

    // now create the flow frapcon operator
    std::shared_ptr<AMP::Database> flowOp_db;
    if ( input_db->getString( "name" ) == "FlowFrapconOperator" ) {
        flowOp_db = input_db;
    } else {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }
    AMP_INSIST( flowOp_db, "Error: The database object for FlowFrapconOperator is NULL" );

    auto flowOpParams    = std::make_shared<FlowFrapconOperatorParameters>( flowOp_db );
    flowOpParams->d_Mesh = mesh;
    return std::make_shared<FlowFrapconOperator>( flowOpParams );
}

Operator::shared_ptr OperatorBuilder::createSubchannelTwoEqLinearOperator(
    std::shared_ptr<AMP::Mesh::Mesh> mesh,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> &elementPhysicsModel )
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
    subchannelParams->d_Mesh = mesh;
    subchannelParams->d_subchannelPhysicsModel = transportModel;

    return std::make_shared<SubchannelTwoEqLinearOperator>( subchannelParams );
}


Operator::shared_ptr OperatorBuilder::createSubchannelTwoEqNonlinearOperator(
    std::shared_ptr<AMP::Mesh::Mesh> mesh,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> &elementPhysicsModel )
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
    subchannelParams->d_Mesh = mesh;
    subchannelParams->d_subchannelPhysicsModel = transportModel;
    return std::make_shared<SubchannelTwoEqNonlinearOperator>( subchannelParams );
}


Operator::shared_ptr OperatorBuilder::createSubchannelFourEqNonlinearOperator(
    std::shared_ptr<AMP::Mesh::Mesh> mesh,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> &elementPhysicsModel )
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
    subchannelParams->d_Mesh = mesh;
    subchannelParams->d_subchannelPhysicsModel = transportModel;
    return std::make_shared<SubchannelFourEqNonlinearOperator>( subchannelParams );
}


Operator::shared_ptr
OperatorBuilder::createNeutronicsRhsOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                              std::shared_ptr<AMP::Database> input_db )
{
    PROFILE( "OperatorBuilder::createNeutronicsRhsOperator" );

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
    neutronicsOpParams->d_Mesh = mesh;
    return std::make_shared<NeutronicsRhs>( neutronicsOpParams );
}


Operator::shared_ptr OperatorBuilder::createLinearDiffusionOperator(
    std::shared_ptr<AMP::Mesh::Mesh> mesh,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> &elementPhysicsModel )
{
    PROFILE( "OperatorBuilder::createLinearDiffusionOperator" );

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


Operator::shared_ptr OperatorBuilder::createVolumeIntegralOperator(
    std::shared_ptr<AMP::Mesh::Mesh> mesh,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> &elementPhysicsModel )
{
    PROFILE( "OperatorBuilder::createVolumeIntegralOperator" );

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
    if ( input_db->getString( "name" ) != "VolumeIntegralOperator" ) {
        AMP_INSIST( input_db->keyExists( "name" ), "Key ''name'' is missing!" );
    }

    auto volumeIntegralParameters = std::make_shared<VolumeIntegralOperatorParameters>( input_db );
    volumeIntegralParameters->d_sourcePhysicsModel = sourcePhysicsModel;
    volumeIntegralParameters->d_elemOp             = sourceNonlinearElem;
    volumeIntegralParameters->d_Mesh               = mesh;
    return std::make_shared<VolumeIntegralOperator>( volumeIntegralParameters );
}


Operator::shared_ptr OperatorBuilder::createNonlinearDiffusionOperator(
    std::shared_ptr<AMP::Mesh::Mesh> mesh,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> &elementPhysicsModel )
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
    auto diffusionNonlinearElem = ElementOperationFactory::createElementOperation(
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
    diffusionNLOpParams->d_Mesh           = mesh;

    // populate the parameters with frozen active variable vectors

    // nullify vectors in parameters
    diffusionNLOpParams->d_FrozenVecs.clear();

    // create variables and vectors for frozen material inputs
    auto NodalScalarDOF = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    for ( auto name : getActiveVariables( diffusionNLinFEOp_db, "ActiveInputVariables" ) ) {
        auto var = std::make_shared<AMP::LinearAlgebra::Variable>( name );
        auto vec = AMP::LinearAlgebra::createVector(
            NodalScalarDOF, var, true, diffusionNLOpParams->d_memory_location );
        if ( diffusionNLinFEOp_db->getWithDefault<bool>( "Freeze" + name, false ) )
            diffusionNLOpParams->d_FrozenVecs[name] = vec;
    }

    // create the nonlinear diffusion operator
    return std::make_shared<DiffusionNonlinearFEOperator>( diffusionNLOpParams );
}

Operator::shared_ptr OperatorBuilder::createNonlinearFickSoretOperator(
    std::shared_ptr<AMP::Mesh::Mesh> mesh,
    std::string operatorName,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> &,
    std::shared_ptr<ElementPhysicsModelFactory> localModelFactory )
{
    std::shared_ptr<Operator> retOperator;
    AMP_INSIST( input_db, "NULL database object passed" );

    auto operator_db = input_db->getDatabase( operatorName );
    AMP_INSIST( operator_db, "NULL database object passed" );

    // operator names
    auto fickOperatorName  = operator_db->getString( "FickOperator" );
    auto soretOperatorName = operator_db->getString( "SoretOperator" );

    // Ensure consistency of operator memory locations
    OperatorBuilder::setNestedOperatorMemoryLocations(
        input_db, operatorName, { fickOperatorName, soretOperatorName } );

    std::shared_ptr<ElementPhysicsModel> fickPhysicsModel;
    std::shared_ptr<ElementPhysicsModel> soretPhysicsModel;

    auto fickOperator = OperatorBuilder::createOperator(
        mesh, fickOperatorName, input_db, fickPhysicsModel, localModelFactory );
    AMP_INSIST(
        fickOperator,
        "Error: unable to create Fick operator in OperatorBuilder::createFickSoretOperator" );

    auto soretOperator = OperatorBuilder::createOperator(
        mesh, soretOperatorName, input_db, soretPhysicsModel, localModelFactory );

    AMP_INSIST(
        soretOperator,
        "Error: unable to create Soret operator in OperatorBuilder::createFickSoretOperator" );

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


Operator::shared_ptr OperatorBuilder::createLinearMechanicsOperator(
    std::shared_ptr<AMP::Mesh::Mesh> mesh,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> &elementPhysicsModel )
{

    // first create a MechanicsMaterialModel
    if ( !elementPhysicsModel ) {
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
    mechanicsOpParams->d_Mesh     = mesh;
    mechanicsOpParams->d_inDofMap = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 3, true );
    mechanicsOpParams->d_outDofMap = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 3, true );

    return std::make_shared<MechanicsLinearFEOperator>( mechanicsOpParams );
}


Operator::shared_ptr OperatorBuilder::createNonlinearMechanicsOperator(
    std::shared_ptr<AMP::Mesh::Mesh> mesh,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> &elementPhysicsModel )
{

    // first create a MechanicsMaterialModel
    if ( !elementPhysicsModel ) {
        AMP_INSIST( input_db->keyExists( "MechanicsMaterialModel" ),
                    "Key ''MechanicsMaterialModel'' is missing!" );

        auto transportModel_db = input_db->getDatabase( "MechanicsMaterialModel" );
        elementPhysicsModel =
            ElementPhysicsModelFactory::createElementPhysicsModel( transportModel_db );
    }
    AMP_INSIST( elementPhysicsModel, "NULL material model" );

    // next create a ElementOperation object
    AMP_INSIST( input_db->keyExists( "MechanicsElement" ), "Key ''MechanicsElement'' is missing!" );
    auto mechanicsElem = ElementOperationFactory::createElementOperation(
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

Operator::shared_ptr OperatorBuilder::createLinearNavierStokesLSWFOperator(
    std::shared_ptr<AMP::Mesh::Mesh> mesh,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> &elementPhysicsModel )
{

    if ( !elementPhysicsModel ) {
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
    flowOpParams->d_Mesh     = mesh;
    flowOpParams->d_inDofMap = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 10, true );
    flowOpParams->d_outDofMap = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 10, true );

    return std::make_shared<NavierStokesLSWFLinearFEOperator>( flowOpParams );
}


Operator::shared_ptr OperatorBuilder::createNonlinearNavierStokesLSWFOperator(
    std::shared_ptr<AMP::Mesh::Mesh> mesh,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> &elementPhysicsModel )
{

    if ( !elementPhysicsModel ) {
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
    flowOpParams->d_Mesh   = mesh;
    flowOpParams->d_dofMap = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 10, true );
    return std::make_shared<NavierStokesLSWFFEOperator>( flowOpParams );
}


Operator::shared_ptr OperatorBuilder::createMassLinearFEOperator(
    std::shared_ptr<AMP::Mesh::Mesh> mesh,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> &elementPhysicsModel )
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
    densityOpParams->d_Mesh         = mesh;
    densityOpParams->d_inDofMap     = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    densityOpParams->d_outDofMap = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    return std::make_shared<MassLinearFEOperator>( densityOpParams );
}


Operator::shared_ptr OperatorBuilder::createLinearBVPOperator(
    std::shared_ptr<AMP::Mesh::Mesh> mesh,
    std::string operatorName,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> &elementPhysicsModel,
    std::shared_ptr<ElementPhysicsModelFactory> localModelFactory )
{
    PROFILE( "OperatorBuilder::createLinearBVPOperator" );

    AMP_INSIST( input_db, "NULL database object passed" );

    auto operator_db = input_db->getDatabase( operatorName );
    AMP_INSIST( operator_db, "NULL database object passed" );

    // names of internal operators
    auto volumeOperatorName   = operator_db->getString( "VolumeOperator" );
    auto boundaryOperatorName = operator_db->getString( "BoundaryOperator" );

    // Ensure consistency of operator memory locations
    OperatorBuilder::setNestedOperatorMemoryLocations(
        input_db, operatorName, { volumeOperatorName, boundaryOperatorName } );

    // create the volume operator
    auto volumeOperator = OperatorBuilder::createOperator(
        mesh, volumeOperatorName, input_db, elementPhysicsModel, localModelFactory );
    auto volumeLinearOp = std::dynamic_pointer_cast<LinearOperator>( volumeOperator );
    AMP_INSIST(
        volumeLinearOp,
        "Error: unable to create linear operator in OperatorBuilder::createLinearBVPOperator" );

    // create the boundary operator
    auto boundaryOperator_db = input_db->getDatabase( boundaryOperatorName );
    AMP_INSIST( boundaryOperator_db, "NULL database object passed for boundary operator" );
    boundaryOperator_db->putScalar(
        "isAttachedToVolumeOperator", true, Units(), Database::Check::Overwrite );

    std::shared_ptr<ElementPhysicsModel> boundaryLocalModel;

    // Use the same local physics model for both boundary and volume operators?
    bool flag = operator_db->getWithDefault<bool>( "useSameLocalModelForVolumeAndBoundaryOperators",
                                                   false );
    if ( flag ) {
        boundaryLocalModel = elementPhysicsModel;
    }

    auto boundaryOperator = OperatorBuilder::createBoundaryOperator( mesh,
                                                                     boundaryOperatorName,
                                                                     input_db,
                                                                     volumeLinearOp,
                                                                     boundaryLocalModel,
                                                                     localModelFactory );

    auto bvpOperatorParams                = std::make_shared<BVPOperatorParameters>( operator_db );
    bvpOperatorParams->d_volumeOperator   = volumeOperator;
    bvpOperatorParams->d_boundaryOperator = boundaryOperator;

    return std::make_shared<LinearBVPOperator>( bvpOperatorParams );
}


Operator::shared_ptr OperatorBuilder::createNonlinearBVPOperator(
    std::shared_ptr<AMP::Mesh::Mesh> mesh,
    std::string operatorName,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<ElementPhysicsModel> &elementPhysicsModel,
    std::shared_ptr<ElementPhysicsModelFactory> localModelFactory )
{
    AMP_INSIST( input_db, "NULL database object passed" );

    auto operator_db = input_db->getDatabase( operatorName );
    AMP_INSIST( operator_db, "NULL database object passed" );

    // operator names
    auto volumeOperatorName   = operator_db->getString( "VolumeOperator" );
    auto boundaryOperatorName = operator_db->getString( "BoundaryOperator" );

    // Ensure consistency of operator memory locations
    OperatorBuilder::setNestedOperatorMemoryLocations(
        input_db, operatorName, { volumeOperatorName, boundaryOperatorName } );

    // create the volume operator
    auto volumeOperator = OperatorBuilder::createOperator(
        mesh, volumeOperatorName, input_db, elementPhysicsModel, localModelFactory );
    AMP_INSIST( volumeOperator,
                "Error: unable to create nonlinear operator in "
                "OperatorBuilder::createNonlinearBVPOperator" );

    // create the boundary operator
    auto boundaryOperator_db = input_db->getDatabase( boundaryOperatorName );
    AMP_INSIST( boundaryOperator_db, "NULL database object passed for boundary operator" );

    boundaryOperator_db->putScalar(
        "isAttachedToVolumeOperator", true, Units(), Database::Check::Overwrite );

    std::shared_ptr<ElementPhysicsModel> boundaryLocalModel;

    // Use the same local physics model for both boundary and volume operators?
    bool flag = operator_db->getWithDefault<bool>( "useSameLocalModelForVolumeAndBoundaryOperators",
                                                   false );
    if ( flag ) {
        boundaryLocalModel = elementPhysicsModel;
    }

    auto boundaryOperator = OperatorBuilder::createBoundaryOperator( mesh,
                                                                     boundaryOperatorName,
                                                                     input_db,
                                                                     volumeOperator,
                                                                     boundaryLocalModel,
                                                                     localModelFactory );

    auto bvpOperatorParams                = std::make_shared<BVPOperatorParameters>( input_db );
    bvpOperatorParams->d_volumeOperator   = volumeOperator;
    bvpOperatorParams->d_boundaryOperator = boundaryOperator;

    return std::make_shared<NonlinearBVPOperator>( bvpOperatorParams );
}


Operator::shared_ptr
OperatorBuilder::createFlowFrapconJacobian( std::shared_ptr<AMP::Mesh::Mesh> mesh,
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
    flowOpParams->d_Mesh = mesh;
    return std::make_shared<FlowFrapconJacobian>( flowOpParams );
}

std::shared_ptr<BoundaryOperator> OperatorBuilder::createBoundaryOperator(
    std::shared_ptr<AMP::Mesh::Mesh> mesh,
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

        auto localModel_db = input_db->getDatabase( localModelName );
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
            mesh, operator_db, volumeOperator, elementPhysicsModel );
    } else if ( boundaryType == "MassMatrixCorrection" ) {
        retOperator =
            createMassMatrixCorrection( mesh, operator_db, volumeOperator, elementPhysicsModel );
    } else if ( boundaryType == "RobinMatrixCorrection" ) {
        // in this case the volume operator has to be a linear operator
        retOperator =
            createRobinMatrixCorrection( mesh, operator_db, volumeOperator, elementPhysicsModel );
    } else if ( boundaryType == "RobinVectorCorrection" ) {
        retOperator =
            createRobinVectorCorrection( mesh, operator_db, volumeOperator, elementPhysicsModel );
    } else if ( boundaryType == "NeumannVectorCorrection" ) {
        retOperator =
            createNeumannVectorCorrection( mesh, operator_db, volumeOperator, elementPhysicsModel );
    } else if ( boundaryType == "DirichletVectorCorrection" ) {
        // in this case the volume operator has to be a nonlinear operator
        retOperator = createDirichletVectorCorrection(
            mesh, operator_db, volumeOperator, elementPhysicsModel );
    } else if ( boundaryType == "PressureBoundaryOperator" ) {
        retOperator = createPressureBoundaryOperator( mesh, operator_db, elementPhysicsModel );
    } else if ( boundaryType == "ColumnBoundaryOperator" ) {
        // note that the global input database is passed here instead of the operator
        // database
        retOperator = createColumnBoundaryOperator( mesh,
                                                    boundaryOperatorName,
                                                    input_db,
                                                    volumeOperator,
                                                    elementPhysicsModel,
                                                    localModelFactory );
    }

    return retOperator;
}


std::shared_ptr<BoundaryOperator> OperatorBuilder::createColumnBoundaryOperator(
    std::shared_ptr<AMP::Mesh::Mesh> mesh,
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

    // Ensure consistency of operator memory locations
    OperatorBuilder::setNestedOperatorMemoryLocations(
        input_db, boundaryOperatorName, boundaryOps );

    auto params    = std::make_shared<OperatorParameters>( operator_db );
    params->d_Mesh = mesh;

    auto columnBoundaryOperator = std::make_shared<ColumnBoundaryOperator>( params );

    for ( int i = 0; i < numberOfBoundaryOperators; i++ ) {
        auto model2     = elementPhysicsModel;
        auto bcOperator = OperatorBuilder::createBoundaryOperator(
            mesh, boundaryOps[i], input_db, volumeOperator, model2, localModelFactory );
        AMP_ASSERT( bcOperator );
        columnBoundaryOperator->append( bcOperator );
    }

    return columnBoundaryOperator;
}


std::shared_ptr<BoundaryOperator>
OperatorBuilder::createDirichletMatrixCorrection( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                                  std::shared_ptr<AMP::Database> input_db,
                                                  Operator::shared_ptr volumeOperator,
                                                  std::shared_ptr<ElementPhysicsModel> & )
{
    std::shared_ptr<BoundaryOperator> retOperator;
    auto linearOperator = std::dynamic_pointer_cast<LinearOperator>( volumeOperator );
    auto matrixCorrectionParameters =
        std::make_shared<DirichletMatrixCorrectionParameters>( input_db );
    matrixCorrectionParameters->d_variable    = linearOperator->getOutputVariable();
    matrixCorrectionParameters->d_inputMatrix = linearOperator->getMatrix();
    matrixCorrectionParameters->d_Mesh        = mesh;

    retOperator.reset( new DirichletMatrixCorrection( matrixCorrectionParameters ) );

    return retOperator;
}


std::shared_ptr<BoundaryOperator>
OperatorBuilder::createMassMatrixCorrection( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                             std::shared_ptr<AMP::Database> input_db,
                                             Operator::shared_ptr volumeOperator,
                                             std::shared_ptr<ElementPhysicsModel> & )
{
    auto linearOperator = std::dynamic_pointer_cast<LinearOperator>( volumeOperator );
    auto matrixCorrectionParameters =
        std::make_shared<DirichletMatrixCorrectionParameters>( input_db );
    matrixCorrectionParameters->d_variable    = linearOperator->getOutputVariable();
    matrixCorrectionParameters->d_inputMatrix = linearOperator->getMatrix();
    matrixCorrectionParameters->d_Mesh        = mesh;

    return std::make_shared<MassMatrixCorrection>( matrixCorrectionParameters );
}


std::shared_ptr<BoundaryOperator> OperatorBuilder::createRobinMatrixCorrection(
    std::shared_ptr<AMP::Mesh::Mesh> mesh,
    std::shared_ptr<AMP::Database> input_db,
    Operator::shared_ptr volumeOperator,
    std::shared_ptr<ElementPhysicsModel> &elementPhysicsModel )
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


std::shared_ptr<BoundaryOperator> OperatorBuilder::createRobinVectorCorrection(
    std::shared_ptr<AMP::Mesh::Mesh> mesh,
    std::shared_ptr<AMP::Database> input_db,
    Operator::shared_ptr volumeOperator,
    std::shared_ptr<ElementPhysicsModel> &elementPhysicsModel )
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


std::shared_ptr<BoundaryOperator> OperatorBuilder::createNeumannVectorCorrection(
    std::shared_ptr<AMP::Mesh::Mesh> mesh,
    std::shared_ptr<AMP::Database> input_db,
    Operator::shared_ptr volumeOperator,
    std::shared_ptr<ElementPhysicsModel> &elementPhysicsModel )
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


std::shared_ptr<BoundaryOperator>
OperatorBuilder::createDirichletVectorCorrection( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                                  std::shared_ptr<AMP::Database> input_db,
                                                  Operator::shared_ptr volumeOperator,
                                                  std::shared_ptr<ElementPhysicsModel> & )
{
    auto vectorCorrectionParameters =
        std::make_shared<DirichletVectorCorrectionParameters>( input_db );
    vectorCorrectionParameters->d_variable = volumeOperator->getOutputVariable();
    vectorCorrectionParameters->d_Mesh     = mesh;
    return std::make_shared<DirichletVectorCorrection>( vectorCorrectionParameters );
}


std::shared_ptr<BoundaryOperator>
OperatorBuilder::createDirichletVectorCorrection( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                                  std::shared_ptr<AMP::Database> input_db,
                                                  std::shared_ptr<ElementPhysicsModel> & )
{
    auto vectorCorrectionParameters =
        std::make_shared<DirichletVectorCorrectionParameters>( input_db );
    vectorCorrectionParameters->d_Mesh = mesh;
    return std::make_shared<DirichletVectorCorrection>( vectorCorrectionParameters );
}


std::shared_ptr<BoundaryOperator>
OperatorBuilder::createPressureBoundaryOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                                 std::shared_ptr<AMP::Database> input_db,
                                                 std::shared_ptr<ElementPhysicsModel> & )
{
    auto params    = std::make_shared<OperatorParameters>( input_db );
    params->d_Mesh = mesh;
    return std::make_shared<PressureBoundaryOperator>( params );
}


std::shared_ptr<Operator> OperatorBuilder::createOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh1,
                                                           std::shared_ptr<AMP::Mesh::Mesh> mesh2,
                                                           const AMP::AMP_MPI &comm,
                                                           std::shared_ptr<AMP::Database> input_db )
{
    std::shared_ptr<Operator> retOperator;
    std::string name = input_db->getString( "name" );
    if ( name == "MapSurface" ) {
        auto mapOperatorParameters       = std::make_shared<MapOperatorParameters>( input_db );
        mapOperatorParameters->d_Mesh    = mesh1;
        mapOperatorParameters->d_MapMesh = mesh2;
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
                                                std::shared_ptr<ElementPhysicsModel> & )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator>
OperatorBuilder::createNonlinearMechanicsOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                                   std::shared_ptr<AMP::Database>,
                                                   std::shared_ptr<ElementPhysicsModel> & )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator>
OperatorBuilder::createLinearDiffusionOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                                std::shared_ptr<AMP::Database>,
                                                std::shared_ptr<ElementPhysicsModel> & )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator>
OperatorBuilder::createNonlinearDiffusionOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                                   std::shared_ptr<AMP::Database>,
                                                   std::shared_ptr<ElementPhysicsModel> & )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator>
OperatorBuilder::createLinearNavierStokesLSWFOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                                       std::shared_ptr<AMP::Database>,
                                                       std::shared_ptr<ElementPhysicsModel> & )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator>
OperatorBuilder::createNonlinearNavierStokesLSWFOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                                          std::shared_ptr<AMP::Database>,
                                                          std::shared_ptr<ElementPhysicsModel> & )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator>
OperatorBuilder::createNonlinearFickSoretOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                                   std::string,
                                                   std::shared_ptr<AMP::Database>,
                                                   std::shared_ptr<ElementPhysicsModel> &,
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
                                                      std::shared_ptr<ElementPhysicsModel> & )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator>
OperatorBuilder::createSubchannelTwoEqNonlinearOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                                         std::shared_ptr<AMP::Database>,
                                                         std::shared_ptr<ElementPhysicsModel> & )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator>
OperatorBuilder::createSubchannelFourEqNonlinearOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                                          std::shared_ptr<AMP::Database>,
                                                          std::shared_ptr<ElementPhysicsModel> & )
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
                                             std::shared_ptr<ElementPhysicsModel> & )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator>
OperatorBuilder::createVolumeIntegralOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                               std::shared_ptr<AMP::Database>,
                                               std::shared_ptr<ElementPhysicsModel> & )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator>
OperatorBuilder::createLinearBVPOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                          std::string,
                                          std::shared_ptr<AMP::Database>,
                                          std::shared_ptr<ElementPhysicsModel> &,
                                          std::shared_ptr<ElementPhysicsModelFactory> )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<Operator>
OperatorBuilder::createNonlinearBVPOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                             std::string,
                                             std::shared_ptr<AMP::Database>,
                                             std::shared_ptr<ElementPhysicsModel> &,
                                             std::shared_ptr<ElementPhysicsModelFactory> )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<BoundaryOperator>
OperatorBuilder::createDirichletVectorCorrection( std::shared_ptr<AMP::Mesh::Mesh>,
                                                  std::shared_ptr<AMP::Database>,
                                                  std::shared_ptr<ElementPhysicsModel> & )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}
std::shared_ptr<BoundaryOperator>
OperatorBuilder::createPressureBoundaryOperator( std::shared_ptr<AMP::Mesh::Mesh>,
                                                 std::shared_ptr<AMP::Database>,
                                                 std::shared_ptr<ElementPhysicsModel> & )
{
    AMP_ERROR( "No libmesh" );
    return nullptr;
}


#endif

} // namespace AMP::Operator
