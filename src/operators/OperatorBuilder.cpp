#include "AMP/operators/OperatorBuilder.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/StructuredMeshHelper.h"
#include "AMP/operators/IdentityOperator.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/OperatorFactory.h"
#include "AMP/operators/boundary/BoundaryOperatorParameters.h"
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
#endif


#include <string>

#include "ProfilerApp.h"


namespace AMP::Operator::OperatorBuilder {


static std::shared_ptr<BoundaryOperator>
createBoundaryOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                        std::string boundaryOperatorName,
                        std::shared_ptr<AMP::Database> input_db,
                        AMP::Operator::Operator::shared_ptr volumeOperator );


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
            AMP_INSIST( inner_db, "OperatorBuilder: inner DB is null (" + innerName + ")" );
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

static Operator::shared_ptr createLinearBVPOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                                     std::string operatorName,
                                                     std::shared_ptr<AMP::Database> input_db )
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
    auto volumeOperator = createOperator( mesh, volumeOperatorName, input_db );
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

static Operator::shared_ptr createNonlinearBVPOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                                        std::string operatorName,
                                                        std::shared_ptr<AMP::Database> input_db )
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
    auto volumeOperator = createOperator( mesh, volumeOperatorName, input_db );

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


#endif


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

    int N = operator_db->getWithDefault<int>( "numberOfBoundaryOperators", 1 );

    auto boundaryOps = operator_db->getVector<std::string>( "boundaryOperators" );
    AMP_ASSERT( N == (int) boundaryOps.size() );

    // Ensure consistency of operator memory locations
    setNestedOperatorMemoryLocations( input_db, boundaryOperatorName, boundaryOps );

    auto params = std::make_shared<OperatorParameters>( operator_db, mesh );

    auto columnBoundaryOperator = std::make_shared<ColumnBoundaryOperator>( params );

    for ( int i = 0; i < N; i++ ) {
        auto bcOperator = createBoundaryOperator( mesh, boundaryOps[i], input_db, volumeOperator );
        AMP_ASSERT( bcOperator );
        columnBoundaryOperator->append( bcOperator );
    }

    return columnBoundaryOperator;
}


/********************************************************
 * createOperator                                        *
 ********************************************************/
std::shared_ptr<Operator> createOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                          const std::string &operatorName,
                                          std::shared_ptr<AMP::Database> input_db )
{
    PROFILE( "createOperator" );
    auto operator_db = input_db;
    if ( input_db->keyExists( operatorName ) )
        operator_db = input_db->getDatabase( operatorName );
    AMP_INSIST( operator_db, "operator database is null" );
    auto operatorType = operator_db->getString( "name" );
    if ( operatorType == "null" ) {
#ifdef AMP_USE_LIBMESH
    } else if ( operatorType == "FickSoretNonlinearFEOperator" ) {
        return createNonlinearFickSoretOperator( mesh, operatorName, input_db );
    } else if ( operatorType == "LinearBVPOperator" ) {
        // note that we pass in the full database here and not the operator db
        return createLinearBVPOperator( mesh, operatorName, input_db );
    } else if ( operatorType == "NonlinearBVPOperator" ) {
        // note that we pass in the full database here and not the operator db
        return createNonlinearBVPOperator( mesh, operatorName, input_db );
#endif
    } else if ( OperatorFactory::exists( operatorType ) ) {
        // Use the OperatorFactory to create the operator
        auto params = std::make_shared<OperatorParameters>( operator_db, mesh );
        return OperatorFactory::create( params );
    }
    AMP_ERROR( "Unable to create operator " + operatorType );
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
    auto boundaryType = operator_db->getString( "name" );
    if ( boundaryType == "ColumnBoundaryOperator" ) {
        // note that the global input database is passed here instead of the operator
        // database
        return createColumnBoundaryOperator( mesh, boundaryOperatorName, input_db, volumeOperator );
    } else if ( OperatorFactory::exists( boundaryType ) ) {
        // Use the OperatorFactory to create the operator
        auto params                  = std::make_shared<BoundaryOperatorParameters>( operator_db );
        params->d_Mesh               = mesh;
        params->d_volumeOperator     = volumeOperator;
        std::shared_ptr<Operator> op = OperatorFactory::create( params );
        auto op2                     = std::dynamic_pointer_cast<BoundaryOperator>( op );
        AMP_ASSERT( op2 );
        return op2;
    }
    return nullptr;
}


} // namespace AMP::Operator::OperatorBuilder
