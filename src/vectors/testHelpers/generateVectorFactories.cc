#include "vectors/testHelpers/generateVectorFactories.h"

#include <vectors/testHelpers/StridedVectorSelector.h>
#include <vectors/testHelpers/VectorFactory.h>

#ifdef USE_EXT_PETSC
#include <vectors/testHelpers/petsc/PetscVectorFactory.h>
#endif

#ifdef USE_EXT_TRILINOS
#ifdef USE_TRILINOS_THYRA
#include <vectors/trilinos/epetra/ManagedEpetraVector.h>
#endif
#ifdef USE_TRILINOS_THYRA
#include <vectors/testHelpers/trilinos/thyra/ThyraVectorFactory.h>
#endif
#endif

#ifdef USE_OPENMP
#include "vectors/operations/OpenMP/VectorOperationsOpenMP.h"
#endif
#ifdef USE_CUDA
#include "vectors/data/cuda/VectorDataGPU.h"
#include "vectors/operations/cuda/VectorOperationsCuda.h"
#endif

#include <string>
#include <vector>


namespace AMP {
namespace LinearAlgebra {


// Trim string
static inline std::string trim( std::string s )
{
    s.erase(
        s.begin(),
        std::find_if( s.begin(), s.end(), std::not1( std::ptr_fun<int, int>( std::isspace ) ) ) );
    s.erase(
        std::find_if( s.rbegin(), s.rend(), std::not1( std::ptr_fun<int, int>( std::isspace ) ) )
            .base(),
        s.end() );
    return s;
}


// Split input arguments
static inline std::vector<std::string> splitArgs( const std::string &input )
{
    auto pos = input.find_first_of( '<' );
    if ( pos == std::string::npos )
        return std::vector<std::string>();
    std::string tmp = trim( input.substr( pos + 1 ) );
    tmp.resize( tmp.size() - 1 );
    tmp = trim( tmp );
    std::vector<std::string> args;
    while ( !tmp.empty() ) {
        pos       = 0;
        int count = 0;
        while ( pos < tmp.length() ) {
            if ( tmp[pos] == '<' )
                count++;
            if ( tmp[pos] == '>' )
                count--;
            if ( tmp[pos] == ',' && count == 0 )
                break;
            pos++;
        }
        args.push_back( tmp.substr( 0, pos ) );
        if ( pos + 1 < tmp.length() )
            tmp = trim( tmp.substr( pos + 1 ) );
        else
            tmp.clear();
    }
    return args;
}


// Convert strings to arguments of the desired type
int to_int( const std::string &s ) { return atoi( s.c_str() ); }
int to_bool( const std::string &s )
{
    if ( s == "false" )
        return false;
    else if ( s == "true" )
        return true;
    else
        AMP_ERROR( "Unknown value for bool" );
    return false;
}


// Generate a SimpleVectorFactory
template<typename TYPE, typename VecOps>
AMP::shared_ptr<VectorFactory>
generateSimpleVectorFactory( int N, bool global, const std::string &data )
{
    AMP::shared_ptr<VectorFactory> factory;
    if ( data == "cpu" ) {
        factory.reset(
            new SimpleVectorFactory<TYPE, VecOps, AMP::LinearAlgebra::VectorDataCPU<TYPE>>(
                N, global ) );
    } else if ( data == "gpu" ) {
#ifdef USE_CUDA
        factory.reset(
            new SimpleVectorFactory<TYPE, VecOps, AMP::LinearAlgebra::VectorDataGPU<TYPE>>(
                N, global ) );
#else
        AMP_ERROR( "gpu data is not supported without CUDA" );
#endif
    } else {
        AMP_ERROR( "Unknown VectorData" );
    }
    return factory;
}
template<typename TYPE>
AMP::shared_ptr<VectorFactory>
generateSimpleVectorFactory( int N, bool global, const std::string &ops, const std::string &data )
{
    AMP::shared_ptr<VectorFactory> factory;
    if ( ops == "default" ) {
        factory =
            generateSimpleVectorFactory<TYPE, AMP::LinearAlgebra::VectorOperationsDefault<TYPE>>(
                N, global, data );
    } else if ( ops == "openmp" ) {
#ifdef USE_OPENMP
        factory =
            generateSimpleVectorFactory<TYPE, AMP::LinearAlgebra::VectorOperationsOpenMP<TYPE>>(
                N, global, data );
#else
        AMP_ERROR( "openmp generators are not supported without OpenMP" );
#endif
    } else if ( ops == "cuda" ) {
#ifdef USE_CUDA
        factory = generateSimpleVectorFactory<TYPE, AMP::LinearAlgebra::VectorOperationsCuda<TYPE>>(
            N, global, data );
#else
        AMP_ERROR( "cuda generators are not supported without CUDA" );
#endif
    } else {
        AMP_ERROR( "Unknown VectorOperations" );
    }
    return factory;
}
AMP::shared_ptr<VectorFactory> generateSimpleVectorFactory(
    int N, bool global, const std::string &type, const std::string &ops, const std::string &data )
{
    AMP::shared_ptr<VectorFactory> factory;
    if ( type == "double" ) {
        factory = generateSimpleVectorFactory<double>( N, global, ops, data );
    } else if ( type == "float" ) {
        factory = generateSimpleVectorFactory<float>( N, global, ops, data );
    } else {
        AMP_ERROR( "Unknown VectorOperations" );
    }
    return factory;
}


AMP::shared_ptr<VectorFactory> generateVectorFactory( const std::string &name )
{
    auto pos                      = name.find_first_of( '<' );
    const std::string factoryName = name.substr( 0, pos );
    auto args                     = splitArgs( name );
    AMP::shared_ptr<VectorFactory> factory;
    if ( factoryName == "SimpleVectorFactory" ) {
        AMP_ASSERT( args.size() >= 2 );
        // Set default arguments
        if ( args.size() < 3 )
            args.emplace_back( "double" );
        if ( args.size() < 4 )
            args.emplace_back( "default" );
        if ( args.size() < 5 )
            args.emplace_back( "cpu" );
        factory = generateSimpleVectorFactory(
            to_int( args[0] ), to_bool( args[1] ), args[2], args[3], args[4] );
    } else if ( factoryName == "ArrayVectorFactory" ) {
        AMP_ASSERT( args.size() >= 3 );
        if ( args.size() == 3 )
            args.emplace_back( "double" );
        if ( args[3] == "double" ) {
            factory.reset( new ArrayVectorFactory<double>(
                to_int( args[0] ), to_int( args[1] ), to_bool( args[2] ) ) );
        } else if ( args[3] == "float" ) {
            factory.reset( new ArrayVectorFactory<float>(
                to_int( args[0] ), to_int( args[1] ), to_bool( args[2] ) ) );
        } else {
            AMP_ERROR( "Unknown type" );
        }
    } else if ( factoryName == "SimplePetscNativeFactory" ) {
        AMP_ASSERT( args.size() == 0 );
#if defined( USE_EXT_PETSC ) && defined( USE_EXT_TRILINOS )
        factory.reset( new SimplePetscNativeFactory() );
#else
        AMP_ERROR( "Generator is not valid without support for PETSc and Trilinos" );
#endif
    } else if ( factoryName == "MultiVectorFactory" ) {
        AMP_ASSERT( args.size() == 4 );
        factory.reset( new MultiVectorFactory( generateVectorFactory( args[0] ),
                                               to_int( args[1] ),
                                               generateVectorFactory( args[2] ),
                                               to_int( args[3] ) ) );
    } else if ( factoryName == "NativePetscVectorFactory" ) {
#if defined( USE_EXT_PETSC )
        factory.reset( new NativePetscVectorFactory() );
#else
        AMP_ERROR( "Generator is not valid without support for PETSc" );
#endif
    } else if ( factoryName == "SimpleManagedVectorFactory" ) {
        AMP_ASSERT( args.size() == 1 );
        if ( args[0] == "ManagedPetscVector" ) {
#if defined( USE_EXT_PETSC ) && defined( USE_EXT_TRILINOS )
            factory.reset( new SimpleManagedVectorFactory<ManagedPetscVector>() );
#else
            AMP_ERROR( "Generator is not valid without support for PETSc and Trilinos" );
#endif
        } else if ( args[0] == "ManagedEpetraVector" ) {
#ifdef USE_EXT_TRILINOS
            factory.reset( new SimpleManagedVectorFactory<ManagedEpetraVector>() );
#else
            AMP_ERROR( "Generator is not valid without support for Petsc" );
#endif
        } else {
            AMP_ERROR( "Unknown template argument for SimpleManagedVectorFactory" );
        }
    } else if ( factoryName == "NativeThyraFactory" ) {
        AMP_ASSERT( args.size() == 0 );
#ifdef USE_TRILINOS_THYRA
        factory.reset( new NativeThyraFactory() );
#else
        AMP_ERROR( "Generator is not valid without support for Thyra" );
#endif
    } else if ( factoryName == "ManagedThyraFactory" ) {
        AMP_ASSERT( args.size() == 1 );
#ifdef USE_TRILINOS_THYRA
        factory.reset( new ManagedThyraFactory( generateVectorFactory( args[0] ) ) );
#else
        AMP_ERROR( "Generator is not valid without support for Thyra" );
#endif
    } else if ( factoryName == "ManagedNativeThyraFactory" ) {
        AMP_ASSERT( args.size() == 1 );
#ifdef USE_TRILINOS_THYRA
        factory.reset( new ManagedNativeThyraFactory( generateVectorFactory( args[0] ) ) );
#else
        AMP_ERROR( "Generator is not valid without support for Thyra" );
#endif
    } else if ( factoryName == "NativeSundialsFactory" ) {
        AMP_ASSERT( args.size() == 0 );
#ifdef EXT_SUNDIALS
        AMP_ERROR( "Not implemented" );
#else
        AMP_ERROR( "Generator is not valid without support for Sundials" );
#endif
    } else if ( factoryName == "ManagedSundialsVectorFactory" ) {
        AMP_ASSERT( args.size() == 0 );
#ifdef EXT_SUNDIALS
        AMP_ERROR( "Not implemented" );
#else
        AMP_ERROR( "Generator is not valid without support for Sundials" );
#endif
    } else if ( factoryName == "ViewFactory" ) {
        AMP_ASSERT( args.size() == 2 );
        auto factory2 = generateVectorFactory( args[1] );
        if ( args[0] == "PetscVector" ) {
#ifdef USE_EXT_PETSC
            factory.reset( new ViewFactory<PetscVector>( factory2 ) );
#else
            AMP_ERROR( "Generator is not valid without support for Petsc" );
#endif
        } else {
            AMP_ERROR( "Unknown template argument for SimpleManagedVectorFactory" );
        }
    } else if ( factoryName == "CloneFactory" ) {
        AMP_ASSERT( args.size() == 1 );
        factory.reset( new CloneFactory( generateVectorFactory( args[0] ) ) );
    } else if ( factoryName == "StridedVectorFactory" ) {
        AMP_ASSERT( args.size() == 1 );
        factory.reset( new CloneFactory( generateVectorFactory( args[0] ) ) );
    } else {
        AMP_ERROR( "Unknown factory" );
    }
    return factory;
}
} // namespace LinearAlgebra
} // namespace AMP
