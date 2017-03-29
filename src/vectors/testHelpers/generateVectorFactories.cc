#include "vectors/testHelpers/generateVectorFactories.h"

#include <vectors/testHelpers/VectorFactory.h>
#include <vectors/testHelpers/StridedVectorSelector.h>

#ifdef USE_EXT_PETSC
#include <vectors/testHelpers/petsc/PetscVectorFactory.h>
#endif

#ifdef USE_EXT_TRILINOS
#include <vectors/trilinos/ManagedEpetraVector.h>
#ifdef USE_TRILINOS_THYRA
#include <vectors/testHelpers/trilinos/thyra/ThyraVectorFactory.h>
#endif
#endif

#include <vector>
#include <string>


namespace AMP {
namespace LinearAlgebra {


// Trim string
static inline std::string trim( std::string s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))));
    s.erase(std::find_if(s.rbegin(), s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}


// Split input arguments
static inline std::vector<std::string> splitArgs( const std::string &input )
{
    auto pos = input.find_first_of('<');
    if ( pos == std::string::npos )
        return std::vector<std::string>();
    std::string tmp = trim( input.substr(pos+1) );
    tmp.resize( tmp.size()-1 );
    tmp = trim(tmp);
    std::vector<std::string> args;
    while ( !tmp.empty() ) {
        pos = 0;
        int count = 0;
        while ( pos < tmp.length() ) {
            if ( tmp[pos] == '<' )
                count++;
            if ( tmp[pos] == '>' )
                count--;
            if ( tmp[pos]==',' && count==0 )
                break;
            pos++;
        }
        args.push_back( tmp.substr(0,pos) );
        if ( pos+1 < tmp.length() )
            tmp = trim( tmp.substr(pos+1) );
        else
            tmp.clear();
    }
    return args;
}


// Convert strings to arguments of the desired type
int to_int( const std::string& s )
{
    return atoi( s.c_str() );
}
int to_bool( const std::string& s )
{
    if ( s == "false" )
        return false;
    else if ( s == "true" )
        return true;
    else
        AMP_ERROR("Unknown value for bool");
    return false;
}


AMP::shared_ptr<VectorFactory> generateVectorFactory( const std::string& name )
{
    auto pos = name.find_first_of('<');
    const std::string factoryName = name.substr( 0, pos );
    auto args = splitArgs( name );
    AMP::shared_ptr<VectorFactory> factory;
    if ( factoryName == "SimpleVectorFactory" ) {
        AMP_ASSERT(args.size()>=2);
        if ( args.size()==2 )
            args.push_back("double");
        if ( args[2] == "double" ) {
            factory.reset( new SimpleVectorFactory<double>( to_int(args[0]), to_bool(args[1]) ) );
        } else if ( args[2] == "float" ) {
            factory.reset( new SimpleVectorFactory<float>( to_int(args[0]), to_bool(args[1]) ) );
        } else {
            AMP_ERROR("Unknown type");
        }
    } else if ( factoryName == "ArrayVectorFactory" ) {
        AMP_ASSERT(args.size()>=3);
        if ( args.size()==3 )
            args.push_back("double");
        if ( args[3] == "double" ) {
            factory.reset( new ArrayVectorFactory<double>( to_int(args[0]), to_int(args[1]), to_bool(args[2]) ) );
        } else if ( args[3] == "float" ) {
            factory.reset( new ArrayVectorFactory<float>( to_int(args[0]), to_int(args[1]), to_bool(args[2]) ) );
        } else {
            AMP_ERROR("Unknown type");
        }
    } else if ( factoryName == "SimplePetscNativeFactory" ) {
        AMP_ASSERT(args.size()==0);
        #if defined(USE_EXT_PETSC) && defined(USE_EXT_TRILINOS)
            factory.reset( new SimplePetscNativeFactory( ) );
        #else
            AMP_ERROR("Generator is not valid without support for PETSc and Trilinos");
        #endif
    } else if ( factoryName == "MultiVectorFactory" ) {
        AMP_ASSERT(args.size()==4);
        factory.reset( new MultiVectorFactory( generateVectorFactory(args[0]), to_int(args[1]),
                                               generateVectorFactory(args[2]), to_int(args[3]) ) );
    } else if ( factoryName == "NativePetscVector" ) {
        AMP_ERROR("Not Finished");
    } else if ( factoryName == "SimpleManagedVectorFactory" ) {
        AMP_ASSERT(args.size()==1);
        if ( args[0] == "ManagedPetscVector" ) {
            #if defined(USE_EXT_PETSC) && defined(USE_EXT_TRILINOS)
                factory.reset( new SimpleManagedVectorFactory<ManagedPetscVector>() );
            #else
                AMP_ERROR("Generator is not valid without support for PETSc and Trilinos");
            #endif
        } else if ( args[0] == "ManagedEpetraVector" ) {
            #ifdef USE_EXT_TRILINOS
                factory.reset( new SimpleManagedVectorFactory<ManagedEpetraVector>() );
            #else
                AMP_ERROR("Generator is not valid without support for Petsc");
            #endif
        } else {
            AMP_ERROR("Unknown template argument for SimpleManagedVectorFactory");
        }
    } else if ( factoryName == "NativeThyraFactory" ) {
        AMP_ASSERT(args.size()==0);
        #ifdef USE_TRILINOS_THYRA
            factory.reset( new NativeThyraFactory() );
        #else
            AMP_ERROR("Generator is not valid without support for Thyra");
        #endif
    } else if ( factoryName == "ManagedThyraFactory" ) {
        AMP_ASSERT(args.size()==1);
        #ifdef USE_TRILINOS_THYRA
            factory.reset( new ManagedThyraFactory( generateVectorFactory(args[0]) ) );
        #else
            AMP_ERROR("Generator is not valid without support for Thyra");
        #endif
    } else if ( factoryName == "ManagedNativeThyraFactory" ) {
        AMP_ASSERT(args.size()==1);
        #ifdef USE_TRILINOS_THYRA
            factory.reset( new ManagedNativeThyraFactory( generateVectorFactory(args[0]) ) );
        #else
            AMP_ERROR("Generator is not valid without support for Thyra");
        #endif
    } else if ( factoryName == "ViewFactory" ) {
        AMP_ASSERT(args.size()==2);
        auto factory2 = generateVectorFactory(args[1]);
        if ( args[0] == "PetscVector" ) {
            #ifdef USE_EXT_PETSC
                factory.reset( new ViewFactory<PetscVector>( factory2 ) );
            #else
                AMP_ERROR("Generator is not valid without support for Petsc");
            #endif
        } else {
            AMP_ERROR("Unknown template argument for SimpleManagedVectorFactory");
        }
    } else if ( factoryName == "CloneFactory" ) {
        AMP_ASSERT(args.size()==1);
        factory.reset( new CloneFactory( generateVectorFactory(args[0]) ) );
    } else if ( factoryName == "StridedVectorFactory" ) {
        AMP_ASSERT(args.size()==1);
        factory.reset( new CloneFactory( generateVectorFactory(args[0]) ) );
    } else {
        AMP_ERROR("Unknown factory");
    }
    return factory;
}

}
}
