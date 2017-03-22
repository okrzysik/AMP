#include <vectors/Vector.h>
#include <vectors/testHelpers/VectorTests.h>
#include <vectors/testHelpers/VectorFactory.h>
#include <vectors/testHelpers/generateVectorFactories.h>
#include <vectors/testHelpers/generateVectorFactories.h>

#include <utils/UnitTest.h>

#include "test_ArrayVector.h"

#ifdef USE_TRILINOS_THYRA
#include <vectors/testHelpers/ThyraVectorFactory.h>
#endif


using namespace AMP::unit_test;
using namespace AMP::LinearAlgebra;
using AMP::LinearAlgebra::VectorTests;
using AMP::LinearAlgebra::generateVectorFactory;


// Define some factories
const std::string SimpleFactory1 = "SimpleVectorFactory<15,false,double>";
const std::string SimpleFactory2 = "SimpleVectorFactory<45, true,double>";
const std::string ArrayFactory1 = "ArrayVectorFactory<4,10,false,double>";
const std::string ArrayFactory2 = "ArrayVectorFactory<4,10, true,double>";
const std::string SNPVFactory = "SimplePetscNativeFactory";
const std::string SMEVFactory = "SimpleManagedVectorFactory<ManagedEpetraVector>";
#ifdef USE_EXT_PETSC
const std::string MVFactory1 = "MultiVectorFactory<"+SMEVFactory+", 1, "+SNPVFactory+", 1>";
const std::string MVFactory2 = "MultiVectorFactory<"+SMEVFactory+", 3, "+SNPVFactory+", 2>";
const std::string MVFactory3 = "MultiVectorFactory<"+MVFactory1+", 2, "+MVFactory2+", 2>";
#else
const std::string MVFactory1 = "MultiVectorFactory<SimpleVectorFactory<15,false>,1,"+SMEVFactory+",1>";
const std::string MVFactory2 = "MultiVectorFactory<SimpleVectorFactory<15,false>,3,"+SMEVFactory+",2>";
const std::string MVFactory3 = "MultiVectorFactory<"+MVFactory1+", 2, "+MVFactory2+", 2>";
#endif
#if defined(USE_EXT_PETSC) && defined(USE_EXT_TRILINOS)
const std::string ViewSMEVFactory = "ViewFactory<PetscVector,SimpleManagedVectorFactory<ManagedEpetraVector>>";
const std::string ViewSNPVFactory = "ViewFactory<PetscVector,SimplePetscNativeFactory>";
const std::string ViewMVFactory1  = "ViewFactory<PetscVector,MultiVectorFactory<"+ViewSMEVFactory+",1,"+ViewSNPVFactory+",1>>";
const std::string ViewMVFactory2  = "ViewFactory<PetscVector,MultiVectorFactory<"+ViewSMEVFactory+",3,"+ViewSNPVFactory+",2>>";
const std::string ViewMVFactory3  = "ViewFactory<PetscVector,MultiVectorFactory<"+ViewMVFactory1+",2,"+ViewMVFactory2+",2>>";
#endif
#if defined(USE_EXT_PETSC) && defined(USE_EXT_TRILINOS)
const std::string CloneSMEVFactory = "CloneFactory<SimpleManagedVectorFactory<ManagedEpetraVector>>";
const std::string CloneSNPVFactory = "CloneFactory<SimplePetscNativeFactory>";
const std::string CloneMVFactory1  = "CloneFactory<MultiVectorFactory<"+CloneSMEVFactory+",1,"+CloneSNPVFactory+",1>>";
const std::string CloneMVFactory2  = "CloneFactory<MultiVectorFactory<"+CloneSMEVFactory+",3,"+CloneSNPVFactory+",2>>";
const std::string CloneMVFactory3  = "CloneFactory<MultiVectorFactory<"+CloneMVFactory1+",2,"+CloneMVFactory2+",2>>";
#endif
#if defined(USE_EXT_PETSC) && defined(USE_EXT_TRILINOS)
const std::string CloneViewSMEVFactory = "CloneFactory<ViewFactory<PetscVector,SimpleManagedVectorFactory<ManagedEpetraVector>>>";
const std::string CloneViewSNPVFactory = "CloneFactory<ViewFactory<PetscVector,SimplePetscNativeFactory>>";
const std::string CloneViewMVFactory1  = "CloneFactory<ViewFactory<PetscVector,MultiVectorFactory<"+CloneViewSMEVFactory+",1,"+CloneViewSNPVFactory+",1>>>";
const std::string CloneViewMVFactory2  = "CloneFactory<ViewFactory<PetscVector,MultiVectorFactory<"+CloneViewSMEVFactory+",3,"+CloneViewSNPVFactory+",2>>>";
const std::string CloneViewMVFactory3  = "CloneFactory<ViewFactory<PetscVector,MultiVectorFactory<"+CloneViewMVFactory1+",2,"+CloneViewMVFactory2+",2>>>";
#endif
#ifdef USE_TRILINOS_THYRA
const std::string NTFactory = "NativeThyraFactory";
const std::string MTFactory = "ManagedThyraFactory<SimpleVectorFactory<45,true,double>>";
const std::string MNTFactory = "ManagedNativeThyraFactory<SimpleVectorFactory<45,true,double>>";
const std::string MNT_MVFactory = "ManagedNativeThyraFactory<"+MVFactory1+">";
#endif
const std::string StridedFactory = "StridedVectorFactory<"+SMEVFactory+">";



// Set the tests
void testNullVector( AMP::UnitTest &ut, const std::string& factoryName )
{
    auto factory = generateVectorFactory( factoryName );
    VectorTests tests( factory );
    tests.testNullVector( &ut );
    AMP::AMP_MPI( AMP_COMM_WORLD ).barrier();
}
void testBasicVector( AMP::UnitTest &ut, const std::string& factoryName )
{
    auto factory = generateVectorFactory( factoryName );
    VectorTests tests( factory );
    tests.testBasicVector( &ut );
    AMP::AMP_MPI( AMP_COMM_WORLD ).barrier();
}
void testManagedVector( AMP::UnitTest &ut, const std::string& factoryName )
{
    auto factory = generateVectorFactory( factoryName );
    VectorTests tests( factory );
    tests.testManagedVector( &ut );
    AMP::AMP_MPI( AMP_COMM_WORLD ).barrier();
}
void VectorIteratorTests( AMP::UnitTest &ut, const std::string& factoryName )
{
    auto factory = generateVectorFactory( factoryName );
    VectorTests tests( factory );
    tests.VectorIteratorTests( &ut );
    AMP::AMP_MPI( AMP_COMM_WORLD ).barrier();
}
void testVectorSelector( AMP::UnitTest &ut, const std::string& factoryName )
{
    auto factory = generateVectorFactory( factoryName );
    VectorTests tests( factory );
    tests.testVectorSelector( &ut );
    AMP::AMP_MPI( AMP_COMM_WORLD ).barrier();
}


