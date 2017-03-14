#ifndef included_test_Vector
#define included_test_Vector

#include "utils/AMP_MPI.h"
#include "vectors/ArrayVector.h"
#include "vectors/ManagedVector.h"
#include "vectors/MultiVariable.h"
#include "vectors/NullVector.h"
#include "vectors/SimpleVector.h"
#include "vectors/Variable.h"
#ifdef USE_EXT_PETSC
#include <vectors/petsc/NativePetscVector.h>
#endif
#ifdef USE_EXT_TRILINOS
#include <vectors/trilinos/EpetraVectorEngine.h>
#ifdef USE_TRILINOS_THYRA
#include "test_ThyraVector.h"
#endif
#endif


namespace AMP {
namespace LinearAlgebra {


template <typename T>
class CloneFactory
{
public:
    typedef typename T::vector vector;

    static AMP::LinearAlgebra::Variable::shared_ptr getVariable()
    {
        return T::getVariable()->cloneVariable( "noname" );
    }

    static AMP::LinearAlgebra::Vector::shared_ptr getVector()
    {
        AMP::LinearAlgebra::Vector::shared_ptr vec = T::getVector();
        return vec->cloneVector();
    }

    static std::string name() { return "CloneFactory<" + T::name() + ">"; }
};


class NullVectorFactory
{
public:
    typedef AMP::LinearAlgebra::NullVector vector;

    static AMP::LinearAlgebra::Variable::shared_ptr getVariable()
    {
        return AMP::LinearAlgebra::Variable::shared_ptr(
            new AMP::LinearAlgebra::Variable( "null" ) );
    }

    static AMP::LinearAlgebra::Vector::shared_ptr getVector()
    {
        return AMP::LinearAlgebra::NullVector::create( "null" );
    }

    static std::string name() { return "NullVectorFactory"; }
};


template <int I, bool GLOBAL, typename T>
class SimpleVectorFactory
{
public:
    typedef AMP::LinearAlgebra::SimpleVector<T> vector;

    static AMP::LinearAlgebra::Variable::shared_ptr getVariable()
    {
        return AMP::LinearAlgebra::Variable::shared_ptr(
            new AMP::LinearAlgebra::Variable( "simple" ) );
    }

    static AMP::LinearAlgebra::Vector::shared_ptr getVector()
    {
        AMP::LinearAlgebra::Vector::shared_ptr vec;
        if ( GLOBAL )
            vec = AMP::LinearAlgebra::SimpleVector<T>::create(
                I, getVariable(), AMP_MPI( AMP_COMM_WORLD ) );
        else
            vec = AMP::LinearAlgebra::SimpleVector<T>::create( I, getVariable() );
        return vec;
    }

    static std::string name() { return "SimpleVectorFactory"; }
};


template <int D, int I, bool GLOBAL, typename T>
class ArrayVectorFactory
{
public:
    typedef AMP::LinearAlgebra::ArrayVector<T> vector;

    static AMP::LinearAlgebra::Variable::shared_ptr getVariable()
    {
        return AMP::LinearAlgebra::Variable::shared_ptr(
            new AMP::LinearAlgebra::Variable( "array" ) );
    }

    static AMP::LinearAlgebra::Vector::shared_ptr getVector()
    {
        AMP::LinearAlgebra::Vector::shared_ptr vec;
        if ( GLOBAL )
            vec = AMP::LinearAlgebra::ArrayVector<T>::create(
                std::vector<size_t>( D, I ), getVariable(), AMP_MPI( AMP_COMM_WORLD ) );
        else
            vec = AMP::LinearAlgebra::ArrayVector<T>::create( std::vector<size_t>( D, I ),
                                                              getVariable() );
        return vec;
    }

    static std::string name() { return "ArrayVectorFactory"; }
};

#ifdef USE_EXT_TRILINOS
template <typename T>
class SimpleManagedVectorFactory
{
public:
    typedef T vector;

    static AMP::LinearAlgebra::Variable::shared_ptr getVariable()
    {
        return AMP::LinearAlgebra::Variable::shared_ptr(
            new AMP::LinearAlgebra::Variable( "..." ) );
    }

    static AMP::shared_ptr<T> getVector()
    {
        const int num_local = 210;
        AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
        AMP::shared_ptr<AMP::LinearAlgebra::EpetraVectorEngineParameters> epetraParams(
            new AMP::LinearAlgebra::EpetraVectorEngineParameters(
                num_local, num_local * globalComm.getSize(), globalComm ) );
        AMP::shared_ptr<AMP::LinearAlgebra::ManagedVectorParameters> managedParams(
            new AMP::LinearAlgebra::ManagedVectorParameters );
        AMP::LinearAlgebra::VectorEngine::BufferPtr buffer( new std::vector<double>( 120 ) );
        managedParams->d_Engine = AMP::LinearAlgebra::VectorEngine::shared_ptr(
            new AMP::LinearAlgebra::EpetraVectorEngine( epetraParams, buffer ) );
        managedParams->d_CommList =
            AMP::LinearAlgebra::CommunicationList::createEmpty( 210, globalComm );
        managedParams->d_DOFManager = AMP::Discretization::DOFManager::shared_ptr(
            new AMP::Discretization::DOFManager( 210, globalComm ) );
        AMP::shared_ptr<T> retval( new T( managedParams ) );
        retval->setVariable( AMP::LinearAlgebra::Variable::shared_ptr(
            new AMP::LinearAlgebra::Variable( "Test Vector" ) ) );
        return retval;
    }

    static std::string name() { return "SimpleManagedVectorFactory"; }
};
#endif


#ifdef USE_EXT_PETSC
template <typename T>
class PetscManagedVectorFactory
{
public:
    typedef T vector;

    static AMP::LinearAlgebra::Variable::shared_ptr getVariable()
    {
        return AMP::LinearAlgebra::Variable::shared_ptr(); // no variable.....
    }

    static AMP::LinearAlgebra::Vector::shared_ptr getVector()
    {
        Vec v;
        AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
        VecCreate( globalComm.getCommunicator(), &v );
        VecSetSizes( v, 15, PETSC_DECIDE );
        AMP::shared_ptr<AMP::LinearAlgebra::NativePetscVectorParameters> npvParams(
            new AMP::LinearAlgebra::NativePetscVectorParameters( v, true ) );
        AMP::LinearAlgebra::NativePetscVector *newVec =
            new AMP::LinearAlgebra::NativePetscVector( npvParams );
        VecSetFromOptions( v );
        newVec->assemble();
        auto p1        = new AMP::LinearAlgebra::ManagedVectorParameters;
        p1->d_Engine   = AMP::LinearAlgebra::VectorEngine::shared_ptr( newVec );
        p1->d_CommList = AMP::LinearAlgebra::CommunicationList::createEmpty( 210, globalComm );
        AMP::LinearAlgebra::Vector::shared_ptr retval(
            new T( AMP::LinearAlgebra::VectorParameters::shared_ptr( p1 ) ) );
        retval->setVariable( AMP::LinearAlgebra::Variable::shared_ptr(
            new AMP::LinearAlgebra::Variable( "Test Vector" ) ) );
        return retval;
    }

    static std::string name() { return "PetscManagedVectorFactory"; }
};
#endif

template <typename TYPE1, typename FACTORY2>
class ViewFactory
{
public:
    typedef TYPE1 vector;

    static AMP::LinearAlgebra::Variable::shared_ptr getVariable()
    {
        return FACTORY2::getVariable();
    }

    static AMP::LinearAlgebra::Vector::shared_ptr getVector()
    {
        return vector::view( FACTORY2::getVector() );
    }

    static std::string name() { return "ViewFactory"; }
};


template <typename FACTORY1, int NUM1, typename FACTORY2, int NUM2>
class MultiVectorFactory
{
public:
    typedef AMP::LinearAlgebra::MultiVector vector;

    static AMP::LinearAlgebra::Variable::shared_ptr getVariable()
    {
        AMP::LinearAlgebra::MultiVariable *newVar = new AMP::LinearAlgebra::MultiVariable( "var1" );
        for ( int i = 0; i != NUM1; i++ )
            newVar->add( FACTORY1::getVariable() );
        for ( int i = 0; i != NUM2; i++ )
            newVar->add( FACTORY2::getVariable() );
        return AMP::LinearAlgebra::Variable::shared_ptr( newVar );
    }

    static AMP::LinearAlgebra::Vector::shared_ptr getVector()
    {
        AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
        AMP::LinearAlgebra::Vector::shared_ptr retVal =
            AMP::LinearAlgebra::MultiVector::create( getVariable(), globalComm );
        for ( int i = 0; i != NUM1; i++ )
            retVal->castTo<AMP::LinearAlgebra::MultiVector>().addVector( FACTORY1::getVector() );
        for ( int i = 0; i != NUM2; i++ )
            retVal->castTo<AMP::LinearAlgebra::MultiVector>().addVector( FACTORY2::getVector() );
        return retVal;
    }

    static std::string name() { return "MultiVectorFactory"; }
};
}
}

/// \endcond

#endif
