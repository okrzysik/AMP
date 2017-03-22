#include "vectors/testHelpers/ThyraVectorFactory.h"

#include "discretization/DOF_Manager.h"
#include "utils/AMP_MPI.h"
#include "vectors/VectorBuilder.h"
#include "vectors/trilinos/thyra/ManagedThyraVector.h"
#include "vectors/trilinos/thyra/NativeThyraVector.h"


// Trilinos includes
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
//#include "Thyra_SpmdVectorBase_def.hpp"
#include "BelosMVOPTester.hpp"
#include "BelosThyraAdapter.hpp"
#include "Thyra_DefaultSpmdVector_def.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#ifdef USE_EXT_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

/// \cond UNDOCUMENTED


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
* NativeThyraFactory                                            *
****************************************************************/
AMP::LinearAlgebra::Variable::shared_ptr NativeThyraFactory::getVariable() const
{
    return AMP::LinearAlgebra::Variable::shared_ptr(
        new AMP::LinearAlgebra::Variable( "thyra" ) );
}
AMP::LinearAlgebra::Vector::shared_ptr NativeThyraFactory::getVector() const
{
    AMP_MPI global_comm( AMP_COMM_WORLD );
    int local_size  = 101;
    int global_size = global_comm.sumReduce( local_size );
// Create an epetra vector
#ifdef USE_EXT_MPI
    Epetra_MpiComm comm = global_comm.getCommunicator();
#else
    Epetra_SerialComm comm;
#endif
    Teuchos::RCP<Epetra_Map> epetra_map( new Epetra_Map( global_size, local_size, 0, comm ) );
    Teuchos::RCP<Epetra_Vector> epetra_v( new Epetra_Vector( *epetra_map, true ) );
    // Create a thyra vector from the epetra vector
    Teuchos::RCP<const Thyra::VectorSpaceBase<double>> space =
        Thyra::create_VectorSpace( epetra_map );
    Teuchos::RCP<Thyra::VectorBase<double>> thyra_v = Thyra::create_Vector( epetra_v, space );
    // Create the NativeThyraVector
    AMP::shared_ptr<AMP::LinearAlgebra::NativeThyraVectorParameters> params(
        new AMP::LinearAlgebra::NativeThyraVectorParameters() );
    params->d_InVec = thyra_v;
    params->d_local = local_size;
    params->d_comm  = global_comm;
    params->d_var   = getVariable();
    AMP::shared_ptr<AMP::LinearAlgebra::NativeThyraVector> vec(
        new AMP::LinearAlgebra::NativeThyraVector( params ) );
    return vec;
}
AMP::Discretization::DOFManager::shared_ptr NativeThyraFactory::getDOFMap() const
{
    AMP_ERROR("Not finished");
    return AMP::Discretization::DOFManager::shared_ptr();
}


/****************************************************************
* ManagedThyraFactory                                           *
****************************************************************/
AMP::LinearAlgebra::Variable::shared_ptr ManagedThyraFactory::getVariable() const
{
    return AMP::LinearAlgebra::Variable::shared_ptr(
        new AMP::LinearAlgebra::Variable( "managed_thyra" ) );
}
AMP::LinearAlgebra::Vector::shared_ptr ManagedThyraFactory::getVector() const
{
    // Create an arbitrary vector
    AMP::LinearAlgebra::Vector::shared_ptr vec1 = d_factory->getVector();
    vec1->setVariable( getVariable() );
    // Create the managed vector
    AMP::LinearAlgebra::Vector::shared_ptr vec2 = AMP::LinearAlgebra::ThyraVector::view( vec1 );
    vec2->setVariable( getVariable() );
    return vec2;
}
AMP::Discretization::DOFManager::shared_ptr ManagedThyraFactory::getDOFMap() const
{
    return d_factory->getDOFMap();
}


/****************************************************************
* ManagedNativeThyraFactory                                     *
****************************************************************/
AMP::LinearAlgebra::Variable::shared_ptr ManagedNativeThyraFactory::getVariable() const
{
    return AMP::LinearAlgebra::Variable::shared_ptr(
        new AMP::LinearAlgebra::Variable( "managed_native_thyra" ) );
}
AMP::LinearAlgebra::Vector::shared_ptr ManagedNativeThyraFactory::getVector() const
{
    // Create an arbitrary vector
    AMP::LinearAlgebra::Vector::shared_ptr vec1 = d_factory->getVector();
    // Create the managed vector
    AMP::shared_ptr<AMP::LinearAlgebra::ManagedThyraVector> vec2 =
        AMP::dynamic_pointer_cast<AMP::LinearAlgebra::ManagedThyraVector>(
            AMP::LinearAlgebra::ThyraVector::view( vec1 ) );
    // Create a native ThyraVector from the managed vector
    AMP::shared_ptr<AMP::LinearAlgebra::NativeThyraVectorParameters> params(
        new AMP::LinearAlgebra::NativeThyraVectorParameters() );
    params->d_InVec = vec2->getVec();
    params->d_local = vec2->getLocalSize();
    params->d_comm  = vec2->getComm();
    params->d_var   = getVariable();
    AMP::shared_ptr<AMP::LinearAlgebra::NativeThyraVector> vec3(
        new AMP::LinearAlgebra::NativeThyraVector( params ) );
    return vec3;
}
AMP::Discretization::DOFManager::shared_ptr ManagedNativeThyraFactory::getDOFMap() const
{
    return d_factory->getDOFMap();
}


/****************************************************************
* testBelosThyraVector                                          *
****************************************************************/
void testBelosThyraVector( AMP::UnitTest &ut, const VectorFactory& factory )
{
    AMP::shared_ptr<AMP::LinearAlgebra::ThyraVector> vector =
        AMP::dynamic_pointer_cast<AMP::LinearAlgebra::ThyraVector>( factory.getVector() );
    typedef Thyra::MultiVectorBase<double> TMVB;
    Teuchos::RCP<Belos::OutputManager<double>> outputmgr =
        Teuchos::rcp( new Belos::OutputManager<double>() );
    bool pass = Belos::TestMultiVecTraits<double, TMVB>( outputmgr, vector->getVec() );
    if ( pass )
        ut.passes( "Belos::TestMultiVecTraits of thyra vector" );
    else
        ut.failure( "Belos::TestMultiVecTraits of thyra vector" );
}




}
}

/// \endcond
