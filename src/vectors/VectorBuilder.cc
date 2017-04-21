#ifdef USE_AMP_DISCRETIZATION

#include "VectorBuilder.h"
#include "discretization/MultiDOF_Manager.h"
#include "vectors/ManagedVector.h"
#include "vectors/MultiVariable.h"
#include "vectors/MultiVector.h"
#include "vectors/SimpleVector.h"
#ifdef USE_EXT_PETSC
#include "vectors/petsc/ManagedPetscVector.h"
#include "vectors/petsc/PetscVector.h"
#endif
#ifdef USE_EXT_TRILINOS
#include "vectors/trilinos/epetra/EpetraVector.h"
#include "vectors/trilinos/epetra/EpetraVectorEngine.h"
#include "vectors/trilinos/ManagedEpetraVector.h"
#endif

#include <iostream>


namespace AMP {
namespace LinearAlgebra {


/********************************************************
* Vector builder                                        *
********************************************************/
AMP::LinearAlgebra::Vector::shared_ptr
createVector( AMP::Discretization::DOFManager::shared_ptr DOFs,
              AMP::LinearAlgebra::Variable::shared_ptr variable,
              bool split )
{
    if ( DOFs.get() == nullptr )
        return AMP::LinearAlgebra::Vector::shared_ptr();
    AMP_ASSERT( variable.get() != nullptr );
    // Check if we are dealing with a multiDOFManager
    AMP::shared_ptr<AMP::Discretization::multiDOFManager> multiDOF;
    if ( split )
        multiDOF = AMP::dynamic_pointer_cast<AMP::Discretization::multiDOFManager>( DOFs );
    // Check if we are dealing with a multiVariable
    AMP::shared_ptr<AMP::LinearAlgebra::MultiVariable> multiVariable =
        AMP::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVariable>( variable );
    if ( multiVariable.get() != nullptr ) {
        // We are dealing with a MultiVariable, first check that there are no duplicate or null
        // variables
        for ( size_t i = 0; i < multiVariable->numVariables(); i++ ) {
            AMP::LinearAlgebra::Variable::shared_ptr var1 = multiVariable->getVariable( i );
            AMP_INSIST( var1.get() != nullptr,
                        "Error using a MultiVariable in createVector, NULL variables detected" );
            for ( size_t j = 0; j < i; j++ ) {
                AMP::LinearAlgebra::Variable::shared_ptr var2 = multiVariable->getVariable( j );
                AMP_INSIST(
                    ( *var1 ) != ( *var2 ),
                    "Error using a MultiVariable in createVector, duplicate variables detected" );
            }
        }
        // Check that all processors have the same number of variables
        size_t N_var  = multiVariable->numVariables();
        size_t N_var0 = DOFs->getComm().bcast( N_var, 0 );
        AMP_INSIST(
            N_var == N_var0,
            "The multivariable has a different number of varaibles on different processors" );
        // Create the Vector for each variable, then combine
        std::vector<AMP::LinearAlgebra::Vector::shared_ptr> vectors;
        for ( auto it = multiVariable->beginVariable(); it != multiVariable->endVariable(); ++it )
            vectors.push_back( createVector( DOFs, *it, split ) );
        // Create the multivector
        AMP_MPI comm = DOFs->getComm();
        AMP_ASSERT( !comm.isNull() );
        comm.barrier();
        AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> multiVector =
            AMP::LinearAlgebra::MultiVector::create( variable, comm );
        multiVector->addVector( vectors );
        return multiVector;
    } else if ( multiDOF.get() != nullptr ) {
        // We are dealing with a multiDOFManager and want to split the vector based on the DOF
        // managers
        std::vector<AMP::Discretization::DOFManager::shared_ptr> subDOFs =
            multiDOF->getDOFManagers();
        // Get the vectors for each DOF manager
        std::vector<AMP::LinearAlgebra::Vector::shared_ptr> vectors( subDOFs.size() );
        for ( size_t i = 0; i < subDOFs.size(); i++ )
            vectors[i] = createVector( subDOFs[i], variable, split );
        // Create the multivector
        AMP_MPI comm = DOFs->getComm();
        AMP_ASSERT( !comm.isNull() );
        comm.barrier();
        AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> multiVector =
            AMP::LinearAlgebra::MultiVector::create( variable, comm );
        multiVector->addVector( vectors );
        return multiVector;
    } else {
        // We are ready to create a single vector
        // Create the communication list
        AMP_MPI comm = DOFs->getComm();
        AMP_ASSERT( !comm.isNull() );
        comm.barrier();
        AMP::LinearAlgebra::CommunicationList::shared_ptr comm_list;
        std::vector<size_t> remote_DOFs = DOFs->getRemoteDOFs();
        bool ghosts                     = comm.anyReduce( !remote_DOFs.empty() );
        if ( !ghosts ) {
            // No need for a communication list
            comm_list = AMP::LinearAlgebra::CommunicationList::createEmpty( DOFs->numLocalDOF(),
                                                                            DOFs->getComm() );
        } else {
            // Construct the communication list
            AMP::LinearAlgebra::CommunicationListParameters::shared_ptr params(
                new AMP::LinearAlgebra::CommunicationListParameters );
            params->d_comm        = comm;
            params->d_localsize   = DOFs->numLocalDOF();
            params->d_remote_DOFs = remote_DOFs;
            comm_list             = AMP::LinearAlgebra::CommunicationList::shared_ptr(
                new AMP::LinearAlgebra::CommunicationList( params ) );
        }
        comm.barrier();
// Create the vector parameters
#if defined( USE_EXT_PETSC ) && defined( USE_EXT_TRILINOS )
        AMP::shared_ptr<AMP::LinearAlgebra::ManagedPetscVectorParameters> mvparams(
            new AMP::LinearAlgebra::ManagedPetscVectorParameters() );
        AMP::shared_ptr<AMP::LinearAlgebra::EpetraVectorEngineParameters> eveparams(
            new AMP::LinearAlgebra::EpetraVectorEngineParameters(
                DOFs->numLocalDOF(), DOFs->numGlobalDOF(), DOFs->getComm() ) );
        comm.barrier();
        AMP::LinearAlgebra::VectorEngine::BufferPtr t_buffer(
            new AMP::LinearAlgebra::VectorEngine::Buffer( DOFs->numLocalDOF() ) );
        AMP_ASSERT( t_buffer->size() == DOFs->numLocalDOF() );
        AMP::LinearAlgebra::VectorEngine::shared_ptr epetra_engine(
            new AMP::LinearAlgebra::EpetraVectorEngine( eveparams, t_buffer ) );
        mvparams->d_Engine     = epetra_engine;
        mvparams->d_Buffer     = t_buffer;
        mvparams->d_CommList   = comm_list;
        mvparams->d_DOFManager = DOFs;
        // Create the vector
        comm.barrier();
        AMP::LinearAlgebra::Vector::shared_ptr vector = AMP::LinearAlgebra::Vector::shared_ptr(
            new AMP::LinearAlgebra::ManagedPetscVector( mvparams ) );
        vector->setVariable( variable );
        comm.barrier();
        return vector;
#elif defined( USE_EXT_TRILINOS )
        AMP::shared_ptr<AMP::LinearAlgebra::ManagedVectorParameters> mvparams(
            new AMP::LinearAlgebra::ManagedVectorParameters() );
        AMP::shared_ptr<AMP::LinearAlgebra::EpetraVectorEngineParameters> eveparams(
            new AMP::LinearAlgebra::EpetraVectorEngineParameters(
                DOFs->numLocalDOF(), DOFs->numGlobalDOF(), DOFs->getComm() ) );
        comm.barrier();
        AMP::LinearAlgebra::VectorEngine::BufferPtr t_buffer(
            new AMP::LinearAlgebra::VectorEngine::Buffer( DOFs->numLocalDOF() ) );
        AMP_ASSERT( t_buffer->size() == DOFs->numLocalDOF() );
        AMP::LinearAlgebra::VectorEngine::shared_ptr epetra_engine(
            new AMP::LinearAlgebra::EpetraVectorEngine( eveparams, t_buffer ) );
        mvparams->d_Engine     = epetra_engine;
        mvparams->d_Buffer     = t_buffer;
        mvparams->d_CommList   = comm_list;
        mvparams->d_DOFManager = DOFs;
        // Create the vector
        comm.barrier();
        AMP::LinearAlgebra::Vector::shared_ptr vector = AMP::LinearAlgebra::Vector::shared_ptr(
            new AMP::LinearAlgebra::ManagedEpetraVector( mvparams ) );
        vector->setVariable( variable );
        comm.barrier();
        return vector;
#else
        AMP::LinearAlgebra::Vector::shared_ptr vector =
            AMP::LinearAlgebra::SimpleVector<double>::create( variable, DOFs, comm_list );
        return vector;
#endif
    }
    return AMP::LinearAlgebra::Vector::shared_ptr();
}
}
}

#endif
