#ifdef USE_AMP_DISCRETIZATION

#include "VectorBuilder.h"

#include "vectors/petsc/ManagedPetscVector.h"
#include "vectors/trilinos/EpetraVectorEngine.h"


namespace AMP {
namespace LinearAlgebra {


/********************************************************
* Vector builder                                        *
********************************************************/
AMP::LinearAlgebra::Vector::shared_ptr  createVector( 
    AMP::Discretization::DOFManager::shared_ptr DOFs, 
    AMP::LinearAlgebra::Variable::shared_ptr variable )
{
    // Create the communication list
    AMP_MPI comm = DOFs->getComm();
    AMP::LinearAlgebra::CommunicationList::shared_ptr comm_list;
    std::vector<size_t> remote_DOFs = DOFs->getRemoteDOFs();
    bool ghosts = comm.maxReduce<char>(remote_DOFs.size()>0)==1;
    if ( !ghosts ) {
        // No need for a communication list
        comm_list = AMP::LinearAlgebra::CommunicationList::createEmpty( DOFs->numLocalDOF(), DOFs->getComm() );
    } else {
        // Construct the communication list
        AMP::LinearAlgebra::CommunicationListParameters::shared_ptr params( new AMP::LinearAlgebra::CommunicationListParameters );
        params->d_comm = comm;
        params->d_localsize = DOFs->numLocalDOF();
        params->d_remote_DOFs = remote_DOFs;
        comm_list = AMP::LinearAlgebra::CommunicationList::shared_ptr( new AMP::LinearAlgebra::CommunicationList(params) );
    }
    // Create the vector parameters
    boost::shared_ptr<AMP::LinearAlgebra::ManagedPetscVectorParameters> mvparams(
        new AMP::LinearAlgebra::ManagedPetscVectorParameters() );
    boost::shared_ptr<AMP::LinearAlgebra::EpetraVectorEngineParameters> eveparams(
        new AMP::LinearAlgebra::EpetraVectorEngineParameters( DOFs->numLocalDOF(), DOFs->numGlobalDOF(), DOFs->getComm() ) );
    int i = 0;
    for (size_t local_start=DOFs->beginDOF(); local_start<DOFs->endDOF(); local_start++, i++ ) {
        eveparams->addMapping( i, local_start );
    }
    AMP::LinearAlgebra::VectorEngine::BufferPtr t_buffer ( new AMP::LinearAlgebra::VectorEngine::Buffer( DOFs->numLocalDOF() ) );
    AMP::LinearAlgebra::VectorEngine::shared_ptr epetra_engine( new AMP::LinearAlgebra::EpetraVectorEngine( eveparams, t_buffer ) );
    mvparams->d_Engine = epetra_engine;
    mvparams->d_CommList = comm_list;
    mvparams->d_DOFManager = DOFs;
    // Create the vector
    AMP::LinearAlgebra::Vector::shared_ptr vector = AMP::LinearAlgebra::Vector::shared_ptr( new AMP::LinearAlgebra::ManagedPetscVector(mvparams) );
    vector->setVariable(variable);
    return vector;
}


}
}

#endif

