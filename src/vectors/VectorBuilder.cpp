#include "AMP/vectors/VectorBuilder.h"
#include "AMP/TPLs.h"
#include "AMP/discretization/MultiDOF_Manager.h"
#include "AMP/vectors/MultiVariable.h"
#include "AMP/vectors/MultiVector.h"

#ifdef USE_EXT_PETSC
    #include "AMP/vectors/petsc/NativePetscVectorData.h"
    #include "AMP/vectors/petsc/NativePetscVectorOperations.h"
    #include "AMP/vectors/petsc/PetscVector.h"
    #include "petscvec.h"
#endif
#ifdef USE_EXT_TRILINOS
    #include "AMP/vectors/trilinos/epetra/EpetraVector.h"
    #include "AMP/vectors/trilinos/epetra/EpetraVectorData.h"
    #include "AMP/vectors/trilinos/epetra/EpetraVectorOperations.h"
    #include "AMP/vectors/trilinos/thyra/NativeThyraVectorData.h"
    #include "AMP/vectors/trilinos/thyra/NativeThyraVectorOperations.h"
#endif

#include <iostream>


namespace AMP::LinearAlgebra {


/********************************************************
 * Vector builder                                        *
 ********************************************************/
Vector::shared_ptr createVector( std::shared_ptr<AMP::Discretization::DOFManager> DOFs,
                                 std::shared_ptr<Variable> variable,
                                 bool split )
{
    if ( DOFs.get() == nullptr )
        return Vector::shared_ptr();
    AMP_ASSERT( variable );
    // Check if we are dealing with a multiDOFManager
    std::shared_ptr<AMP::Discretization::multiDOFManager> multiDOF;
    if ( split )
        multiDOF = std::dynamic_pointer_cast<AMP::Discretization::multiDOFManager>( DOFs );
    // Check if we are dealing with a multiVariable
    auto multiVariable = std::dynamic_pointer_cast<MultiVariable>( variable );
    if ( multiVariable ) {
        // We are dealing with a MultiVariable, first check that there are no duplicate or null
        // variables
        for ( size_t i = 0; i < multiVariable->numVariables(); i++ ) {
            auto var1 = multiVariable->getVariable( i );
            AMP_INSIST( var1,
                        "Error using a MultiVariable in createVector, NULL variables detected" );
            for ( size_t j = 0; j < i; j++ ) {
                auto var2 = multiVariable->getVariable( j );
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
        std::vector<Vector::shared_ptr> vectors;
        for ( auto var : *multiVariable )
            vectors.push_back( createVector( DOFs, var, split ) );
        // Create the multivector
        AMP_MPI comm = DOFs->getComm();
        AMP_ASSERT( !comm.isNull() );
        comm.barrier();
        auto multiVector = MultiVector::create( variable, comm );
        multiVector->addVector( vectors );
        return multiVector;
    } else if ( multiDOF ) {
        // We are dealing with a multiDOFManager and want to split the vector based on the DOF
        // managers
        auto subDOFs = multiDOF->getDOFManagers();
        // Get the vectors for each DOF manager
        std::vector<Vector::shared_ptr> vectors( subDOFs.size() );
        for ( size_t i = 0; i < subDOFs.size(); i++ )
            vectors[i] = createVector( subDOFs[i], variable, split );
        // Create the multivector
        AMP_MPI comm = DOFs->getComm();
        AMP_ASSERT( !comm.isNull() );
        comm.barrier();
        auto multiVector = MultiVector::create( variable, comm );
        multiVector->addVector( vectors );
        return multiVector;
    } else {
        // We are ready to create a single vector
        // Create the communication list
        AMP_MPI comm = DOFs->getComm();
        AMP_ASSERT( !comm.isNull() );
        comm.barrier();
        std::shared_ptr<CommunicationList> comm_list;
        auto remote_DOFs = DOFs->getRemoteDOFs();
        bool ghosts      = comm.anyReduce( !remote_DOFs.empty() );
        if ( !ghosts ) {
            // No need for a communication list
            comm_list = CommunicationList::createEmpty( DOFs->numLocalDOF(), DOFs->getComm() );
        } else {
            // Construct the communication list
            auto params           = std::make_shared<CommunicationListParameters>();
            params->d_comm        = comm;
            params->d_localsize   = DOFs->numLocalDOF();
            params->d_remote_DOFs = remote_DOFs;
            comm_list             = std::make_shared<CommunicationList>( params );
        }
        comm.barrier();
        // Create the vector
        auto vector = createSimpleVector<double>( variable, DOFs, comm_list );
        return vector;
    }
    return Vector::shared_ptr();
}


/********************************************************
 * create vector from PETSc Vec                          *
 ********************************************************/
#if defined( USE_EXT_PETSC )
std::shared_ptr<Vector>
createVector( Vec v, bool deleteable, AMP_MPI comm, std::shared_ptr<Variable> var )
{
    if ( !var )
        var = std::make_shared<Variable>( "vec" );
    auto ops  = std::make_shared<NativePetscVectorOperations>();
    auto data = std::make_shared<NativePetscVectorData>( v, deleteable, comm );
    return std::make_shared<Vector>( data, ops, var, nullptr );
}
#endif


/********************************************************
 * create vector from Trilinos Thyra vector              *
 ********************************************************/
#if defined( USE_EXT_TRILINOS ) && defined( USE_TRILINOS_THYRA )
std::shared_ptr<Vector> createVector( Teuchos::RCP<Thyra::VectorBase<double>> vec,
                                      size_t local,
                                      AMP_MPI comm,
                                      std::shared_ptr<Variable> var )
{
    if ( !var )
        var = std::make_shared<Variable>( "vec" );
    auto ops  = std::make_shared<NativeThyraVectorOperations>();
    auto data = std::make_shared<NativeThyraVectorData>( vec, local, comm );
    return std::make_shared<Vector>( data, ops, var, nullptr );
}
#endif


/********************************************************
 * create Trilinos Epetra vector                         *
 ********************************************************/
#if defined( USE_EXT_TRILINOS ) && defined( USE_TRILINOS_EPETRA )
std::shared_ptr<Vector> createEpetraVector( std::shared_ptr<CommunicationList> commList,
                                            std::shared_ptr<AMP::Discretization::DOFManager> DOFs,
                                            std::shared_ptr<VectorData> buf )
{
    auto var    = std::make_shared<Variable>( "vec" );
    auto ops    = std::make_shared<EpetraVectorOperations>();
    auto params = std::make_shared<EpetraVectorEngineParameters>(
        DOFs->numLocalDOF(), DOFs->getComm(), commList );
    auto data = EpetraVectorData::create( params, buf );
    return std::make_shared<Vector>( data, ops, var, DOFs );
}
#endif


} // namespace AMP::LinearAlgebra
