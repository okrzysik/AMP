#include "AMP/vectors/VectorBuilder.h"
#include "AMP/discretization/MultiDOF_Manager.h"
#include "AMP/vectors/MultiVariable.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/VectorBuilder.hpp"
#include "AMP/utils/memory.h"
#include <iostream>


namespace AMP::LinearAlgebra {


/********************************************************
 * Vector builder                                        *
 ********************************************************/
Vector::shared_ptr createVector( std::shared_ptr<AMP::Discretization::DOFManager> DOFs,
                                 std::shared_ptr<Variable> variable,
                                 bool split )
{
    if ( !DOFs )
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
            comm_list = std::make_shared<CommunicationList>( DOFs->numLocalDOF(), DOFs->getComm() );
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
 *  Explicit instantiations                              *
 ********************************************************/
#define INSTANTIATE_VECTOR( TYPE )                                                                 \
    template Vector::shared_ptr createSimpleVector<TYPE>( size_t, const std::string & );           \
    template Vector::shared_ptr createSimpleVector<TYPE>( size_t, std::shared_ptr<Variable> );     \
    template Vector::shared_ptr createSimpleVector<TYPE>(                                          \
        size_t, std::shared_ptr<Variable>, AMP_MPI );                                              \
    template Vector::shared_ptr createSimpleVector<TYPE>(                                          \
        std::shared_ptr<Variable>,                                                                 \
        std::shared_ptr<AMP::Discretization::DOFManager>,                                          \
        std::shared_ptr<CommunicationList> );                                                      \
    template Vector::shared_ptr createArrayVector<TYPE>( const ArraySize &, const std::string & ); \
    template Vector::shared_ptr createArrayVector<TYPE>( const ArraySize &,                        \
                                                         std::shared_ptr<Variable> );              \
    template Vector::shared_ptr createArrayVector<TYPE>(                                           \
        const ArraySize &, const ArraySize &, const AMP_MPI &, std::shared_ptr<Variable> );        \
    template Vector::shared_ptr createVectorAdaptor<TYPE>(                                         \
        const std::string &, std::shared_ptr<AMP::Discretization::DOFManager>, TYPE * );
INSTANTIATE_VECTOR( double )
INSTANTIATE_VECTOR( float )

#ifdef USE_CUDA
template Vector::shared_ptr
    createSimpleVector<double,
                       VectorOperationsDefault<double>,
                       VectorDataDefault<double, AMP::ManagedAllocator<double>>>(
        std::shared_ptr<Variable>,
        std::shared_ptr<AMP::Discretization::DOFManager>,
        std::shared_ptr<CommunicationList> );
#endif
#ifdef USE_HIP
template Vector::shared_ptr
    createSimpleVector<double,
                       VectorOperationsDefault<double>,
                       VectorDataDefault<double, AMP::ManagedAllocator<double>>>(
        std::shared_ptr<Variable>,
        std::shared_ptr<AMP::Discretization::DOFManager>,
        std::shared_ptr<CommunicationList> );
#endif

} // namespace AMP::LinearAlgebra
