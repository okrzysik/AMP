#include "AMP/vectors/MultiVector.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/MultiDOF_Manager.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/ManagedVector.h"
#include "AMP/vectors/MultiVariable.h"

#include "ProfilerApp.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
MultiVector::MultiVector( const std::string &name )
{
    d_pVariable.reset( new MultiVariable( name ) );
    d_CommCreated = false;
}
AMP::shared_ptr<MultiVector> MultiVector::create( Variable::shared_ptr variable,
                                                  AMP_MPI comm,
                                                  const std::vector<Vector::shared_ptr> &vecs )
{
    auto retval    = AMP::shared_ptr<MultiVector>( new MultiVector( variable->getName() ) );
    retval->d_Comm = comm;
    retval->addVector( vecs );
    return retval;
}
AMP::shared_ptr<MultiVector> MultiVector::create( const std::string &name,
                                                  AMP_MPI comm,
                                                  const std::vector<Vector::shared_ptr> &vecs )
{
    Variable::shared_ptr variable( new Variable( name ) );
    return MultiVector::create( variable, comm, vecs );
}
AMP::shared_ptr<const MultiVector> MultiVector::const_create(
    Variable::shared_ptr variable, AMP_MPI comm, const std::vector<Vector::const_shared_ptr> &vecs )
{
    std::vector<Vector::shared_ptr> vecs2( vecs.size() );
    for ( size_t i = 0; i < vecs.size(); i++ )
        vecs2[i] = AMP::const_pointer_cast<Vector>( vecs[i] );
    return MultiVector::create( variable, comm, vecs2 );
}
AMP::shared_ptr<const MultiVector> MultiVector::const_create(
    const std::string &name, AMP_MPI comm, const std::vector<Vector::const_shared_ptr> &vecs )
{
    Variable::shared_ptr variable( new Variable( name ) );
    std::vector<Vector::shared_ptr> vecs2( vecs.size() );
    for ( size_t i = 0; i < vecs.size(); i++ )
        vecs2[i] = AMP::const_pointer_cast<Vector>( vecs[i] );
    return MultiVector::create( variable, comm, vecs2 );
}
AMP::shared_ptr<MultiVector> MultiVector::encapsulate( Vector::shared_ptr vec, AMP_MPI comm )
{
    auto multivec = AMP::dynamic_pointer_cast<MultiVector>( vec );
    if ( multivec ) {
        if ( !comm.isNull() )
            AMP_ASSERT( comm.compare( vec->getComm() ) != 0 );
        return multivec;
    }
    if ( comm.isNull() )
        comm = vec->getComm();
    multivec =
        create( vec->getVariable()->getName(), comm, std::vector<Vector::shared_ptr>( 1, vec ) );
    auto firer = dynamic_pointer_cast<DataChangeFirer>( vec );
    if ( firer )
        firer->registerListener( multivec.get() );
    return multivec;
}
AMP::shared_ptr<MultiVector> MultiVector::view( Vector::shared_ptr vec, AMP_MPI comm )
{
    // Check to see if this is a multivector
    auto multivec = AMP::dynamic_pointer_cast<MultiVector>( vec );
    if ( multivec ) {
        if ( !comm.isNull() )
            AMP_ASSERT( comm.compare( vec->getComm() ) != 0 );
    }
    // Check to see if the engine is a multivector
    auto managed = AMP::dynamic_pointer_cast<ManagedVector>( vec );
    if ( managed ) {
        auto vec2 = AMP::dynamic_pointer_cast<MultiVector>( managed->getVectorEngine() );
        if ( vec2 ) {
            if ( !comm.isNull() )
                AMP_ASSERT( comm.compare( vec->getComm() ) != 0 );
            multivec = vec2;
        }
    }
    // If still don't have a multivector, make one
    if ( !multivec ) {
        if ( comm.isNull() )
            comm = vec->getComm();
        multivec = create(
            vec->getVariable()->getName(), comm, std::vector<Vector::shared_ptr>( 1, vec ) );
    }
    return multivec;
}
AMP::shared_ptr<const MultiVector> MultiVector::constView( Vector::const_shared_ptr vec,
                                                           AMP_MPI comm )
{
    // Check to see if this is a multivector
    auto multivec = AMP::dynamic_pointer_cast<const MultiVector>( vec );
    if ( multivec ) {
        if ( !comm.isNull() )
            AMP_ASSERT( comm.compare( vec->getComm() ) != 0 );
    }
    // Check to see if the engine is a multivector
    auto managed = AMP::dynamic_pointer_cast<const ManagedVector>( vec );
    if ( managed ) {
        auto vec2 = AMP::dynamic_pointer_cast<const MultiVector>( managed->getVectorEngine() );
        if ( vec2 ) {
            if ( !comm.isNull() )
                AMP_ASSERT( comm.compare( vec->getComm() ) != 0 );
            multivec = vec2;
        }
    }
    // If still don't have a multivector, make one
    if ( !multivec ) {
        if ( comm.isNull() )
            comm = vec->getComm();
        multivec = const_create(
            vec->getVariable()->getName(), comm, std::vector<Vector::const_shared_ptr>( 1, vec ) );
    }
    return multivec;
}


/****************************************************************
 * Functions to add/remove vectors                               *
 ****************************************************************/
void MultiVector::addVector( Vector::shared_ptr v )
{
    this->addVector( std::vector<Vector::shared_ptr>( 1, v ) );
}
void MultiVector::addVector( std::vector<Vector::shared_ptr> v )
{
    // Add the vectors
    for ( auto &elem : v )
        addVectorHelper( elem );
    // Create a new multiDOFManager for the multivector
    std::vector<AMP::Discretization::DOFManager::shared_ptr> managers( d_vVectors.size() );
    for ( size_t i = 0; i < d_vVectors.size(); i++ ) {
        AMP_ASSERT( d_vVectors[i].get() != nullptr );
        managers[i] = d_vVectors[i]->getDOFManager();
        AMP_INSIST( managers[i].get() != nullptr,
                    "All vectors must have a DOFManager for MultiVector to work properly" );
    }
    d_DOFManager = AMP::Discretization::DOFManager::shared_ptr(
        new AMP::Discretization::multiDOFManager( d_Comm, managers ) );
    // Create a new communication list
    std::vector<size_t> remote_DOFs = d_DOFManager->getRemoteDOFs();
    bool ghosts                     = d_Comm.anyReduce( !remote_DOFs.empty() );
    if ( !ghosts ) {
        d_CommList = AMP::LinearAlgebra::CommunicationList::createEmpty(
            d_DOFManager->numLocalDOF(), d_Comm );
    } else {
        AMP::LinearAlgebra::CommunicationListParameters::shared_ptr params(
            new AMP::LinearAlgebra::CommunicationListParameters );
        params->d_comm        = d_Comm;
        params->d_localsize   = d_DOFManager->numLocalDOF();
        params->d_remote_DOFs = remote_DOFs;
        d_CommList            = AMP::make_shared<AMP::LinearAlgebra::CommunicationList>( params );
    }
    // Set the vector data and operations
    updateVectorData();
    updateVectorOperations();
}
void MultiVector::updateVectorOperations()
{
    d_operations.resize( d_vVectors.size() );
    for ( size_t i = 0; i < d_vVectors.size(); i++ )
        d_operations[i] = d_vVectors[i].get();
}
void MultiVector::updateVectorData()
{
    d_data.resize( d_vVectors.size() );
    for ( size_t i = 0; i < d_vVectors.size(); i++ )
        d_data[i] = d_vVectors[i].get();
    auto *manager = dynamic_cast<AMP::Discretization::multiDOFManager *>( d_DOFManager.get() );
    AMP_ASSERT( manager != nullptr );
    d_globalDOFManager = manager;
    auto subManagers   = manager->getDOFManagers();
    AMP_ASSERT( subManagers.size() == d_vVectors.size() );
    d_subDOFManager.resize( d_vVectors.size() );
    for ( size_t i = 0; i < d_vVectors.size(); i++ )
        d_subDOFManager[i] = subManagers[i].get();
}
void MultiVector::addVectorHelper( Vector::shared_ptr vec )
{
    if ( vec.get() == nullptr )
        return;
    auto id = vec->getDataID();
    if ( id == 0 ) {
        // We are dealing with a multivector, empty vector, or something special
        if ( vec->getGlobalSize() == 0 )
            return;
        auto multivec = AMP::dynamic_pointer_cast<MultiVector>( vec );
        if ( multivec == nullptr ) {
            auto managed = AMP::dynamic_pointer_cast<ManagedVector>( vec );
            if ( managed )
                multivec = AMP::dynamic_pointer_cast<MultiVector>( managed->getVectorEngine() );
        }
        if ( multivec.get() != nullptr ) {
            for ( size_t i = 0; i != multivec->getNumberOfSubvectors(); i++ )
                addVectorHelper( multivec->getVector( i ) );
        } else {
            AMP_ERROR( "Not finished" );
            d_vVectors.push_back( vec );
            auto firer = AMP::dynamic_pointer_cast<DataChangeFirer>( vec );
            if ( firer )
                firer->registerListener( this );
        }
    } else {
        // We are dealing with a single vector, check if it is already added in some form
        int index = -1;
        for ( size_t j = 0; j < d_vVectors.size(); j++ ) {
            if ( id == d_vVectors[j]->getDataID() )
                index = static_cast<int>( j );
        }
        if ( index == -1 ) {
            // Add the vector
            d_vVectors.push_back( vec );
            auto firer = AMP::dynamic_pointer_cast<DataChangeFirer>( vec );
            if ( firer )
                firer->registerListener( this );
        } else {
            // the vector exists, which vector (or both) do we keep?
            auto dof1 = vec->getDOFManager();
            auto dof2 = d_vVectors[index]->getDOFManager();
            if ( dof1 == dof2 ) {
                // We are dealing with the same vector, no need to do anything
            } else {
                AMP_ERROR( "Not finished" );
            }
        }
    }
    // Append the variable if we have a multivariable
    auto multiVar = AMP::dynamic_pointer_cast<MultiVariable>( d_pVariable );
    if ( multiVar != nullptr )
        multiVar->add( vec->getVariable() );
}
void MultiVector::eraseVector( Vector::shared_ptr ) { AMP_ERROR( "Needs to be fixed" ); }
void MultiVector::replaceSubVector( Vector::shared_ptr oldVec, Vector::shared_ptr newVec )
{
    AMP_ASSERT( oldVec.get() != nullptr );
    AMP_ASSERT( newVec.get() != nullptr );
    AMP_INSIST( oldVec->getDOFManager() == newVec->getDOFManager(),
                "oldVec and newVec must chare the same DOFManager" );
    int pos = -1;
    for ( size_t i = 0; i < d_vVectors.size(); ++i ) {
        if ( d_vVectors[i] == oldVec ) {
            pos = i;
            break;
        }
    } // end for i
    AMP_INSIST( pos != -1, "oldVec was not found" );
    d_vVectors[pos] = newVec;
    updateVectorData();
    updateVectorOperations();
}


/****************************************************************
 * Select into the vector                                        *
 ****************************************************************/
Vector::shared_ptr MultiVector::selectInto( const VectorSelector &s )
{
    // Subset each vector
    std::vector<Vector::shared_ptr> subvectors;
    for ( size_t i = 0; i != d_vVectors.size(); i++ ) {
        // Subset the individual vector
        Vector::shared_ptr retVec2 = d_vVectors[i]->selectInto( s );
        AMP_ASSERT( AMP::dynamic_pointer_cast<MultiVector>( retVec2 ) == nullptr );
        if ( retVec2 != nullptr ) {
            if ( retVec2->getDOFManager()->numGlobalDOF() > 0 )
                subvectors.push_back( retVec2 );
        }
    }
    // Add the subsets to the multivector
    AMP_MPI comm = s.communicator( shared_from_this() );
    auto retVec  = MultiVector::create( "tmp_vector", comm, subvectors );
    if ( retVec->getDOFManager()->numGlobalDOF() == 0 )
        retVec.reset();
    return retVec;
}
Vector::const_shared_ptr MultiVector::selectInto( const VectorSelector &s ) const
{
    // Subset each vector
    std::vector<Vector::const_shared_ptr> subvectors;
    for ( size_t i = 0; i != d_vVectors.size(); i++ ) {
        // Subset the individual vector
        Vector::const_shared_ptr retVec2 = d_vVectors[i]->selectInto( s );
        AMP_ASSERT( AMP::dynamic_pointer_cast<const MultiVector>( retVec2 ) == nullptr );
        if ( retVec2 != nullptr ) {
            if ( retVec2->getDOFManager()->numGlobalDOF() > 0 )
                subvectors.push_back( retVec2 );
        }
    }
    // Add the subsets to the multivector
    AMP_MPI comm = s.communicator( shared_from_this() );
    auto retVec  = MultiVector::const_create( "tmp_vector", comm, subvectors );
    if ( retVec->getDOFManager()->numGlobalDOF() == 0 )
        retVec.reset();
    return retVec;
}


/****************************************************************
 * Other functions                                               *
 ****************************************************************/
bool MultiVector::containsPointer( const Vector::shared_ptr p ) const
{
    for ( size_t i = 0; i != d_vVectors.size(); i++ ) {
        if ( d_vVectors[i].get() == p.get() ) {
            return true;
        }
    }
    return false;
}


/****************************************************************
 * Functions to print the data                                   *
 ****************************************************************/
void MultiVector::dumpOwnedData( std::ostream &out, size_t GIDoffset, size_t LIDoffset ) const
{
    size_t localOffset = 0;
    auto *manager      = (AMP::Discretization::multiDOFManager *) d_DOFManager.get();
    AMP_ASSERT( manager->getDOFManagers().size() == d_vVectors.size() );
    for ( size_t i = 0; i != d_vVectors.size(); i++ ) {
        if ( d_vVectors[i]->getVariable() )
            out << "[ " << d_vVectors[i]->getVariable()->getName() << " ]\n";
        AMP::Discretization::DOFManager::shared_ptr subManager = d_vVectors[i]->getDOFManager();
        std::vector<size_t> subStartDOF( 1, subManager->beginDOF() );
        std::vector<size_t> globalStartDOF = manager->getGlobalDOF( i, subStartDOF );
        size_t globalOffset                = globalStartDOF[0] - subStartDOF[0];
        d_vVectors[i]->dumpOwnedData( out, GIDoffset + globalOffset, LIDoffset + localOffset );
        localOffset += d_vVectors[i]->getLocalSize();
    }
}
void MultiVector::dumpGhostedData( std::ostream &out, size_t offset ) const
{
    auto *manager = (AMP::Discretization::multiDOFManager *) d_DOFManager.get();
    AMP_ASSERT( manager->getDOFManagers().size() == d_vVectors.size() );
    for ( size_t i = 0; i != d_vVectors.size(); i++ ) {
        if ( d_vVectors[i]->getVariable() )
            out << "[ " << d_vVectors[i]->getVariable()->getName() << " ]\n";
        AMP::Discretization::DOFManager::shared_ptr subManager = d_vVectors[i]->getDOFManager();
        std::vector<size_t> subStartDOF( 1, subManager->beginDOF() );
        std::vector<size_t> globalStartDOF = manager->getGlobalDOF( i, subStartDOF );
        size_t globalOffset                = globalStartDOF[0] - subStartDOF[0];
        d_vVectors[i]->dumpGhostedData( out, offset + globalOffset );
    }
}


/****************************************************************
 * Subset                                                        *
 ****************************************************************/
Vector::shared_ptr MultiVector::subsetVectorForVariable( Variable::const_shared_ptr name )
{
    // Subset a multivector for a variable
    /* A variable used to contain a mesh and a name, now it only contains a name
     * as a result we need to subset for the variable name (there may be many)
     * and then create a new multivector if necessary */
    AMP_ASSERT( name.get() != nullptr );

    // Check if the variable matches the variable of the multivector
    if ( *d_pVariable == *name ) {
        return shared_from_this();
    }

    // Get a list of all sub vectors matching the variable
    std::vector<Vector::shared_ptr> subvectors;
    for ( size_t i = 0; i != d_vVectors.size(); i++ ) {
        Vector::shared_ptr subset = d_vVectors[i]->subsetVectorForVariable( name );
        if ( subset.get() != nullptr )
            subvectors.push_back( subset );
    }

    // If no vectors were found, check if the variable is actually a multivariable and subset on it
    const AMP_MPI &comm = getComm();
    if ( comm.sumReduce( subvectors.size() ) == 0 ) {
        auto multivariable = AMP::dynamic_pointer_cast<const MultiVariable>( name );
        if ( multivariable.get() != nullptr ) {
            bool all_found = true;
            std::vector<Vector::shared_ptr> sub_subvectors( multivariable->numVariables() );
            for ( size_t i = 0; i != multivariable->numVariables(); i++ ) {
                Vector::shared_ptr t = subsetVectorForVariable( multivariable->getVariable( i ) );
                if ( !t ) {
                    all_found = false;
                    break;
                }
                sub_subvectors[i] = t;
            }
            if ( all_found ) {
                for ( size_t i = 0; i != multivariable->numVariables(); i++ )
                    subvectors.push_back( sub_subvectors[i] );
            }
        }
    }

    // Create the new vector
    auto N_procs = comm.sumReduce<int>( subvectors.empty() ? 0 : 1 );
    if ( N_procs == 0 )
        return Vector::shared_ptr();
    auto variable = name->cloneVariable( name->getName() );
    AMP::shared_ptr<MultiVector> retVal;
    if ( N_procs == comm.getSize() ) {
        // All processor have a variable
        retVal = create( variable, getComm() );
        retVal->addVector( subvectors );
    } else {
        // Only a subset of processors have a variable
        AMP_MPI new_comm = comm.split( subvectors.empty() ? -1 : 0, comm.getRank() );
        retVal           = create( variable, new_comm );
        retVal->addVector( subvectors );
    }
    return retVal;
}
Vector::const_shared_ptr
MultiVector::constSubsetVectorForVariable( Variable::const_shared_ptr name ) const
{
    auto *tmp = const_cast<MultiVector *>( this );
    return tmp->subsetVectorForVariable( name );
}


/****************************************************************
 * Misc functions                                                *
 ****************************************************************/
void MultiVector::assemble()
{
    for ( size_t i = 0; i != d_vVectors.size(); i++ )
        d_vVectors[i]->assemble();
}
void MultiVector::swapVectors( Vector &other )
{
    for ( size_t i = 0; i != d_vVectors.size(); i++ )
        d_vVectors[i]->swapVectors( getVector( other, i ) );
}
void MultiVector::aliasVector( Vector &other )
{
    for ( size_t i = 0; i != d_vVectors.size(); i++ )
        d_vVectors[i]->aliasVector( getVector( other, i ) );
}
AMP::shared_ptr<VectorData> MultiVector::getNewBuffer() { return AMP::shared_ptr<VectorData>(); }
bool MultiVector::sameEngine( VectorEngine &rhs ) const
{
    return dynamic_cast<MultiVector *>( &rhs ) != nullptr;
}
void MultiVector::swapEngines( AMP::shared_ptr<VectorEngine> p )
{
    auto vec = dynamic_pointer_cast<MultiVector>( p );
    AMP_INSIST( vec != nullptr, "Cannot swap with a non-MulitVector" );
    return swapVectors( *vec );
}
AMP::shared_ptr<VectorEngine> MultiVector::cloneEngine( AMP::shared_ptr<VectorData> ) const
{
    return AMP::dynamic_pointer_cast<VectorEngine>( Vector::cloneVector( "engine_clone" ) );
}
Vector::shared_ptr MultiVector::cloneVector( const Variable::shared_ptr name ) const
{
    auto retVec          = AMP::shared_ptr<MultiVector>( new MultiVector( name->getName() ) );
    retVec->d_Comm       = d_Comm;
    retVec->d_DOFManager = d_DOFManager;
    retVec->d_CommList   = d_CommList;
    retVec->d_vVectors.resize( d_vVectors.size() );
    for ( size_t i = 0; i != d_vVectors.size(); i++ )
        retVec->d_vVectors[i] = d_vVectors[i]->cloneVector();
    retVec->updateVectorData();
    retVec->updateVectorOperations();
    return retVec;
}


} // namespace LinearAlgebra
} // namespace AMP
