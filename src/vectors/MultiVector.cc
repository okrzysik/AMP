#include "vectors/MultiVector.h"
#include "discretization/DOF_Manager.h"
#include "discretization/MultiDOF_Manager.h"
#include "utils/Utilities.h"
#include "vectors/ManagedVector.h"
#include "vectors/MultiVariable.h"

#include "ProfilerApp.h"

#include <algorithm>
#include <math.h>
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
    AMP::shared_ptr<MultiVector> retval( new MultiVector( variable->getName() ) );
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
        vecs2[i]   = AMP::const_pointer_cast<Vector>( vecs[i] );
    return MultiVector::create( variable, comm, vecs2 );
}
AMP::shared_ptr<const MultiVector> MultiVector::const_create(
    const std::string &name, AMP_MPI comm, const std::vector<Vector::const_shared_ptr> &vecs )
{
    Variable::shared_ptr variable( new Variable( name ) );
    std::vector<Vector::shared_ptr> vecs2( vecs.size() );
    for ( size_t i = 0; i < vecs.size(); i++ )
        vecs2[i]   = AMP::const_pointer_cast<Vector>( vecs[i] );
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
        d_CommList            = AMP::LinearAlgebra::CommunicationList::shared_ptr(
            new AMP::LinearAlgebra::CommunicationList( params ) );
    }
    // Set the vector operations
    updateVectorOperations();
}
void MultiVector::updateVectorOperations()
{
    d_operations.resize( d_vVectors.size() );
    for ( size_t i      = 0; i < d_vVectors.size(); i++ )
        d_operations[i] = d_vVectors[i].get();
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
    AMP::shared_ptr<MultiVariable> multiVar =
        AMP::dynamic_pointer_cast<MultiVariable>( d_pVariable );
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
    AMP_MPI comm                        = s.communicator( shared_from_this() );
    AMP::shared_ptr<MultiVector> retVec = MultiVector::create( "tmp_vector", comm, subvectors );
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
    AMP::shared_ptr<const MultiVector> retVec =
        MultiVector::const_create( "tmp_vector", comm, subvectors );
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
* makeConsistent                                                *
****************************************************************/
void MultiVector::makeConsistent( ScatterType t )
{
    for ( size_t i = 0; i != d_vVectors.size(); i++ )
        d_vVectors[i]->makeConsistent( t );
    *d_UpdateState = Vector::UpdateState::UNCHANGED;
}


/****************************************************************
* query basic info                                              *
****************************************************************/
size_t MultiVector::getLocalSize() const
{
    size_t ans = 0;
    for ( size_t i = 0; i != d_vVectors.size(); i++ )
        ans += d_vVectors[i]->getLocalSize();
    AMP_ASSERT( ans == d_DOFManager->numLocalDOF() );
    return ans;
}
size_t MultiVector::getGlobalSize() const { return d_DOFManager->numGlobalDOF(); }
size_t MultiVector::getGhostSize() const
{
    size_t ans = 0;
    for ( size_t i = 0; i != d_vVectors.size(); i++ )
        ans += d_vVectors[i]->getGhostSize();
    return ans;
}


/****************************************************************
* Functions to get access to the data                           *
****************************************************************/
void MultiVector::putRawData( const double *in )
{
    int cur_off = 0;
    for ( size_t i = 0; i != d_vVectors.size(); i++ ) {
        d_vVectors[i]->putRawData( in + cur_off );
        cur_off += d_vVectors[i]->getLocalSize();
    }
}
size_t MultiVector::numberOfDataBlocks() const
{
    size_t ans = 0;
    for ( size_t i = 0; i != d_vVectors.size(); i++ )
        ans += d_vVectors[i]->numberOfDataBlocks();
    return ans;
}
size_t MultiVector::sizeOfDataBlock( size_t i ) const
{
    size_t retVal = 0;
    size_t rightOffset, leftOffset;
    rightOffset = leftOffset = 0;
    for ( size_t j = 0; j != d_vVectors.size(); j++ ) {
        rightOffset += d_vVectors[j]->numberOfDataBlocks();
        if ( i < rightOffset ) {
            retVal = d_vVectors[j]->sizeOfDataBlock( i - leftOffset );
            break;
        }
        leftOffset = rightOffset;
    }
    return retVal;
}
void MultiVector::copyOutRawData( double *out ) const
{
    size_t curOffset = 0;
    for ( size_t j = 0; j != d_vVectors.size(); j++ ) {
        d_vVectors[j]->copyOutRawData( &out[curOffset] );
        curOffset += d_vVectors[j]->getLocalSize();
    }
}
void *MultiVector::getRawDataBlockAsVoid( size_t i )
{
    size_t curOffset = 0;
    for ( size_t j = 0; j != d_vVectors.size(); j++ ) {
        curOffset += d_vVectors[j]->numberOfDataBlocks();
        if ( i < curOffset ) {
            size_t index = i + d_vVectors[j]->numberOfDataBlocks() - curOffset;
            return d_vVectors[j]->getRawDataBlock<double>( index );
        }
    }
    return nullptr;
}
const void *MultiVector::getRawDataBlockAsVoid( size_t i ) const
{
    size_t curOffset = 0;
    for ( size_t j = 0; j != d_vVectors.size(); j++ ) {
        curOffset += d_vVectors[j]->numberOfDataBlocks();
        if ( i < curOffset ) {
            size_t index = i + d_vVectors[j]->numberOfDataBlocks() - curOffset;
            return d_vVectors[j]->getRawDataBlock<double>( index );
        }
    }
    return nullptr;
}
size_t MultiVector::sizeofDataBlockType( size_t block ) const
{
    size_t curOffset = 0;
    for ( size_t j = 0; j != d_vVectors.size(); j++ ) {
        curOffset += d_vVectors[j]->numberOfDataBlocks();
        if ( block < curOffset ) {
            size_t index = block + d_vVectors[j]->numberOfDataBlocks() - curOffset;
            return d_vVectors[j]->sizeofDataBlockType( index );
        }
    }
    return 0;
}
bool MultiVector::isTypeId( size_t hash, size_t block ) const
{
    size_t curOffset = 0;
    for ( size_t j = 0; j != d_vVectors.size(); j++ ) {
        curOffset += d_vVectors[j]->numberOfDataBlocks();
        if ( block < curOffset ) {
            size_t index = block + d_vVectors[j]->numberOfDataBlocks() - curOffset;
            return d_vVectors[j]->isTypeId( hash, index );
        }
    }
    return false;
}
const void *MultiVector::getDataBlock( size_t i ) const { return getRawDataBlockAsVoid( i ); }
void *MultiVector::getDataBlock( size_t i ) { return getRawDataBlockAsVoid( i ); }


/****************************************************************
* Functions to print the data                                   *
****************************************************************/
void MultiVector::dumpOwnedData( std::ostream &out, size_t GIDoffset, size_t LIDoffset ) const
{
    size_t localOffset = 0;
    AMP::Discretization::multiDOFManager *manager =
        (AMP::Discretization::multiDOFManager *) d_DOFManager.get();
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
    AMP::Discretization::multiDOFManager *manager =
        (AMP::Discretization::multiDOFManager *) d_DOFManager.get();
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
    int N_procs = comm.sumReduce<int>( subvectors.empty() ? 0 : 1 );
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
    MultiVector *tmp = const_cast<MultiVector *>( this );
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


Vector::UpdateState MultiVector::getUpdateStatus() const
{
    Vector::UpdateState state = *d_UpdateState;
    for ( size_t i = 0; i != d_vVectors.size(); i++ ) {
        Vector::UpdateState sub_state = d_vVectors[i]->getUpdateStatus();
        if ( sub_state == UpdateState::UNCHANGED ) {
            continue;
        } else if ( sub_state == UpdateState::LOCAL_CHANGED && state == UpdateState::UNCHANGED ) {
            state = UpdateState::LOCAL_CHANGED;
        } else if ( sub_state == UpdateState::LOCAL_CHANGED ) {
            continue;
        } else if ( sub_state == UpdateState::ADDING &&
                    ( state == UpdateState::UNCHANGED || state == UpdateState::LOCAL_CHANGED ||
                      state == UpdateState::ADDING ) ) {
            state = UpdateState::ADDING;
        } else if ( sub_state == UpdateState::SETTING &&
                    ( state == UpdateState::UNCHANGED || state == UpdateState::LOCAL_CHANGED ||
                      state == UpdateState::SETTING ) ) {
            state = UpdateState::SETTING;
        } else {
            state = UpdateState::MIXED;
        }
    }
    return state;
}


void MultiVector::setUpdateStatus( UpdateState state )
{
    *d_UpdateState = state;
    for ( size_t i = 0; i != d_vVectors.size(); i++ )
        d_vVectors[i]->setUpdateStatus( state );
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


VectorEngine::BufferPtr MultiVector::getNewBuffer() { return VectorEngine::BufferPtr(); }


bool MultiVector::sameEngine( VectorEngine &rhs ) const
{
    return dynamic_cast<MultiVector *>( &rhs ) != nullptr;
}


void MultiVector::swapEngines( VectorEngine::shared_ptr p )
{
    auto vec = dynamic_pointer_cast<MultiVector>( p );
    AMP_INSIST( vec != nullptr, "Cannot swap with a non-MulitVector" );
    return swapVectors( *vec );
}


VectorEngine::shared_ptr MultiVector::cloneEngine( VectorEngine::BufferPtr ) const
{
    return AMP::dynamic_pointer_cast<VectorEngine>( Vector::cloneVector( "engine_clone" ) );
}


Vector::shared_ptr MultiVector::cloneVector( const Variable::shared_ptr name ) const
{
    AMP::shared_ptr<MultiVector> retVec( new MultiVector( name->getName() ) );
    retVec->d_Comm       = d_Comm;
    retVec->d_DOFManager = d_DOFManager;
    retVec->d_CommList   = d_CommList;
    retVec->d_vVectors.resize( d_vVectors.size() );
    for ( size_t i            = 0; i != d_vVectors.size(); i++ )
        retVec->d_vVectors[i] = d_vVectors[i]->cloneVector();
    retVec->updateVectorOperations();
    return retVec;
}


/****************************************************************
* Functions to access data by ID                                *
****************************************************************/
void MultiVector::setValuesByLocalID( int num, size_t *indices, const double *in_vals )
{
    if ( num == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<double>> vals;
    partitionLocalValues( num, indices, in_vals, ndxs, vals );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_vVectors[i]->setValuesByLocalID( ndxs[i].size(), &( ndxs[i][0] ), &( vals[i][0] ) );
    }
}
void MultiVector::setLocalValuesByGlobalID( int num, size_t *indices, const double *in_vals )
{
    if ( num == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<double>> vals;
    partitionGlobalValues( num, indices, in_vals, ndxs, vals );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_vVectors[i]->setLocalValuesByGlobalID(
                ndxs[i].size(), &( ndxs[i][0] ), &( vals[i][0] ) );
    }
}
void MultiVector::setGhostValuesByGlobalID( int num, size_t *indices, const double *in_vals )
{
    if ( num == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<double>> vals;
    partitionGlobalValues( num, indices, in_vals, ndxs, vals );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_vVectors[i]->setGhostValuesByGlobalID(
                ndxs[i].size(), &( ndxs[i][0] ), &( vals[i][0] ) );
    }
}
void MultiVector::setValuesByGlobalID( int num, size_t *indices, const double *in_vals )
{
    if ( num == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<double>> vals;
    partitionGlobalValues( num, indices, in_vals, ndxs, vals );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_vVectors[i]->setValuesByGlobalID( ndxs[i].size(), &( ndxs[i][0] ), &( vals[i][0] ) );
    }
}
void MultiVector::addValuesByLocalID( int num, size_t *indices, const double *in_vals )
{
    if ( num == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<double>> vals;
    partitionLocalValues( num, indices, in_vals, ndxs, vals );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_vVectors[i]->addValuesByLocalID( ndxs[i].size(), &( ndxs[i][0] ), &( vals[i][0] ) );
    }
}
void MultiVector::addLocalValuesByGlobalID( int num, size_t *indices, const double *in_vals )
{
    if ( num == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<double>> vals;
    partitionGlobalValues( num, indices, in_vals, ndxs, vals );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_vVectors[i]->addLocalValuesByGlobalID(
                ndxs[i].size(), &( ndxs[i][0] ), &( vals[i][0] ) );
    }
}
void MultiVector::addValuesByGlobalID( int num, size_t *indices, const double *in_vals )
{
    if ( num == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<double>> vals;
    partitionGlobalValues( num, indices, in_vals, ndxs, vals );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_vVectors[i]->addValuesByGlobalID( ndxs[i].size(), &( ndxs[i][0] ), &( vals[i][0] ) );
    }
}
void MultiVector::getValuesByGlobalID( int num, size_t *indices, double *out_vals ) const
{
    if ( num == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<double>> vals;
    std::vector<std::vector<int>> remap;
    partitionGlobalValues( num, indices, out_vals, ndxs, vals, &remap );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_vVectors[i]->getValuesByGlobalID( ndxs[i].size(), &( ndxs[i][0] ), &( vals[i][0] ) );
    }
    for ( size_t i = 0; i != remap.size(); i++ ) {
        for ( size_t j            = 0; j != remap[i].size(); j++ )
            out_vals[remap[i][j]] = vals[i][j];
    }
}
void MultiVector::getLocalValuesByGlobalID( int num, size_t *indices, double *out_vals ) const
{
    if ( num == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<double>> vals;
    std::vector<std::vector<int>> remap;
    partitionGlobalValues( num, indices, out_vals, ndxs, vals, &remap );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_vVectors[i]->getLocalValuesByGlobalID(
                ndxs[i].size(), &( ndxs[i][0] ), &( vals[i][0] ) );
    }
    for ( size_t i = 0; i != remap.size(); i++ ) {
        for ( size_t j            = 0; j != remap[i].size(); j++ )
            out_vals[remap[i][j]] = vals[i][j];
    }
}
void MultiVector::getGhostValuesByGlobalID( int num, size_t *indices, double *out_vals ) const
{
    if ( num == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<double>> vals;
    std::vector<std::vector<int>> remap;
    partitionGlobalValues( num, indices, out_vals, ndxs, vals, &remap );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_vVectors[i]->getGhostValuesByGlobalID(
                ndxs[i].size(), &( ndxs[i][0] ), &( vals[i][0] ) );
    }
    for ( size_t i = 0; i != remap.size(); i++ ) {
        for ( size_t j            = 0; j != remap[i].size(); j++ )
            out_vals[remap[i][j]] = vals[i][j];
    }
}
void MultiVector::getValuesByLocalID( int num, size_t *indices, double *out_vals ) const
{
    if ( num == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<double>> vals;
    std::vector<std::vector<int>> remap;
    partitionLocalValues( num, indices, out_vals, ndxs, vals, &remap );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_vVectors[i]->getValuesByLocalID( ndxs[i].size(), &( ndxs[i][0] ), &( vals[i][0] ) );
    }
    for ( size_t i = 0; i != remap.size(); i++ ) {
        for ( size_t j            = 0; j != remap[i].size(); j++ )
            out_vals[remap[i][j]] = vals[i][j];
    }
}


/****************************************************************
* Function to partition the global ids by the sub vectors       *
****************************************************************/
void MultiVector::partitionGlobalValues( const int num,
                                         const size_t *indices,
                                         const double *vals,
                                         std::vector<std::vector<size_t>> &out_indices,
                                         std::vector<std::vector<double>> &out_vals,
                                         std::vector<std::vector<int>> *remap ) const
{
    PROFILE_START( "partitionGlobalValues", 2 );
    const size_t neg_one = ~( (size_t) 0 );
    std::vector<size_t> globalDOFs( num, neg_one );
    for ( int i       = 0; i < num; i++ )
        globalDOFs[i] = indices[i];
    out_indices.resize( d_vVectors.size() );
    out_vals.resize( d_vVectors.size() );
    if ( remap != nullptr )
        remap->resize( d_vVectors.size() );
    AMP::Discretization::multiDOFManager *manager =
        (AMP::Discretization::multiDOFManager *) d_DOFManager.get();
    std::vector<AMP::Discretization::DOFManager::shared_ptr> DOFManagers =
        manager->getDOFManagers();
    AMP_ASSERT( DOFManagers.size() == d_vVectors.size() );
    for ( size_t i = 0; i < d_vVectors.size(); i++ ) {
        AMP_ASSERT( d_vVectors[i]->getDOFManager() == DOFManagers[i] );
        std::vector<size_t> subDOFs = manager->getSubDOF( i, globalDOFs );
        size_t count                = 0;
        for ( auto &subDOF : subDOFs ) {
            if ( subDOF != neg_one )
                count++;
        }
        out_indices[i] = std::vector<size_t>( count, neg_one );
        out_vals[i]    = std::vector<double>( count, 0.0 );
        if ( remap != nullptr )
            remap->operator[]( i ) = std::vector<int>( count, -1 );
        count                      = 0;
        for ( size_t j = 0; j < subDOFs.size(); j++ ) {
            if ( subDOFs[j] != neg_one ) {
                out_indices[i][count] = subDOFs[j];
                out_vals[i][count]    = vals[j];
                if ( remap != nullptr )
                    remap->operator[]( i )[count] = j;
                count++;
            }
        }
    }
    PROFILE_STOP( "partitionGlobalValues", 2 );
}


/****************************************************************
* Function to partition the local ids by the sub vectors       *
****************************************************************/
void MultiVector::partitionLocalValues( const int num,
                                        const size_t *indices,
                                        const double *vals,
                                        std::vector<std::vector<size_t>> &out_indices,
                                        std::vector<std::vector<double>> &out_vals,
                                        std::vector<std::vector<int>> *remap ) const
{
    if ( num == 0 )
        return;
    PROFILE_START( "partitionLocalValues", 2 );
    // Convert the local ids to global ids
    size_t begin_DOF = d_DOFManager->beginDOF();
    size_t end_DOF   = d_DOFManager->endDOF();
    std::vector<size_t> global_indices( num );
    for ( int i = 0; i < num; i++ ) {
        AMP_INSIST( indices[i] < end_DOF, "Invalid local id" );
        global_indices[i] = indices[i] + begin_DOF;
    }
    // Partition based on the global ids
    partitionGlobalValues( num, &global_indices[0], vals, out_indices, out_vals, remap );
    // Convert the new global ids back to local ids
    const size_t neg_one = ~( (size_t) 0 );
    for ( size_t i = 0; i < d_vVectors.size(); i++ ) {
        if ( out_indices[i].size() == 0 )
            continue;
        begin_DOF = d_vVectors[i]->getDOFManager()->beginDOF();
        end_DOF   = d_vVectors[i]->getDOFManager()->endDOF();
        for ( auto &elem : out_indices[i] ) {
            AMP_ASSERT( elem != neg_one );
            elem -= begin_DOF;
            AMP_ASSERT( elem < end_DOF );
        }
    }
    PROFILE_STOP( "partitionLocalValues", 2 );
}


} // LinearAlgebra namespace
} // AMP namespace
