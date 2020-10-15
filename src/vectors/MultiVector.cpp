#include "AMP/vectors/MultiVector.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/MultiDOF_Manager.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/ManagedVector.h"
#include "AMP/vectors/MultiVariable.h"
#include "AMP/vectors/VectorSelector.h"
#include "AMP/vectors/data/MultiVectorData.h"
#include "AMP/vectors/operations/MultiVectorOperations.h"

#include "ProfilerApp.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>


namespace AMP {
namespace LinearAlgebra {

/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
MultiVector::MultiVector( const std::string &name, const AMP_MPI &comm ) : Vector()
{
    d_VectorOps  = std::make_shared<MultiVectorOperations>();
    d_VectorData = std::make_shared<MultiVectorData>( comm );
    d_pVariable.reset( new MultiVariable( name ) );
}

std::shared_ptr<MultiVector> MultiVector::create( Variable::shared_ptr variable,
                                                  const AMP_MPI &comm,
                                                  const std::vector<Vector::shared_ptr> &vecs )
{
    std::shared_ptr<MultiVector> retval( new MultiVector( variable->getName(), comm ) );
    retval->addVector( vecs );
    return retval;
}
std::shared_ptr<MultiVector> MultiVector::create( const std::string &name,
                                                  const AMP_MPI &comm,
                                                  const std::vector<Vector::shared_ptr> &vecs )
{
    auto variable = std::make_shared<Variable>( name );
    return MultiVector::create( variable, comm, vecs );
}
std::shared_ptr<const MultiVector>
MultiVector::const_create( Variable::shared_ptr variable,
                           const AMP_MPI &comm,
                           const std::vector<Vector::const_shared_ptr> &vecs )
{
    std::vector<Vector::shared_ptr> vecs2( vecs.size() );
    for ( size_t i = 0; i < vecs.size(); i++ )
        vecs2[i] = std::const_pointer_cast<Vector>( vecs[i] );
    return MultiVector::create( variable, comm, vecs2 );
}
std::shared_ptr<const MultiVector>
MultiVector::const_create( const std::string &name,
                           const AMP_MPI &comm,
                           const std::vector<Vector::const_shared_ptr> &vecs )
{
    auto variable = std::make_shared<Variable>( name );
    std::vector<Vector::shared_ptr> vecs2( vecs.size() );
    for ( size_t i = 0; i < vecs.size(); i++ )
        vecs2[i] = std::const_pointer_cast<Vector>( vecs[i] );
    return MultiVector::create( variable, comm, vecs2 );
}
std::shared_ptr<MultiVector> MultiVector::encapsulate( Vector::shared_ptr vec,
                                                       const AMP_MPI &comm_in )
{
    auto comm     = comm_in;
    auto multivec = std::dynamic_pointer_cast<MultiVector>( vec );
    if ( multivec ) {
        if ( !comm.isNull() )
            AMP_ASSERT( comm.compare( vec->getComm() ) != 0 );
        return multivec;
    }
    if ( comm.isNull() )
        comm = vec->getComm();
    multivec  = create( vec->getVariable()->getName(), comm, { vec } );
    auto data = multivec->Vector::getVectorData();
    AMP_ASSERT( data );
    auto listener = std::dynamic_pointer_cast<DataChangeListener>( data );
    vec->getVectorData()->registerListener( listener );
    return multivec;
}
std::shared_ptr<MultiVector> MultiVector::view( Vector::shared_ptr vec, const AMP_MPI &comm_in )
{
    auto comm = comm_in;
    // Check to see if this is a multivector
    auto multivec = std::dynamic_pointer_cast<MultiVector>( vec );
    if ( multivec ) {
        AMP_ASSERT( comm.isNull() || comm.compare( vec->getComm() ) != 0 );
        return multivec;
    }
    // Check to see if the managed vector engine is a multivector
    auto managed = std::dynamic_pointer_cast<ManagedVector>( vec );
    if ( managed )
        multivec = std::dynamic_pointer_cast<MultiVector>( managed->getVectorEngine() );
    // If still don't have a multivector, make one
    if ( !multivec ) {
        if ( comm.isNull() )
            comm = vec->getComm();
        multivec = create( vec->getVariable()->getName(), comm, { vec } );
    }
    return multivec;
}
std::shared_ptr<const MultiVector> MultiVector::constView( Vector::const_shared_ptr vec,
                                                           const AMP_MPI &comm_in )
{
    auto comm = comm_in;
    // Check to see if this is a multivector
    auto multivec = std::dynamic_pointer_cast<const MultiVector>( vec );
    if ( multivec ) {
        AMP_ASSERT( comm.isNull() || comm.compare( vec->getComm() ) != 0 );
        return multivec;
    }
    // Check to see if the manged vector engine is a multivector
    auto managed = std::dynamic_pointer_cast<const ManagedVector>( vec );
    if ( managed )
        multivec = std::dynamic_pointer_cast<const MultiVector>( managed->getVectorEngine() );
    // If still don't have a multivector, make one
    if ( !multivec ) {
        if ( comm.isNull() )
            comm = vec->getComm();
        multivec = const_create( vec->getVariable()->getName(), comm, { vec } );
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
    // Set the vector data and operations
    resetVectorData();
    resetVectorOperations();
}
void MultiVector::resetVectorOperations()
{
    std::vector<std::shared_ptr<VectorOperations>> operations( d_vVectors.size() );
    for ( size_t i = 0; i < d_vVectors.size(); i++ )
        operations[i] = d_vVectors[i]->getVectorOperations();
    AMP_ASSERT( d_VectorOps );
    auto mvOps = std::dynamic_pointer_cast<MultiVectorOperations>( d_VectorOps );
    AMP_ASSERT( mvOps );
    mvOps->resetVectorOperations( operations );
}
void MultiVector::resetVectorData()
{
    // Create a new multiDOFManager for the multivector
    std::vector<AMP::Discretization::DOFManager::shared_ptr> managers( d_vVectors.size() );
    for ( size_t i = 0; i < d_vVectors.size(); i++ ) {
        AMP_ASSERT( d_vVectors[i].get() != nullptr );
        managers[i] = d_vVectors[i]->getDOFManager();
        AMP_INSIST( managers[i].get() != nullptr,
                    "All vectors must have a DOFManager for MultiVector to work properly" );
    }
    d_DOFManager = std::make_shared<AMP::Discretization::multiDOFManager>( getComm(), managers );

    auto data   = Vector::getVectorData();
    auto mvData = std::dynamic_pointer_cast<MultiVectorData>( data );
    AMP_ASSERT( mvData );
    std::vector<VectorData *> dataComponents( d_vVectors.size() );
    for ( size_t i = 0; i < d_vVectors.size(); i++ )
        dataComponents[i] = d_vVectors[i]->Vector::getVectorData().get();
    mvData->resetMultiVectorData( d_DOFManager.get(), dataComponents );
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
        auto multivec = std::dynamic_pointer_cast<MultiVector>( vec );
        if ( !multivec ) {
            auto managed = std::dynamic_pointer_cast<ManagedVector>( vec );
            if ( managed )
                multivec = std::dynamic_pointer_cast<MultiVector>( managed->getVectorEngine() );
        }
        if ( multivec.get() != nullptr ) {
            for ( size_t i = 0; i != multivec->getNumberOfSubvectors(); i++ )
                addVectorHelper( multivec->getVector( i ) );
        } else {
            AMP_ERROR( "Not finished" );
            d_vVectors.push_back( vec );
            auto data = Vector::getVectorData();
            AMP_ASSERT( data );
            auto listener = std::dynamic_pointer_cast<DataChangeListener>( data );
            vec->getVectorData()->registerListener( listener );
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
            auto data = Vector::getVectorData();
            AMP_ASSERT( data );
            auto listener = std::dynamic_pointer_cast<DataChangeListener>( data );
            vec->getVectorData()->registerListener( listener );
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
    auto multiVar = std::dynamic_pointer_cast<MultiVariable>( d_pVariable );
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
    if ( pos >= 0 )
        d_vVectors[pos] = newVec;
    resetVectorData();
    resetVectorOperations();
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
        auto retVec2 = d_vVectors[i]->selectInto( s );
        AMP_ASSERT( std::dynamic_pointer_cast<MultiVector>( retVec2 ) == nullptr );
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
        auto retVec2 = d_vVectors[i]->selectInto( s );
        AMP_ASSERT( std::dynamic_pointer_cast<const MultiVector>( retVec2 ) == nullptr );
        if ( retVec2 != nullptr ) {
            if ( retVec2->getDOFManager()->numGlobalDOF() > 0 )
                subvectors.push_back( retVec2 );
        }
    }
    // Add the subsets to the multivector
    auto comm   = s.communicator( shared_from_this() );
    auto retVec = MultiVector::const_create( "tmp_vector", comm, subvectors );
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
        auto multivariable = std::dynamic_pointer_cast<const MultiVariable>( name );
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
    std::shared_ptr<MultiVector> retVal;
    if ( N_procs == comm.getSize() ) {
        // All processor have a variable
        retVal = create( variable, getComm() );
        retVal->addVector( subvectors );
    } else {
        // Only a subset of processors have a variable
        auto new_comm = comm.split( subvectors.empty() ? -1 : 0, comm.getRank() );
        retVal        = create( variable, new_comm );
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
void MultiVector::swapVectors( Vector &other )
{
    for ( size_t i = 0; i != d_vVectors.size(); i++ )
        d_vVectors[i]->swapVectors( getVector( other, i ) );
}
Vector::shared_ptr MultiVector::cloneVector( const Variable::shared_ptr name ) const
{
    std::shared_ptr<MultiVector> retVec( new MultiVector( name->getName(), getComm() ) );
    retVec->d_DOFManager = d_DOFManager;
    retVec->setCommunicationList( getCommunicationList() );
    retVec->d_vVectors.resize( d_vVectors.size() );
    for ( size_t i = 0; i != d_vVectors.size(); i++ )
        retVec->d_vVectors[i] = d_vVectors[i]->cloneVector();
    retVec->resetVectorData();
    retVec->resetVectorOperations();
    return retVec;
}

Vector::shared_ptr MultiVector::getVector( size_t i ) { return d_vVectors[i]; }

Vector::const_shared_ptr MultiVector::getVector( size_t i ) const { return d_vVectors[i]; }

size_t MultiVector::getNumberOfSubvectors() const { return d_vVectors.size(); }

std::string MultiVector::type() const { return "MultiVector"; }

MultiVector::~MultiVector() {}

const Vector::shared_ptr &MultiVector::getVector( const Vector &rhs, size_t which ) const
{
    auto x = dynamic_cast<const MultiVector *>( &rhs );
    AMP_ASSERT( x != nullptr );
    AMP_ASSERT( which < x->d_vVectors.size() );
    return x->d_vVectors[which];
}

Vector::shared_ptr &MultiVector::getVector( Vector &rhs, size_t which ) const
{
    auto x = dynamic_cast<MultiVector *>( &rhs );
    AMP_ASSERT( x != nullptr );
    AMP_ASSERT( which < x->d_vVectors.size() );
    return x->d_vVectors[which];
}


} // namespace LinearAlgebra
} // namespace AMP
