#include "MultiVector.h"
#include "ManagedVector.h"
#include <stdexcept>
#include <algorithm>
#include <math.h>

#include "utils/Utilities.h"
#include "ProfilerApp.h"
#include "discretization/DOF_Manager.h"
#include "discretization/MultiDOF_Manager.h"


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
* Constructors                                                  *
****************************************************************/
MultiVector::MultiVector ( const std::string& name )
{
    d_pVariable.reset( new MultiVariable( name ) );
    d_CommCreated = false;
}
AMP::shared_ptr<MultiVector>  MultiVector::create( Variable::shared_ptr variable, AMP_MPI comm, 
    const std::vector<Vector::shared_ptr>& vecs )
{
    AMP::shared_ptr<MultiVector>  retval( new MultiVector(variable->getName()) );
    retval->d_Comm = comm;
    retval->addVector(vecs);
    return retval;
}
AMP::shared_ptr<MultiVector>  MultiVector::create( const std::string &name, AMP_MPI comm,
    const std::vector<Vector::shared_ptr>& vecs )
{
    Variable::shared_ptr  variable( new Variable( name ) );
    return MultiVector::create( variable, comm, vecs );
}
AMP::shared_ptr<const MultiVector>  MultiVector::const_create( Variable::shared_ptr variable, AMP_MPI comm, 
    const std::vector<Vector::const_shared_ptr>& vecs )
{
    std::vector<Vector::shared_ptr> vecs2(vecs.size());
    for (size_t i=0; i<vecs.size(); i++)
        vecs2[i] = AMP::const_pointer_cast<Vector>(vecs[i]);
    return MultiVector::create( variable, comm, vecs2 );
}
AMP::shared_ptr<const MultiVector>  MultiVector::const_create( const std::string &name, AMP_MPI comm,
    const std::vector<Vector::const_shared_ptr>& vecs )
{
    Variable::shared_ptr  variable( new Variable( name ) );
    std::vector<Vector::shared_ptr> vecs2(vecs.size());
    for (size_t i=0; i<vecs.size(); i++)
        vecs2[i] = AMP::const_pointer_cast<Vector>(vecs[i]);
    return MultiVector::create( variable, comm, vecs2 );
}
AMP::shared_ptr<MultiVector>  MultiVector::encapsulate ( Vector::shared_ptr &vec, AMP_MPI comm )
{
    if ( vec->isA<MultiVector>() ) {
        if ( !comm.isNull() )
            AMP_ASSERT( comm.compare(vec->getComm()) != 0 );
        return AMP::dynamic_pointer_cast<MultiVector>( vec );
    }
    if ( comm.isNull() )
        comm = vec->getComm();
    AMP::shared_ptr<MultiVector>  retval = create( vec->getVariable()->getName(), comm, std::vector<Vector::shared_ptr>(1,vec) );
    if ( vec->isA<DataChangeFirer>() )
        vec->castTo<DataChangeFirer>().registerListener ( &(retval->castTo<DataChangeListener>()) );
    return retval;
}
AMP::shared_ptr<MultiVector>  MultiVector::view ( Vector::shared_ptr &vec, AMP_MPI comm )
{
    AMP::shared_ptr<MultiVector>  retval;
    // Check to see if this is a multivector
    if ( vec->isA<MultiVector>() ) {
        if ( !comm.isNull() )
            AMP_ASSERT( comm.compare(vec->getComm()) != 0 );
        retval = AMP::dynamic_pointer_cast<MultiVector>( vec );
    }
    // Check to see if the engine is a multivector
    if ( vec->isA<ManagedVector>() ) {
        if ( vec->castTo<ManagedVector>().getVectorEngine()->isA<MultiVector>() ) {
            if ( !comm.isNull() )
                AMP_ASSERT( comm.compare(vec->getComm()) != 0 );
            retval = AMP::dynamic_pointer_cast<MultiVector> ( vec->castTo<ManagedVector>().getVectorEngine() );
        }
    }
    // If still don't have a multivector, make one
    if ( !retval ) {
        if ( comm.isNull() )
            comm = vec->getComm();
        retval = create ( vec->getVariable()->getName(), comm, std::vector<Vector::shared_ptr>(1,vec) );
    }
    return retval;
}
AMP::shared_ptr<const MultiVector>  MultiVector::view ( Vector::const_shared_ptr &vec, AMP_MPI comm )
{
    AMP::shared_ptr<const MultiVector>  retval;
    // Check to see if this is a multivector
    if ( vec->isA<const MultiVector>() ) {
        if ( !comm.isNull() )
            AMP_ASSERT( comm.compare(vec->getComm()) != 0 );
        retval = AMP::dynamic_pointer_cast<const MultiVector>( vec );
    }
    // Check to see if the engine is a multivector
    if ( vec->isA<const ManagedVector>() ) {
        if ( vec->castTo<const ManagedVector>().getVectorEngine()->isA<const MultiVector>() ) {
            if ( !comm.isNull() )
                AMP_ASSERT( comm.compare(vec->getComm()) != 0 );
            retval = AMP::dynamic_pointer_cast<const MultiVector>( vec->castTo<const ManagedVector>().getVectorEngine() );
        }
    }
    // If still don't have a multivector, make one
    if ( !retval ) {
        if ( comm.isNull() )
            comm = vec->getComm();
        retval = const_create( vec->getVariable()->getName(), comm, std::vector<Vector::const_shared_ptr>(1,vec) );
    }
    return retval;
}


/****************************************************************
* Functions to add/remove vectors                               *
****************************************************************/
void MultiVector::addVector ( Vector::shared_ptr  v )
{
    this->addVector ( std::vector<Vector::shared_ptr>(1,v) );
}
void MultiVector::addVector ( std::vector<Vector::shared_ptr> v )
{
    // Add the vectors
    for (size_t i=0; i<v.size(); i++) {
        if ( v[i].get() != NULL ) {
            AMP::shared_ptr<MultiVector> vec;
            if ( v[i]->isA<MultiVector>() ) {
                vec = AMP::dynamic_pointer_cast<MultiVector>( v[i] );
            } else if ( v[i]->isA<ManagedVector>() ) {
                if ( v[i]->castTo<ManagedVector>().getVectorEngine()->isA<MultiVector>() ) {
                    vec = AMP::dynamic_pointer_cast<MultiVector> ( v[i]->castTo<ManagedVector>().getVectorEngine() );
                }
            }
            if ( vec.get() != NULL ) {
                for (size_t i = 0 ; i != vec->castTo<MultiVector>().getNumberOfSubvectors() ; i++ ) {
                    Vector::shared_ptr  curvec = vec->castTo<MultiVector>().getVector ( i );
                    d_vVectors.push_back ( curvec );
                    if ( curvec->isA<DataChangeFirer>() ) {
                        curvec->castTo<DataChangeFirer>().registerListener ( this );
                    }
                }
            } else {
                d_vVectors.push_back ( v[i] );
                if ( v[i]->isA<DataChangeFirer>() ) {
                    v[i]->castTo<DataChangeFirer>().registerListener ( this );
                }
            }
            // Append the variable if we have a multivariable
            AMP::shared_ptr<MultiVariable> multiVar = AMP::dynamic_pointer_cast<MultiVariable>(d_pVariable);
            if ( multiVar!=NULL )
                multiVar->add( v[i]->getVariable() );
        }
    }
    // Create a new multiDOFManager for the multivector
    std::vector<AMP::Discretization::DOFManager::shared_ptr> managers(d_vVectors.size());
    for (size_t i=0; i<d_vVectors.size(); i++) {
        AMP_ASSERT(d_vVectors[i].get()!=NULL);
        managers[i] = d_vVectors[i]->getDOFManager();
        AMP_INSIST(managers[i].get()!=NULL,"All vectors must have a DOFManager for MultiVector to work properly");
    }
    d_DOFManager = AMP::Discretization::DOFManager::shared_ptr( new AMP::Discretization::multiDOFManager( d_Comm, managers ) );
    // Create a new communication list
    std::vector<size_t> remote_DOFs = d_DOFManager->getRemoteDOFs();
    bool ghosts = d_Comm.anyReduce(!remote_DOFs.empty());
    if ( !ghosts ) {
         d_CommList = AMP::LinearAlgebra::CommunicationList::createEmpty( d_DOFManager->numLocalDOF(), d_Comm );
    } else {
         AMP::LinearAlgebra::CommunicationListParameters::shared_ptr params( new AMP::LinearAlgebra::CommunicationListParameters );
         params->d_comm = d_Comm;
         params->d_localsize = d_DOFManager->numLocalDOF();
         params->d_remote_DOFs = remote_DOFs;
         d_CommList = AMP::LinearAlgebra::CommunicationList::shared_ptr( new AMP::LinearAlgebra::CommunicationList(params) );
    }
}
void  MultiVector::eraseVector ( Vector::shared_ptr  v )
{
    AMP_ERROR("Needs to be fixed");
}
void MultiVector::replaceSubVector(Vector::shared_ptr oldVec, Vector::shared_ptr newVec) 
{
    AMP_ASSERT(oldVec.get()!=NULL);
    AMP_ASSERT(newVec.get()!=NULL);
    AMP_INSIST(oldVec->getDOFManager()==newVec->getDOFManager(),"oldVec and newVec must chare the same DOFManager");
    int pos = -1;
    for(size_t i = 0; i < d_vVectors.size(); ++i) {
        if(d_vVectors[i] == oldVec) {
            pos = i;
            break;
        }
    }//end for i
    AMP_INSIST(pos!=-1,"oldVec was not found");
    d_vVectors[pos] = newVec;
}


/****************************************************************
* Select into the vector                                        *
****************************************************************/
Vector::shared_ptr MultiVector::selectInto ( const VectorSelector &s )
{
    // Subset each vector
    std::vector<Vector::shared_ptr> subvectors;
    for (size_t i=0; i!=d_vVectors.size(); i++) {
        // Subset the individual vector
        Vector::shared_ptr  retVec2 = d_vVectors[i]->selectInto(s);
        AMP_ASSERT(AMP::dynamic_pointer_cast<MultiVector>(retVec2)==NULL);
        if ( retVec2 != NULL ) {
            if ( retVec2->getDOFManager()->numGlobalDOF() > 0 )
                subvectors.push_back( retVec2 );
        }
    }
    // Add the subsets to the multivector
    AMP_MPI comm = s.communicator( shared_from_this() );
    AMP::shared_ptr<MultiVector>  retVec = MultiVector::create( "tmp_vector", comm, subvectors );
    if ( retVec->getDOFManager()->numGlobalDOF() == 0 )
        retVec.reset();
    return retVec;
}
Vector::const_shared_ptr MultiVector::selectInto ( const VectorSelector &s ) const
{
    // Subset each vector
    std::vector<Vector::const_shared_ptr> subvectors;
    for (size_t i=0; i!=d_vVectors.size(); i++) {
        // Subset the individual vector
        Vector::const_shared_ptr  retVec2 = d_vVectors[i]->selectInto(s);
        AMP_ASSERT(AMP::dynamic_pointer_cast<const MultiVector>(retVec2)==NULL);
        if ( retVec2 != NULL ) {
            if ( retVec2->getDOFManager()->numGlobalDOF() > 0 )
                subvectors.push_back( retVec2 );
        }
    }
    // Add the subsets to the multivector
    AMP_MPI comm = s.communicator( shared_from_this() );
    AMP::shared_ptr<const MultiVector>  retVec = MultiVector::const_create( "tmp_vector", comm, subvectors );
    if ( retVec->getDOFManager()->numGlobalDOF() == 0 )
        retVec.reset();
    return retVec;
}


/****************************************************************
* Other functions                                               *
****************************************************************/
bool MultiVector::containsPointer ( const Vector::shared_ptr p ) const
{
    for (size_t i=0; i!=d_vVectors.size(); i++ )
    {
      if ( d_vVectors[i].get() == p.get() )
      {
        return true;
      }
    }
    return false;
}


/****************************************************************
* Basic linear algebra                                          *
****************************************************************/
void MultiVector::subtract ( const VectorOperations &x, const VectorOperations &y )
{
    if ( !x.isA<MultiVector>() || !y.isA<MultiVector>() )
        AMP_ERROR("x or y is not a multivector and this is");
    for (size_t i=0; i!=d_vVectors.size(); i++ )
        d_vVectors[i]->subtract ( getVector( x, i ), getVector ( y, i ) );
}
void MultiVector::multiply ( const VectorOperations &x, const VectorOperations &y )
{
    if ( !x.isA<MultiVector>() || !y.isA<MultiVector>() )
        AMP_ERROR("x or y is not a multivector and this is");
    for (size_t i=0; i!=d_vVectors.size(); i++ )
        d_vVectors[i]->multiply ( getVector( x, i ), getVector ( y, i ) );
}
void MultiVector::divide ( const VectorOperations &x, const VectorOperations &y )
{
    if ( !x.isA<MultiVector>() || !y.isA<MultiVector>() )
        AMP_ERROR("x or y is not a multivector and this is");
    for (size_t i=0; i!=d_vVectors.size(); i++ )
        d_vVectors[i]->divide( getVector( x, i ), getVector ( y, i ) );
}
void MultiVector::reciprocal ( const VectorOperations &x )
{
    if ( !x.isA<MultiVector>() )
        AMP_ERROR("x is not a multivector and this is");
    for (size_t i=0; i!=d_vVectors.size(); i++ )
        d_vVectors[i]->reciprocal( getVector( x, i ) );
}
void MultiVector::linearSum(double alpha, const VectorOperations &x, double beta, const VectorOperations &y)
{
    if ( !x.isA<MultiVector>() || !y.isA<MultiVector>() )
        AMP_ERROR("x or y is not a multivector and this is");
    for (size_t i=0; i!=d_vVectors.size(); i++ )
        d_vVectors[i]->linearSum( alpha, getVector( x, i ), beta, getVector ( y, i ) );
}
void MultiVector::axpy(double alpha, const VectorOperations &x, const VectorOperations &y)
{
    if ( !x.isA<MultiVector>() || !y.isA<MultiVector>() )
        AMP_ERROR("x or y is not a multivector and this is");
    for (size_t i=0; i!=d_vVectors.size(); i++ )
        d_vVectors[i]->axpy ( alpha, getVector( x, i ), getVector ( y, i ) );
}
void MultiVector::axpby(double alpha, double beta, const VectorOperations &x)
{
    if ( x.isA<MultiVector>() ) {
        // this and x are both multivectors
        for (size_t i=0; i!=d_vVectors.size(); i++ )
            d_vVectors[i]->axpby( alpha, beta, getVector( x, i ) );
    } else if ( d_vVectors.size()==1 ) {
        // x is not a multivector, but this only contains one vector, try comparing
        d_vVectors[0]->axpby( alpha, beta, x );
    } else {
        AMP_ERROR("x is not a multivector and this is");
    }
}
void MultiVector::abs ( const VectorOperations &x )
{
    if ( !x.isA<MultiVector>() )
        AMP_ERROR("x is not a multivector and this is");
    for (size_t i=0; i!=d_vVectors.size(); i++ )
        d_vVectors[i]->abs( getVector( x, i ) );
}
void MultiVector::setRandomValues ()
{
    for (size_t i=0; i!=d_vVectors.size(); i++ )
        d_vVectors[i]->setRandomValues ();
}
void MultiVector::setToScalar ( double alpha )
{
    for (size_t i=0; i!=d_vVectors.size(); i++ )
      d_vVectors[i]->setToScalar ( alpha );
}
void MultiVector::scale ( double alpha, const VectorOperations &x )
{
    if ( !x.isA<MultiVector>() )
        AMP_ERROR("x is not a multivector and this is");
    for (size_t i=0; i!=d_vVectors.size(); i++ )
      d_vVectors[i]->scale ( alpha, getVector( x, i ) );
}
void MultiVector::scale ( double alpha )
{
    for (size_t i=0; i!=d_vVectors.size(); i++ )
      d_vVectors[i]->scale ( alpha );
}
void MultiVector::add ( const VectorOperations &x, const VectorOperations &y )
{
    if ( !x.isA<MultiVector>() || !y.isA<MultiVector>() )
        AMP_ERROR("x or y is not a multivector and this is");
    for (size_t i=0; i!=d_vVectors.size(); i++ )
      d_vVectors[i]->add ( getVector( x, i ), getVector ( y, i ) );
}


/****************************************************************
* min, max, norms, etc.                                         *
* Note: these routines require communication                    *
****************************************************************/
double MultiVector::min(void) const
{
    double ans = 1e300;
    for (size_t i=0; i<d_vVectors.size(); i++)
        ans = std::min ( ans, d_vVectors[i]->localMin() );
    ans = getComm().minReduce(ans);
    return ans;
}
double MultiVector::max(void) const
{
    double ans = -1e300;
    for (size_t i=0; i<d_vVectors.size(); i++)
        ans = std::max ( ans, d_vVectors[i]->localMax() );
    ans = getComm().maxReduce(ans);
    return ans;
}
double MultiVector::L1Norm () const
{
    double ans = 0.0;
    for (size_t i=0; i<d_vVectors.size(); i++)
        ans += d_vVectors[i]->localL1Norm();
    ans = getComm().sumReduce(ans);
    return ans;
}
double MultiVector::L2Norm () const
{
    double ans = 0.0;
    for (size_t i=0; i<d_vVectors.size(); i++) {
        double tmp = d_vVectors[i]->localL2Norm();
        ans += tmp*tmp;
    }
    ans = getComm().sumReduce(ans);
    return sqrt( ans );
}
double MultiVector::maxNorm () const
{
    double ans = 0.0;
    for (size_t i=0; i<d_vVectors.size(); i++)
        ans = std::max ( ans, d_vVectors[i]->localMaxNorm() );
    ans = getComm().maxReduce(ans);
    return ans;
}
double MultiVector::dot ( const VectorOperations &x ) const
{
    if ( !x.isA<MultiVector>() )
        AMP_ERROR("x is not a multivector and this is");
    double ans = 0.0;
    for (size_t i=0; i<d_vVectors.size(); i++)
        ans += d_vVectors[i]->localDot( getVector( x, i ) );
    ans = getComm().sumReduce(ans);
    return ans;
}


/****************************************************************
* makeConsistent                                                *
****************************************************************/
void MultiVector::makeConsistent ( ScatterType t )
{
    for (size_t i=0; i!=d_vVectors.size(); i++ )
        d_vVectors[i]->makeConsistent ( t );
    *d_UpdateState = Vector::UNCHANGED;
}


/****************************************************************
* query basic info                                              *
****************************************************************/
size_t MultiVector::getLocalSize () const
{
    size_t ans = 0;
    for (size_t i=0; i!=d_vVectors.size(); i++ )
      ans += d_vVectors[i]->getLocalSize();
    AMP_ASSERT(ans==d_DOFManager->numLocalDOF());
    return ans;
}
size_t MultiVector::getGlobalSize () const
{
    return d_DOFManager->numGlobalDOF();
}
size_t MultiVector::getGhostSize () const
{
    size_t ans = 0;
    for (size_t i=0; i!=d_vVectors.size(); i++ )
      ans += d_vVectors[i]->getGhostSize();
    return ans;
}


/****************************************************************
* Functions to get access to the data                           *
****************************************************************/
void  MultiVector::putRawData ( const double *in )
{
    int cur_off = 0;
    for (size_t i=0; i!=d_vVectors.size(); i++ )
    {
      d_vVectors[i]->putRawData ( in + cur_off );
      cur_off += d_vVectors[i]->getLocalSize();
    }
}
size_t MultiVector::numberOfDataBlocks () const
{
    size_t ans = 0;
    for (size_t i=0; i!=d_vVectors.size(); i++ )
      ans += d_vVectors[i]->numberOfDataBlocks();
    return ans;
}
size_t MultiVector::sizeOfDataBlock ( size_t i ) const
{
    size_t retVal = 0;
    size_t rightOffset, leftOffset;
    rightOffset = leftOffset = 0;
    for (size_t j = 0 ; j != d_vVectors.size() ; j++ )
    {
      rightOffset += d_vVectors[j]->numberOfDataBlocks();
      if ( i < rightOffset )
      {
        retVal = d_vVectors[j]->sizeOfDataBlock ( i - leftOffset );
        break;
      }
     leftOffset = rightOffset;
    }
    return retVal;
}
void  MultiVector::copyOutRawData ( double * out ) const
{
    size_t curOffset = 0;
    for (size_t j = 0 ; j != d_vVectors.size() ; j++ )
    {
        d_vVectors[j]->copyOutRawData( &out[curOffset] );
        curOffset += d_vVectors[j]->getLocalSize();
    }
}
void * MultiVector::getRawDataBlockAsVoid ( size_t i )
{
    size_t curOffset = 0;
    for (size_t j = 0 ; j != d_vVectors.size() ; j++ )
    {
      curOffset += d_vVectors[j]->numberOfDataBlocks();
      if ( i < curOffset )
      {
        return d_vVectors[j]->getRawDataBlock<double> ( i - curOffset + d_vVectors[j]->numberOfDataBlocks () );
      }
    }
    return 0;
}
const void * MultiVector::getRawDataBlockAsVoid ( size_t i ) const
{
    size_t curOffset = 0;
    for (size_t j = 0 ; j != d_vVectors.size() ; j++ )
    {
      curOffset += d_vVectors[j]->numberOfDataBlocks();
      if ( i < curOffset )
      {
        return d_vVectors[j]->getRawDataBlock<double> ( i - curOffset + d_vVectors[j]->numberOfDataBlocks () );
      }
    }
    return 0;
}
const void* MultiVector::getDataBlock ( size_t i ) const
{
    return getRawDataBlockAsVoid ( i );
}
void* MultiVector::getDataBlock ( size_t i )
{
    return getRawDataBlockAsVoid ( i );
}


/****************************************************************
* Functions to print the data                                   *
****************************************************************/
void MultiVector::dumpOwnedData ( std::ostream &out, size_t GIDoffset, size_t LIDoffset ) const
{
    size_t localOffset = 0;
    AMP::Discretization::multiDOFManager* manager = (AMP::Discretization::multiDOFManager*) d_DOFManager.get();
    AMP_ASSERT(manager->getDOFManagers().size()==d_vVectors.size());
    for (size_t i=0; i!=d_vVectors.size(); i++ ) {
        if ( d_vVectors[i]->getVariable() )
            out << "[ " << d_vVectors[i]->getVariable()->getName() << " ]\n";
        AMP::Discretization::DOFManager::shared_ptr subManager = d_vVectors[i]->getDOFManager();
        std::vector<size_t> subStartDOF(1,subManager->beginDOF());
        std::vector<size_t> globalStartDOF = manager->getGlobalDOF(i,subStartDOF);
        size_t globalOffset = globalStartDOF[0]-subStartDOF[0];
        d_vVectors[i]->dumpOwnedData( out, GIDoffset+globalOffset, LIDoffset+localOffset );
        localOffset += d_vVectors[i]->getLocalSize();
    }
}
void MultiVector::dumpGhostedData ( std::ostream &out, size_t offset ) const
{
    AMP::Discretization::multiDOFManager* manager = (AMP::Discretization::multiDOFManager*) d_DOFManager.get();
    AMP_ASSERT(manager->getDOFManagers().size()==d_vVectors.size());
    for (size_t i=0; i!=d_vVectors.size(); i++ ) {
        if ( d_vVectors[i]->getVariable() )
            out << "[ " << d_vVectors[i]->getVariable()->getName() << " ]\n";
        AMP::Discretization::DOFManager::shared_ptr subManager = d_vVectors[i]->getDOFManager();
        std::vector<size_t> subStartDOF(1,subManager->beginDOF());
        std::vector<size_t> globalStartDOF = manager->getGlobalDOF(i,subStartDOF);
        size_t globalOffset = globalStartDOF[0]-subStartDOF[0];
        d_vVectors[i]->dumpGhostedData( out, offset+globalOffset );
    }
}


/****************************************************************
* Subset                                                        *
****************************************************************/
Vector::shared_ptr  MultiVector::subsetVectorForVariable ( const Variable::shared_ptr  &name )
{
    // Subset a multivector for a variable
    /* A variable used to contain a mesh and a name, now it only contains a name
     * as a result we need to subset for the variable name (there may be many)
     * and then create a new multivector if necessary */
    AMP_ASSERT ( name.get()!=NULL );

    // Check if the variable matches the variable of the multivector
    if ( *d_pVariable == *name ) {
        return shared_from_this ();
    }

    // Get a list of all sub vectors matching the variable
    std::vector< Vector::shared_ptr > subvectors;
    for (size_t i=0; i!=d_vVectors.size(); i++) {
        Vector::shared_ptr subset = d_vVectors[i]->subsetVectorForVariable( name );
        if ( subset.get() != NULL )
            subvectors.push_back( subset );
    }

    // If no vectors were found, check if the variable is actually a multivariable and subset on it
    const AMP_MPI &comm = getComm();
    if ( comm.sumReduce(subvectors.size())==0 ) {
        AMP::shared_ptr<MultiVariable> multivariable = AMP::dynamic_pointer_cast<MultiVariable>(name);
        if ( multivariable.get() != NULL ) {
            bool all_found = true;
            std::vector<Vector::shared_ptr> sub_subvectors(multivariable->numVariables());
            for (size_t i=0; i!=multivariable->numVariables(); i++) {
                Vector::shared_ptr  t = subsetVectorForVariable ( multivariable->getVariable ( i ) );
                if ( !t ) {
                    all_found = false;
                    break;
                }
                sub_subvectors[i] = t;
            }
            if ( all_found ) {
                for (size_t i=0; i!=multivariable->numVariables(); i++)
                    subvectors.push_back( sub_subvectors[i] );
            }
        }
    }

    // Create the new vector
    int N_procs = comm.sumReduce<int>(subvectors.empty()?0:1);
    if ( N_procs==0 )
        return Vector::shared_ptr();
    AMP::shared_ptr<MultiVector> retVal;
    if ( N_procs == comm.getSize() ) {
        // All processor have a variable
        retVal = create ( name, getComm() );
        retVal->addVector( subvectors );
    } else {
        // Only a subset of processors have a variable
        AMP_MPI new_comm = comm.split( subvectors.empty()?-1:0, comm.getRank() );
        retVal = create ( name, new_comm );
        retVal->addVector( subvectors );
    }
    return retVal;
}
Vector::const_shared_ptr  MultiVector::constSubsetVectorForVariable ( const Variable::shared_ptr  &name ) const
{
    MultiVector *tmp = const_cast<MultiVector*>(this);
    return tmp->subsetVectorForVariable(name);
}


/****************************************************************
* Misc functions                                                *
****************************************************************/
void MultiVector::assemble ()
{
    for (size_t i=0; i!=d_vVectors.size(); i++ )
      d_vVectors[i]->assemble ();
}


Vector::UpdateState  MultiVector::getUpdateStatus () const
{
    Vector::UpdateState state = *d_UpdateState;
    for (size_t i=0; i!=d_vVectors.size(); i++) {
        Vector::UpdateState  sub_state = d_vVectors[i]->getUpdateStatus();
        if ( sub_state==UNCHANGED ) {
            continue;
        } else if ( sub_state==LOCAL_CHANGED && state==UNCHANGED ) {
            state = LOCAL_CHANGED;
        } else if ( sub_state==LOCAL_CHANGED ) {
            continue;
        } else if ( sub_state==ADDING && ( state==UNCHANGED || state==LOCAL_CHANGED || state==ADDING ) ) {
            state = ADDING;
        } else if ( sub_state==SETTING && ( state==UNCHANGED || state==LOCAL_CHANGED || state==SETTING ) ) {
            state = SETTING;
        } else {
            state = MIXED;
        }
    }
    return state;
}



void MultiVector::setUpdateStatus ( UpdateState state ) 
{
    *d_UpdateState = state;
    for (size_t i=0; i!=d_vVectors.size(); i++)
        d_vVectors[i]->setUpdateStatus(state);
}


void MultiVector::copyVector( Vector::const_shared_ptr src )
{
    AMP::shared_ptr<const MultiVector> rhs = AMP::dynamic_pointer_cast<const MultiVector>(src);
    if ( rhs.get()!=NULL ) {
        // We are dealing with 2 multivectors
        AMP_ASSERT(rhs->d_vVectors.size()==d_vVectors.size());
        for (size_t i=0; i!=d_vVectors.size(); i++ )
          d_vVectors[i]->copyVector( rhs->d_vVectors[i] );
        *d_UpdateState = *(rhs->getUpdateStatusPtr());
    } else if ( d_vVectors.size()==1 ) {
        // We have a multivector of a single vector
        d_vVectors[0]->copyVector( src );
    } else if ( *getDOFManager()==*(src->getDOFManager()) ) {
        // The two DOFManagers are compatible, we can perform a basic copy
        VectorDataIterator dst_it = Vector::begin();
        for (ConstVectorDataIterator src_it=src->begin(); src_it!=src->end(); ++dst_it, ++src_it)
            *dst_it = *src_it;
    } else {
        AMP_ERROR("Unable to copy vector");
    }
}


void MultiVector::swapVectors(Vector &other)
{
    for (size_t i=0; i!=d_vVectors.size(); i++ )
      d_vVectors[i]->swapVectors ( getVector ( other, i ) );
}


void MultiVector::aliasVector (Vector &other)
{
    for (size_t i=0; i!=d_vVectors.size(); i++ )
      d_vVectors[i]->aliasVector ( getVector ( other, i ) );
}



VectorEngine::BufferPtr  MultiVector::getNewBuffer()
{
    return VectorEngine::BufferPtr ();
}


bool                     MultiVector::sameEngine ( VectorEngine &rhs ) const
{
    return rhs.isA<MultiVector>();
}


void                     MultiVector::swapEngines ( VectorEngine::shared_ptr p )
{
    return swapVectors ( p->castTo<MultiVector>() );
}


VectorEngine::shared_ptr MultiVector::cloneEngine ( VectorEngine::BufferPtr  ) const
{
    return AMP::dynamic_pointer_cast<VectorEngine> ( Vector::cloneVector ( "engine_clone" ) );
}



Vector::shared_ptr MultiVector::cloneVector(const Variable::shared_ptr name) const
{
    AMP::shared_ptr<MultiVector> retVec( new MultiVector( name->getName() ) );
    retVec->d_Comm = d_Comm;
    retVec->d_DOFManager = d_DOFManager;
    retVec->d_CommList = d_CommList;
    retVec->d_vVectors.resize ( d_vVectors.size() );
    for (size_t i=0; i!=d_vVectors.size(); i++ )
        retVec->d_vVectors[i] = d_vVectors[i]->cloneVector ();
    return retVec;
}


/****************************************************************
* Functions to access data by ID                                *
****************************************************************/
void MultiVector::setValuesByLocalID ( int num, size_t *indices, const double *in_vals )
{
    if ( num==0 )
        return;
    INCREMENT_COUNT("Virtual");
    std::vector<std::vector<size_t> >  ndxs;
    std::vector<std::vector<double> >  vals;
    partitionLocalValues( num, indices, in_vals, ndxs, vals );
    for (size_t i = 0 ; i != ndxs.size() ; i++ ) {
        if ( ndxs[i].size() )
            d_vVectors[i]->setValuesByLocalID ( ndxs[i].size(), &(ndxs[i][0]), &(vals[i][0]) );
    }
}
void MultiVector::setLocalValuesByGlobalID ( int num, size_t *indices, const double *in_vals )
{
    if ( num==0 )
        return;
    INCREMENT_COUNT("Virtual");
    std::vector<std::vector<size_t> >  ndxs;
    std::vector<std::vector<double> >  vals;
    partitionGlobalValues( num, indices, in_vals, ndxs, vals );
    for (size_t i = 0 ; i != ndxs.size() ; i++ ) {
        if ( ndxs[i].size() )
            d_vVectors[i]->setLocalValuesByGlobalID ( ndxs[i].size(), &(ndxs[i][0]), &(vals[i][0]) );
    }
}
void MultiVector::setGhostValuesByGlobalID ( int num, size_t *indices, const double *in_vals )
{
    if ( num==0 )
        return;
    INCREMENT_COUNT("Virtual");
    std::vector<std::vector<size_t> >  ndxs;
    std::vector<std::vector<double> >  vals;
    partitionGlobalValues( num, indices, in_vals, ndxs, vals );
    for (size_t i = 0 ; i != ndxs.size() ; i++ ) {
        if ( ndxs[i].size() )
            d_vVectors[i]->setGhostValuesByGlobalID ( ndxs[i].size(), &(ndxs[i][0]), &(vals[i][0]) );
    }
}
void MultiVector::setValuesByGlobalID ( int num, size_t *indices, const double *in_vals )
{
    if ( num==0 )
        return;
    INCREMENT_COUNT("Virtual");
    std::vector<std::vector<size_t> >  ndxs;
    std::vector<std::vector<double> >  vals;
    partitionGlobalValues( num, indices, in_vals, ndxs, vals );
    for (size_t i = 0 ; i != ndxs.size() ; i++ ) {
        if ( ndxs[i].size() )
            d_vVectors[i]->setValuesByGlobalID ( ndxs[i].size(), &(ndxs[i][0]), &(vals[i][0]) );
    }
}
void MultiVector::addValuesByLocalID ( int num, size_t *indices, const double *in_vals )
{
    if ( num==0 )
        return;
    INCREMENT_COUNT("Virtual");
    std::vector<std::vector<size_t> >  ndxs;
    std::vector<std::vector<double> >  vals;
    partitionLocalValues( num, indices, in_vals, ndxs, vals );
    for (size_t i = 0 ; i != ndxs.size() ; i++ ) {
        if ( ndxs[i].size() )
            d_vVectors[i]->addValuesByLocalID ( ndxs[i].size(), &(ndxs[i][0]), &(vals[i][0]) );
    }
}
void MultiVector::addLocalValuesByGlobalID ( int num, size_t *indices, const double *in_vals )
{
    if ( num==0 )
        return;
    INCREMENT_COUNT("Virtual");
    std::vector<std::vector<size_t> >  ndxs;
    std::vector<std::vector<double> >  vals;
    partitionGlobalValues( num, indices, in_vals, ndxs, vals );
    for (size_t i = 0 ; i != ndxs.size() ; i++ ) {
        if ( ndxs[i].size() )
            d_vVectors[i]->addLocalValuesByGlobalID ( ndxs[i].size(), &(ndxs[i][0]), &(vals[i][0]) );
    }
}
void MultiVector::addValuesByGlobalID ( int num, size_t *indices, const double *in_vals )
{
    if ( num==0 )
        return;
    INCREMENT_COUNT("Virtual");
    std::vector<std::vector<size_t> >  ndxs;
    std::vector<std::vector<double> >  vals;
    partitionGlobalValues( num, indices, in_vals, ndxs, vals );
    for (size_t i = 0 ; i != ndxs.size() ; i++ ) {
        if ( ndxs[i].size() )
            d_vVectors[i]->addValuesByGlobalID ( ndxs[i].size(), &(ndxs[i][0]), &(vals[i][0]) );
    }
}
void MultiVector::getValuesByGlobalID ( int num, size_t *indices, double *out_vals ) const
{
    if ( num==0 )
        return;
    INCREMENT_COUNT("Virtual");
    std::vector<std::vector<size_t> >  ndxs;
    std::vector<std::vector<double> >  vals;
    std::vector<std::vector<int> >  remap;
    partitionGlobalValues ( num, indices, out_vals, ndxs, vals, &remap );
    for (size_t i = 0 ; i != ndxs.size() ; i++ ) {
        if ( ndxs[i].size() )
            d_vVectors[i]->getValuesByGlobalID ( ndxs[i].size(), &(ndxs[i][0]), &(vals[i][0]) );
    }
    for (size_t i = 0 ; i != remap.size() ; i++ ) {
        for (size_t j = 0 ; j != remap[i].size() ; j++ )
            out_vals[remap[i][j]] = vals[i][j];
    }
}
void MultiVector::getLocalValuesByGlobalID ( int num, size_t *indices, double *out_vals ) const
{
    if ( num==0 )
        return;
    INCREMENT_COUNT("Virtual");
    std::vector<std::vector<size_t> >  ndxs;
    std::vector<std::vector<double> >  vals;
    std::vector<std::vector<int> >  remap;
    partitionGlobalValues ( num, indices, out_vals, ndxs, vals, &remap );
    for (size_t i = 0 ; i != ndxs.size() ; i++ ) {
        if ( ndxs[i].size() )
            d_vVectors[i]->getLocalValuesByGlobalID ( ndxs[i].size(), &(ndxs[i][0]), &(vals[i][0]) );
    }
    for (size_t i = 0 ; i != remap.size() ; i++ ) {
        for (size_t j = 0 ; j != remap[i].size() ; j++ )
            out_vals[remap[i][j]] = vals[i][j];
    }
}
void MultiVector::getGhostValuesByGlobalID ( int num, size_t *indices, double *out_vals ) const
{
    if ( num==0 )
        return;
    INCREMENT_COUNT("Virtual");
    std::vector<std::vector<size_t> >  ndxs;
    std::vector<std::vector<double> >  vals;
    std::vector<std::vector<int> >  remap;
    partitionGlobalValues ( num, indices, out_vals, ndxs, vals, &remap );
    for (size_t i = 0 ; i != ndxs.size() ; i++ ) {
        if ( ndxs[i].size() )
            d_vVectors[i]->getGhostValuesByGlobalID ( ndxs[i].size(), &(ndxs[i][0]), &(vals[i][0]) );
    }
    for (size_t i = 0 ; i != remap.size() ; i++ ) {
        for (size_t j = 0 ; j != remap[i].size() ; j++ )
            out_vals[remap[i][j]] = vals[i][j];
    }
}
void MultiVector::getValuesByLocalID ( int num, size_t *indices, double *out_vals ) const
{
    if ( num==0 )
        return;
    INCREMENT_COUNT("Virtual");
    std::vector<std::vector<size_t> >  ndxs;
    std::vector<std::vector<double> >  vals;
    std::vector<std::vector<int> >  remap;
    partitionLocalValues ( num, indices, out_vals, ndxs, vals, &remap );
    for (size_t i = 0 ; i != ndxs.size() ; i++ ) {
        if ( ndxs[i].size() )
            d_vVectors[i]->getValuesByLocalID ( ndxs[i].size(), &(ndxs[i][0]), &(vals[i][0]) );
    }
    for (size_t i = 0 ; i != remap.size() ; i++ ) {
        for (size_t j = 0 ; j != remap[i].size() ; j++ )
            out_vals[remap[i][j]] = vals[i][j];
    }
}


/****************************************************************
* Function to partition the global ids by the sub vectors       *
****************************************************************/
void  MultiVector::partitionGlobalValues ( const int num, const size_t *indices, const double *vals,
    std::vector<std::vector<size_t> > &out_indices, std::vector<std::vector<double> > &out_vals, std::vector<std::vector<int> > *remap ) const
{
    PROFILE_START("partitionGlobalValues",2);
    const size_t neg_one = ~((size_t)0);
    std::vector<size_t> globalDOFs(num,neg_one);
    for (int i=0; i<num; i++)
        globalDOFs[i] = indices[i];
    out_indices.resize( d_vVectors.size() );
    out_vals.resize( d_vVectors.size() );
    if ( remap != NULL ) 
        remap->resize( d_vVectors.size() );
    AMP::Discretization::multiDOFManager* manager = (AMP::Discretization::multiDOFManager*) d_DOFManager.get();
    std::vector<AMP::Discretization::DOFManager::shared_ptr> DOFManagers = manager->getDOFManagers();
    AMP_ASSERT(DOFManagers.size()==d_vVectors.size());
    for (size_t i=0; i<d_vVectors.size(); i++) {
        AMP_ASSERT(d_vVectors[i]->getDOFManager()==DOFManagers[i]);
        std::vector<size_t> subDOFs = manager->getSubDOF( i, globalDOFs );
        size_t count = 0; 
        for (size_t j=0; j<subDOFs.size(); j++) {
            if ( subDOFs[j] != neg_one )
                count++;
        }
        out_indices[i] = std::vector<size_t>(count,neg_one);
        out_vals[i] = std::vector<double>(count,0.0);
        if ( remap != NULL ) 
            remap->operator[](i) = std::vector<int>(count,neg_one);
        count = 0; 
        for (size_t j=0; j<subDOFs.size(); j++) {
            if ( subDOFs[j] != neg_one ) {
                out_indices[i][count] = subDOFs[j];
                out_vals[i][count] = vals[j];
                if ( remap != NULL )
                    remap->operator[](i)[count] = j;
                count++;
            }
        }
    }
    PROFILE_STOP("partitionGlobalValues",2);
}


/****************************************************************
* Function to partition the local ids by the sub vectors       *
****************************************************************/
void  MultiVector::partitionLocalValues ( const int num, const size_t *indices, const double *vals,
    std::vector<std::vector<size_t> > &out_indices, std::vector<std::vector<double> > &out_vals, std::vector<std::vector<int> > *remap ) const
{
    if ( num == 0 )
        return;
    PROFILE_START("partitionLocalValues",2);
    // Convert the local ids to global ids
    size_t begin_DOF = d_DOFManager->beginDOF();
    size_t end_DOF = d_DOFManager->endDOF();
    std::vector<size_t> global_indices(num);
    for (int i=0; i<num; i++) {
        AMP_INSIST(indices[i]<end_DOF,"Invalid local id");
        global_indices[i] = indices[i] + begin_DOF;
    }
    // Partition based on the global ids
    partitionGlobalValues( num, &global_indices[0], vals, out_indices, out_vals, remap );
    // Convert the new global ids back to local ids
    const size_t neg_one = ~((size_t)0);
    for (size_t i=0; i<d_vVectors.size(); i++) {
        if ( out_indices[i].size()==0 )
            continue;
        begin_DOF = d_vVectors[i]->getDOFManager()->beginDOF();
        end_DOF = d_vVectors[i]->getDOFManager()->endDOF();
        for (size_t j=0; j<out_indices[i].size(); j++) {
            AMP_ASSERT(out_indices[i][j]!=neg_one);
            out_indices[i][j] -= begin_DOF;
            AMP_ASSERT(out_indices[i][j]<end_DOF);
        }
    }
    PROFILE_STOP("partitionLocalValues",2);
}


}
}

