#include "Vector.h"


#include <typeinfo>
#include <float.h>
#include <math.h>

#include "utils/AMP_MPI.h"
#include "utils/Utilities.h"

#include "DataChangeFirer.h"

#include "MultiVector.h"
#include "VectorSelector.h"

namespace AMP {
namespace LinearAlgebra {


RNG::shared_ptr  Vector::d_DefaultRNG;


void Vector::selectInto ( const VectorSelector &s , Vector::shared_ptr retVal )
{
    if ( s.isSelected ( shared_from_this () ) )
    {
      retVal->castTo<MultiVector> ().addVector ( s.getSubsetter()->subset ( shared_from_this () ) );
    }
}


Vector::shared_ptr  Vector::select ( const VectorSelector &s , const std::string &name )
{
    Vector::shared_ptr  retVal = MultiVector::create ( name , getComm() );
    selectInto ( s , retVal );
    if ( retVal->numberOfDataBlocks() )
      return retVal;
    return Vector::shared_ptr ();
}


void Vector::registerView ( Vector::shared_ptr v )
{
    for ( size_t i = 0 ; i != d_Views.size() ; i++ )
      if ( d_Views[i].lock() == v )
        return;

    d_Views.push_back ( v );
}


void Vector::setCommunicationList ( CommunicationList::shared_ptr  comm )
{
    AMP_ASSERT ( comm );
    d_CommList = comm;
    if ( comm )
    {
      addCommunicationListToParameters ( comm );
      d_Ghosts = boost::shared_ptr<std::vector<double> > ( new std::vector<double> ( d_CommList->getVectorReceiveBufferSize() ) );
      d_AddBuffer = boost::shared_ptr<std::vector<double> > ( new std::vector<double> ( d_CommList->getVectorReceiveBufferSize() ) );
    }
}


Vector::shared_ptr  Vector::subsetVectorForVariable ( const Variable::shared_ptr  &name )
{
    Vector::shared_ptr retVal;
    if ( d_pVariable )  // If there is a variable...
    {
      if ( *d_pVariable == *name )
      {
        retVal = shared_from_this ();
      }
    }
    return retVal;
}


void Vector::copyVector ( const Vector &other )
{
    if ( other.getLocalSize() != getLocalSize() )
    {  // Another error condition
      AMP_ERROR( "Destination vector and source vector not the same size" );
    }
    VectorDataIterator  cur1 = begin();
    VectorDataIterator  end1 = end();
    ConstVectorDataIterator  cur2 = other.begin();
    while ( cur1 != end1 )
    {
      *cur1 = *cur2;
      cur1++; cur2++;
    }
    if ( isA<DataChangeFirer>() )
    {
      castTo<DataChangeFirer>().fireDataChange();
    }
    copyGhostValues ( other );
}


Vector::Vector ()
{
    d_Ghosts = boost::shared_ptr<std::vector<double> > ( new std::vector<double> );
    d_AddBuffer = boost::shared_ptr<std::vector<double> > ( new std::vector<double> );
    d_UpdateState.reset( new UpdateState );
    *d_UpdateState = NOT_UPDATING;
}


Vector::shared_ptr Vector::cloneVector ( const std::string &name ) const
{
    Vector::shared_ptr retVal;
    if ( getVariable() )
    {
      retVal = cloneVector ( getVariable()->cloneVariable ( name ) );
    }
    else
    {
      retVal = cloneVector ();
    }
    return retVal;
}


bool Vector::equals ( Vector &rhs , double tol )
{
    int RetVal = 0;
    if (( getGlobalSize() == rhs.getGlobalSize() ) && ( getLocalSize() == rhs.getLocalSize() ))
    {
      VectorDataIterator cur1 = begin();
      VectorDataIterator cur2 = rhs.begin();
      VectorDataIterator last = end();
      bool failed = false;
      while ( cur1 != last )
      {
        double v1 = *cur1;
        double v2 = *cur2;
        if ( fabs ( v1 - v2 ) > tol )
        {
          failed = true;
          break;
        }
        ++cur1; ++cur2;
      }
      if ( !failed )
      {
        RetVal = 1;
      }
    }

    int ans;
    if ( d_CommList ) {
        ans = d_CommList->getComm().minReduce(RetVal);
    } else {
        ans = RetVal;
    }

    return ans == 1 ? true : false;
}


// The following two functions are Jungho's
// This will fail when y_i = 0... Needs to be fixed
double Vector::minQuotient(const VectorOperations &x, const VectorOperations &y)
{
    const Vector &x_vec = x.castTo<Vector>();
    const Vector &y_vec = y.castTo<Vector>();
    Vector::const_iterator curx = x_vec.begin();
    Vector::const_iterator cury = y_vec.begin();
    Vector::const_iterator endx = x_vec.end();
    while ( curx != endx )
    {
      if ( *cury != 0.0 ) break;
      curx++;
      cury++;
    }
    //Probably should do a test across processors, not just the one
    AMP_INSIST ( curx != endx , "denominator is the zero vector on an entire process" );
    double myRetVal = (*curx) / (*cury);
    while ( curx != endx )
    {
      if ( *cury != 0.0 )
      {
        myRetVal = std::min ( myRetVal , (*curx)/(*cury) );
      }
      curx++;
      cury++;
    }
    double retVal = x_vec.getComm().minReduce(myRetVal);
    return retVal;
}


double Vector::wrmsNorm(const VectorOperations &x, const VectorOperations &y)
{
    double dot_prod=0.0;
    int global_size=0;
    Vector::shared_ptr  temp_vec = x.castTo<Vector>().cloneVector ();
    temp_vec->multiply(x, y);

    dot_prod = temp_vec->dot(temp_vec);
    global_size = temp_vec->getGlobalSize();

    return(sqrt(dot_prod/global_size));
}


double Vector::wrmsNormMask ( const VectorOperations &x , const VectorOperations &y , const VectorOperations &mask )
{
    double dot_prod = 0.0;
    const Vector &x_vec = x.castTo<Vector> ();
    const Vector &y_vec = y.castTo<Vector> ();
    const Vector &m_vec = mask.castTo<Vector> ();

    Vector::const_iterator curx = x_vec.begin();
    Vector::const_iterator endx = x_vec.end();
    Vector::const_iterator cury = y_vec.begin();
    Vector::const_iterator curm = m_vec.begin();
    while ( curx != endx )
    {
      if ( *curm > 0.0 )
      {
        dot_prod += (*curx)*(*curx)*(*cury)*(*cury);
      }
      curx++; cury++; curm++;
    }
    AMP_ASSERT ( cury == y_vec.end() );
    AMP_ASSERT ( curm == m_vec.end() );
    double all_dot = x_vec.getComm().sumReduce(dot_prod);
    return sqrt ( all_dot / (double) x_vec.getGlobalSize() );
}


#define DESCRIPTOR_ID_ARRAY_SCRATCH_SPACE (10)


/************************************************************************
*                                                                       *
* The constructor for Vector<DIM> objects initializes vector structure. *
*                                                                       *
************************************************************************/
Vector::Vector( VectorParameters::shared_ptr  parameters)
{

    // Set default output stream
    d_output_stream = &AMP::plog;

    setCommunicationList ( parameters->d_CommList );
    d_UpdateState.reset( new UpdateState );
    *d_UpdateState = NOT_UPDATING;
    d_DOFManager = parameters->d_DOFManager;
}


void Vector::makeConsistent ( ScatterType  t )
{
    if ( t == CONSISTENT_ADD )
    {
      AMP_ASSERT ( *d_UpdateState != SETTING );
      std::vector<double>  send_vec_add ( d_CommList->getVectorReceiveBufferSize() );
      std::vector<double>  recv_vec_add ( d_CommList->getVectorSendBufferSize() );
      d_CommList->packReceiveBuffer ( send_vec_add , *this );
      d_CommList->scatter_add ( send_vec_add , recv_vec_add );
      d_CommList->unpackSendBufferAdd ( recv_vec_add , *this );
      for ( std::vector<double>::iterator curAdd = d_AddBuffer->begin() ;
            curAdd != d_AddBuffer->end() ;
            curAdd++ )
      {
        *curAdd = 0.0;
      }
    }
    *d_UpdateState = SETTING;
    std::vector<double>  send_vec ( d_CommList->getVectorSendBufferSize() );
    std::vector<double>  recv_vec ( d_CommList->getVectorReceiveBufferSize() );
    d_CommList->packSendBuffer ( send_vec , *this );
    d_CommList->scatter_set ( send_vec , recv_vec );
    d_CommList->unpackReceiveBufferSet ( recv_vec , *this );
    *d_UpdateState = NOT_UPDATING;
}


void Vector::setValuesByGlobalID ( int numVals , int *ndx , const double *vals )
{
    INCREMENT_COUNT("Virtual");
    AMP_ASSERT ( *d_UpdateState != ADDING );
    *d_UpdateState = SETTING;
    for ( unsigned int i = 0 ; i != (unsigned int) numVals ; i++ )
    {
      if ( ( ((unsigned int) ndx[i]) < getLocalStartID() ) ||
           ( ((unsigned int) ndx[i]) >= (getLocalStartID() + getLocalMaxID()) ) )
      {
        (*d_Ghosts)[d_CommList->getLocalGhostID ( ndx[i] )] = vals[i];
      }
      else
      {
        setLocalValuesByGlobalID ( 1 , ndx+i , vals+i );
      }
    }
}


void Vector::addValuesByGlobalID ( int numVals , int *ndx , const double *vals )
{
    INCREMENT_COUNT("Virtual");
    AMP_ASSERT ( *d_UpdateState != SETTING );
    *d_UpdateState = ADDING;
    for ( unsigned int i = 0 ; i != (unsigned int) numVals ; i++ )
    {
      if ( ( ((unsigned int) ndx[i]) < getLocalStartID() ) ||
           ( ((unsigned int) ndx[i]) >= (getLocalStartID() + getLocalMaxID()) ) )
      {
        (*d_AddBuffer)[d_CommList->getLocalGhostID ( ndx[i] )] += vals[i];
      }
      else
      {
        addLocalValuesByGlobalID ( 1 , ndx+i , vals+i );
      }
    }
}


void Vector::getValuesByLocalID ( int num , int *ndx , double *vals ) const
{
    INCREMENT_COUNT("Virtual");
    for ( int i = 0 ; i != num ; i++ )
    {
      int block_number = 0;
      int offset = ndx[i];
      while ( offset >= (int)sizeOfDataBlock ( block_number ) )
      {
        offset -= (int)sizeOfDataBlock ( block_number );
        block_number++;
        if ( block_number >= (int)numberOfDataBlocks() )
        {
          AMP_ERROR( "Bad local id!" );
        }
      }
      vals[i] = getRawDataBlock<double> ( block_number ) [ offset ];
    }
}


void Vector::getValuesByGlobalID ( int numVals , int *ndx , double *vals ) const
{
    INCREMENT_COUNT("Virtual");
    for ( unsigned int i = 0 ; i != (unsigned int) numVals ; i++ )
    {
      if ( ( ((unsigned int) ndx[i]) < getLocalStartID() ) ||
           ( ((unsigned int) ndx[i]) >= (getLocalStartID() + getLocalMaxID()) ) )
      {
        vals[i] = (*d_Ghosts)[d_CommList->getLocalGhostID ( ndx[i] )] +
                  (*d_AddBuffer)[d_CommList->getLocalGhostID ( ndx[i] )];
      }
      else
      {
        getLocalValuesByGlobalID ( 1 , ndx+i , vals+i );
      }
    }
}


void Vector::getGhostAddValuesByGlobalID ( int numVals , int *ndx , double *vals ) const
{
    INCREMENT_COUNT("Virtual");
    for ( unsigned int i = 0 ; i != (unsigned int) numVals ; i++ )
    {
      if ( ( ((unsigned int) ndx[i]) < getLocalStartID() ) ||
           ( ((unsigned int) ndx[i]) >= (getLocalStartID() + getLocalMaxID()) ) )
      {
        vals[i] = (*d_AddBuffer)[d_CommList->getLocalGhostID ( ndx[i] )];
      }
      else
      {
        AMP_ERROR( "Tried to get a non-ghost ghost value" );
      }
    }
}


void  Vector::dumpOwnedData ( std::ostream &out , size_t GIDoffset , size_t LIDoffset ) const
{
    const_iterator  curElement = begin();
    size_t gid = GIDoffset;
    if ( getCommunicationList() )
      gid += getCommunicationList()->getStartGID();
    size_t lid = LIDoffset;
    while ( curElement != end() )
    {
      out << "  GID: " << gid
          << "  LID: " << lid
          << "  Value: " << *curElement << "\n";
      curElement++;
      gid++;
      lid++;
    }
}

  
void  Vector::dumpGhostedData ( std::ostream &out , size_t offset ) const
{
    if ( !getCommunicationList() ) return;
    const std::vector<unsigned int> &ghosts = getCommunicationList()->getGhostIDList();
    std::vector<unsigned int>::const_iterator curGhost = ghosts.begin();
    std::vector<double>::iterator curVal = d_Ghosts->begin();
    while ( curVal != d_Ghosts->end() )
    {
      out << "  GID: " << (*curGhost + offset ) << "  Value: " << (*curVal) << "\n";
      curGhost++;
      curVal++;
    }
}

  
std::ostream &operator << ( std::ostream &out , const Vector &v )
{
    out << "Vector type: " << v.type() << "\n";
    if ( v.getVariable() )
    {
      out << "Variable name: " << v.getVariable()->getName() << "\n";
    }
    if ( v.getCommunicationList() )
    {
      int rank = v.getComm().getRank();
      out << "Processor: " << rank << "\n";
    }
    out << "\n"
        << "Number of owned elements: " << v.getLocalSize() << "\n"
        << "Number of ghosted elements: " << v.getGhostSize() << "\n";
    if ( v.numberOfDataBlocks() > 1 )
      out << "Number of sub-vectors: " << v.numberOfDataBlocks() << "\n";
    out << "\n"
        << "Data block pointers: \n";

    for ( size_t i = 0 ; i != v.numberOfDataBlocks() ; i++ )
      out << "                     " << v.getRawDataBlock<double> ( i ) << "  ( " << v.sizeOfDataBlock(i) << " elements )\n";

    out << "\n"
        << "Local Data:\n";
    v.dumpOwnedData ( out );

    out << "\n"
        << "Ghosted Data:\n";
    v.dumpGhostedData ( out );

    return out;
}

  
void Vector::scale(double alpha, const VectorOperations &x)
{
    const Vector &t_c = x.castTo<Vector>();
    AMP_ASSERT ( getLocalSize() == t_c.getLocalSize() );
    iterator  curMe = begin();
    iterator  last = end();
    const_iterator  curRhs = t_c.begin();
    while ( curMe != last )
    {
      *curMe = alpha * *curRhs;
      curRhs++;
      curMe++;
    }
    dataChanged ();
}


void Vector::add(const VectorOperations &x, const VectorOperations &y)
{
    const Vector &t_x = x.castTo<const Vector>();
    const Vector &t_y = y.castTo<const Vector>();
    AMP_ASSERT ( getLocalSize() == t_x.getLocalSize() );
    AMP_ASSERT ( getLocalSize() == t_y.getLocalSize() );
    iterator  curMe = begin();
    iterator  last = end();
    const_iterator  curXRhs = t_x.begin();
    const_iterator  curYRhs = t_y.begin();
    while ( curMe != last )
    {
      *curMe = *curXRhs + *curYRhs;
      curXRhs++;
      curYRhs++;
      curMe++;
    }
    dataChanged ();
}


void Vector::subtract(const VectorOperations &x, const VectorOperations &y)
{
    const Vector &t_x = x.castTo<const Vector>();
    const Vector &t_y = y.castTo<const Vector>();
    AMP_ASSERT ( getLocalSize() == t_x.getLocalSize() );
    AMP_ASSERT ( getLocalSize() == t_y.getLocalSize() );
    iterator  curMe = begin();
    iterator  last = end();
    const_iterator  curXRhs = t_x.begin();
    const_iterator  curYRhs = t_y.begin();
    while ( curMe != last )
    {
      *curMe = *curXRhs - *curYRhs;
      curXRhs++;
      curYRhs++;
      curMe++;
    }
    dataChanged ();
}

void Vector::multiply( const VectorOperations &x, const VectorOperations &y)
{
    const Vector &t_x = x.castTo<const Vector>();
    const Vector &t_y = y.castTo<const Vector>();
    AMP_ASSERT ( getLocalSize() == t_x.getLocalSize() );
    AMP_ASSERT ( getLocalSize() == t_y.getLocalSize() );
    iterator  curMe = begin();
    iterator  last = end();
    const_iterator  curXRhs = t_x.begin();
    const_iterator  curYRhs = t_y.begin();
    while ( curMe != last )
    {
      *curMe = *curXRhs * *curYRhs;
      curXRhs++;
      curYRhs++;
      curMe++;
    }
    dataChanged ();
}


void Vector::divide( const VectorOperations &x, const VectorOperations &y)
{
    const Vector &t_x = x.castTo<const Vector>();
    const Vector &t_y = y.castTo<const Vector>();
    AMP_ASSERT ( getLocalSize() == t_x.getLocalSize() );
    AMP_ASSERT ( getLocalSize() == t_y.getLocalSize() );
    iterator  curMe = begin();
    iterator  last = end();
    const_iterator  curXRhs = t_x.begin();
    const_iterator  curYRhs = t_y.begin();
    while ( curMe != last )
    {
      *curMe = *curXRhs / *curYRhs;
      curXRhs++;
      curYRhs++;
      curMe++;
    }
    dataChanged ();
}


void Vector::reciprocal(const VectorOperations &x)
{
    const Vector &t_c = x.castTo<Vector>();
    AMP_ASSERT ( getLocalSize() == t_c.getLocalSize() );
    iterator  curMe = begin();
    iterator  last = end();
    const_iterator  curRhs = t_c.begin();
    while ( curMe != last )
    {
      *curMe = 1. / *curRhs;
      curRhs++;
      curMe++;
    }
    dataChanged ();
}


void Vector::linearSum(double alpha, const VectorOperations &x, double beta, const VectorOperations &y)
{
    const Vector &t_x = x.castTo<const Vector>();
    const Vector &t_y = y.castTo<const Vector>();
    AMP_ASSERT ( getLocalSize() == t_x.getLocalSize() );
    AMP_ASSERT ( getLocalSize() == t_y.getLocalSize() );
    iterator  curMe = begin();
    iterator  last = end();
    const_iterator  curXRhs = t_x.begin();
    const_iterator  curYRhs = t_y.begin();
    while ( curMe != last )
    {
      *curMe = alpha * *curXRhs + beta * *curYRhs;
      curXRhs++;
      curYRhs++;
      curMe++;
    }
    dataChanged ();
}


void Vector::axpy(double alpha, const VectorOperations &x, const VectorOperations &y)
{
    const Vector &t_x = x.castTo<const Vector>();
    const Vector &t_y = y.castTo<const Vector>();
    AMP_ASSERT ( getLocalSize() == t_x.getLocalSize() );
    AMP_ASSERT ( getLocalSize() == t_y.getLocalSize() );
    iterator  curMe = begin();
    iterator  last = end();
    const_iterator  curXRhs = t_x.begin();
    const_iterator  curYRhs = t_y.begin();
    while ( curMe != last )
    {
      *curMe = alpha * *curXRhs + *curYRhs;
      curXRhs++;
      curYRhs++;
      curMe++;
    }
    dataChanged ();
}


void Vector::axpby(double alpha, double beta, const VectorOperations &x)
{
    const Vector &t_x = x.castTo<const Vector>();
    AMP_ASSERT ( getLocalSize() == t_x.getLocalSize() );
    iterator  curMe = begin();
    iterator  last = end();
    const_iterator  curXRhs = t_x.begin();
    while ( curMe != last )
    {
      *curMe = alpha * *curXRhs + beta * *curMe;
      curXRhs++;
      curMe++;
    }
    dataChanged ();
}


void Vector::abs(const VectorOperations &x)
{
    const Vector &t_x = x.castTo<const Vector>();
    AMP_ASSERT ( getLocalSize() == t_x.getLocalSize() );
    iterator  curMe = begin();
    iterator  last = end();
    const_iterator  curXRhs = t_x.begin();
    while ( curMe != last )
    {
      *curMe = fabs ( *curXRhs );
      curXRhs++;
      curMe++;
    }
    dataChanged ();
}


void  Vector::setRandomValues ()
{
    RandomVariable<double>  r ( 0. , 1. , getDefaultRNG() );
    iterator  curMe = begin();
    iterator  last = end();
    while ( curMe != last )
    {
      double curRand = r;
      *curMe = curRand;
      curMe++;
    }
    dataChanged ();
}


void  Vector::setRandomValues ( RNG::shared_ptr rng )
{
    RandomVariable<double>  r ( 0. , 1. , rng );
    iterator  curMe = begin();
    iterator  last = end();
    while ( curMe != last )
    {
      *curMe = r;
      curMe++;
    }
    dataChanged ();
}


double Vector::dot(const VectorOperations &x) const
{
    const Vector &t_x = x.castTo<const Vector>();
    AMP_ASSERT ( getLocalSize() == t_x.getLocalSize() );
    const_iterator  curMe = begin();
    const_iterator  last = end();
    const_iterator  curXRhs = t_x.begin();
    double ans = 0.;
    while ( curMe != last )
    {
      ans += *curMe * *curXRhs;
      curXRhs++;
      curMe++;
    }
    double retVal;
    if ( getCommunicationList() ) {
        retVal = getComm().sumReduce(ans);
    } else {
        retVal = ans;
    }
    return retVal;
}


double Vector::min () const
{
    const_iterator  curMe = begin();
    const_iterator  last = end();
    double ans = *curMe;
    curMe++;
    while ( curMe != last )
    {
      ans = std::min ( *curMe , ans );
      curMe++;
    }
    double retVal;
    if ( getCommunicationList() ) {
        retVal = getComm().minReduce(ans);
    } else {
        retVal = ans;
    }
    return retVal;
}


double Vector::max () const
{
    const_iterator  curMe = begin();
    const_iterator  last = end();
    double ans = *curMe;
    curMe++;
    while ( curMe != last )
    {
      ans = std::max ( *curMe , ans );
      curMe++;
    }
    double retVal;
    if ( getCommunicationList() ) {
        retVal = getComm().maxReduce(ans);
    } else {
        retVal = ans;
    }
    return retVal;
}


double Vector::L1Norm(void) const
{
    const_iterator  curMe = begin();
    const_iterator  last = end();
    double ans = 0.;
    while ( curMe != last )
    {
      ans += fabs ( *curMe );
      curMe++;
    }
    double retVal;
    if ( getCommunicationList() ) {
        retVal = getComm().sumReduce(ans);
    } else {
        retVal = ans;
    }
    return retVal;
}


double Vector::maxNorm(void) const
{
    const_iterator  curMe = begin();
    const_iterator  last = end();
    double ans = 0.;
    while ( curMe != last )
    {
      ans = std::max ( ans , fabs ( *curMe ) );
      curMe++;
    }
    double retVal;
    if ( getCommunicationList() ) {
        retVal = getComm().maxReduce(ans);
    } else {
        retVal = ans;
    }
    return retVal;
}


double Vector::L2Norm(void) const
{
    const_iterator  curMe = begin();
    const_iterator  last = end();
    double ans = 0.;
    while ( curMe != last )
    {
      ans += *curMe * *curMe;
      curMe++;
    }
    double retVal;
    if ( getCommunicationList() ) {
        retVal = getComm().sumReduce(ans);
    } else {
        retVal = ans;
    }
    return sqrt ( retVal );
}


void Vector::scale(double alpha)
{
    iterator  curMe = begin();
    iterator  last = end();
    while ( curMe != last )
    {
      *curMe *= alpha;
      curMe++;
    }
    dataChanged();
}


void Vector::setToScalar(double alpha)
{
    iterator  curMe = begin();
    iterator  last = end();
    while ( curMe != last )
    {
      *curMe = alpha;
      curMe++;
    }
    dataChanged();
}


}
}

