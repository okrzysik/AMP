#include "CommCollectVector.h"
#include "ManagedVector.h"


// This needs a smarter communication pattern...
namespace AMP {
namespace LinearAlgebra {


void  CommCollectVector::makeRankComm ()
{
    int mySmallRank = d_SmallCommVector->getComm().getRank();
    int myBigRank = d_LargeComm.getRank();
    d_RankComm = d_LargeComm.split(mySmallRank,myBigRank);
}


const VectorEngine::shared_ptr  CommCollectVector::getVectorEngine () const
{
    VectorEngine::shared_ptr retval;
    if ( d_SmallCommVector->isA<VectorEngine> () )
        retval = boost::dynamic_pointer_cast<VectorEngine> ( d_SmallCommVector );
    else if ( d_SmallCommVector->isA<ManagedVector> () )
        retval = d_SmallCommVector->castTo<ManagedVector>().getVectorEngine();
    AMP_ASSERT ( retval );
    return retval;
}


VectorEngine::shared_ptr  CommCollectVector::getVectorEngine ()
{
    VectorEngine::shared_ptr retval;
    if ( d_SmallCommVector->isA<VectorEngine> () )
        retval = boost::dynamic_pointer_cast<VectorEngine> ( d_SmallCommVector );
    else if ( d_SmallCommVector->isA<ManagedVector> () )
        retval = d_SmallCommVector->castTo<ManagedVector>().getVectorEngine();
    AMP_ASSERT ( retval );
    return retval;
}


double CommCollectVector::min(void) const
{
    double t = d_SmallCommVector->min ();
    double retval = d_RankComm.minReduce(t);
    retval = d_SmallCommVector->getComm().bcast(retval,0);
    return retval;
}


double CommCollectVector::max(void) const
{
    double t = d_SmallCommVector->max ();
    double retval = d_RankComm.maxReduce(t);
    retval = d_SmallCommVector->getComm().bcast(retval,0);
    return retval;
}


double CommCollectVector::L1Norm(void) const
{
    double t = d_SmallCommVector->L1Norm ();
    double retval = d_RankComm.sumReduce(t);
    retval = d_SmallCommVector->getComm().bcast(retval,0);
    return retval;
}


double CommCollectVector::L2Norm(void) const
{
    double t = d_SmallCommVector->L2Norm ();
    t *= t;
    double retval = d_RankComm.sumReduce(t);
    retval = sqrt( retval );
    retval = d_SmallCommVector->getComm().bcast(retval,0);
    return retval;
}


double CommCollectVector::maxNorm(void) const
{
    double t = d_SmallCommVector->maxNorm ();
    double retval = d_RankComm.maxReduce(t);
    retval = d_SmallCommVector->getComm().bcast(retval,0);
    return retval;
}


size_t CommCollectVector::getGlobalSize() const
{
    int t = d_SmallCommVector->getGlobalSize();
    int retval = d_RankComm.sumReduce(t);
    retval = d_SmallCommVector->getComm().bcast(retval,0);
    return retval;
}


double CommCollectVector::dot(const VectorOperations &x) const
{
    if ( x.isA<CommCollectVector> () ) {
        double t = d_SmallCommVector->dot ( x.castTo<CommCollectVector>().d_SmallCommVector );
        double retval = d_RankComm.sumReduce(t);
        retval = d_SmallCommVector->getComm().bcast(retval,0);
        return retval;
    }
    return d_SmallCommVector->dot ( x );
}


} // LinearAlgebra
} // AMP

