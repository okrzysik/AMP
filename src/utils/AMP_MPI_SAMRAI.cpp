#include "AMP/AMP_TPLs.h"
#include "AMP/utils/AMP_MPI.h"

#ifdef AMP_USE_SAMRAI
    #include "SAMRAI/tbox/SAMRAI_MPI.h"


    #ifdef AMP_USE_MPI

/****************************************************************************
 * Interface to SAMRAI (with MPI)                                            *
 ****************************************************************************/
AMP::AMP_MPI::AMP_MPI( const SAMRAI::tbox::SAMRAI_MPI &comm )
    : AMP_MPI( comm.getCommunicator(), false )
{
}
AMP::AMP_MPI::operator SAMRAI::tbox::SAMRAI_MPI() const
{
    return SAMRAI::tbox::SAMRAI_MPI( d_comm );
}

    #else

/****************************************************************************
 * Interface to SAMRAI (without MPI)                                         *
 ****************************************************************************/
AMP::AMP_MPI::AMP_MPI( const SAMRAI::tbox::SAMRAI_MPI &comm ) : AMP_MPI()
{
    if ( comm.hasNullCommunicator() ) {
        *this = AMP_MPI( AMP_COMM_NULL );
    } else if ( comm.getCommunicator() == MPI_COMM_WORLD ) {
        *this = AMP_MPI( AMP_COMM_WORLD );
    } else if ( comm == comm.getSAMRAIWorld() ) {
        *this = AMP_MPI( AMP_COMM_WORLD );
    } else if ( comm.getCommunicator() == MPI_COMM_SELF ) {
        *this = AMP_MPI( AMP_COMM_SELF );
    } else {
        *this = AMP_MPI( comm.getCommunicator() );
    }
}
AMP::AMP_MPI::operator SAMRAI::tbox::SAMRAI_MPI() const
{
    if ( isNull() )
        return SAMRAI::tbox::SAMRAI_MPI( MPI_COMM_NULL );
    return SAMRAI::tbox::SAMRAI_MPI( MPI_COMM_WORLD );
}

    #endif

#endif
