#include "AMP/utils/AMP_MPI.I"
#include "AMP/utils/Database.hpp"


/****************************************************************************
 * pack/unpack functions                                                     *
 ****************************************************************************/
template<>
size_t packSize<AMP::Database>( AMP::Database const & )
{
    AMP_ERROR( "Not finished" );
    return 0;
}
template<>
size_t pack<AMP::Database>( const AMP::Database &, std::byte * )
{
    AMP_ERROR( "Not finished" );
    return 0;
}
template<>
size_t unpack<AMP::Database>( AMP::Database &, const std::byte * )
{
    AMP_ERROR( "Not finished" );
    return 0;
}


/****************************************************************************
 * Explicit instantiation                                                    *
 ****************************************************************************/
INSTANTIATE_MPI_BCAST( std::shared_ptr<AMP::Database> );
INSTANTIATE_MPI_SENDRECV( std::shared_ptr<AMP::Database> );
