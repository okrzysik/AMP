// This file contains the default instantiations for templated operations
// Note: Intel compilers need definitions before all default instantions to compile correctly
#ifndef included_AMP_MPI_I
#define included_AMP_MPI_I

#include "AMP/AMP_TPLs.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/TypeTraits.h"
#include "AMP/utils/UtilityMacros.h"

#include "ProfilerApp.h"

#include <complex>
#include <typeinfo>

#ifdef AMP_USE_MPI
    #ifndef USE_MPI
        #define USE_MPI
    #endif
#else
    #undef USE_MPI
#endif


#define MPI_CLASS AMP_MPI
#define MPI_CLASS_ERROR AMP_ERROR
#define MPI_CLASS_ASSERT AMP_ASSERT
#define MPI_CLASS_INSIST AMP_INSIST
#define MPI_CLASS_WARNING AMP_WARNING


// Include mpi.h (or define MPI objects)
// clang-format off
#ifdef USE_MPI
    // Building with mpi
    #include "mpi.h"
#elif defined( USE_SAMRAI ) && defined( USE_PETSC )
    // Building with SAMRAI and PETSc and without MPI is complicated
    #include "petsc/mpiuni/mpi.h"
    #undef MPI_REQUEST_NULL
    #define HAVE_MPI
    #ifndef included_tbox_SAMRAI_MPI
        extern int MPI_REQUEST_NULL;
        extern int MPI_ERR_IN_STATUS;
    #endif
    #ifdef INCLUDED_SAMRAI_CONFIG_H
        #include "SAMRAI/tbox/SAMRAI_MPI.h"
    #else
        #define INCLUDED_SAMRAI_CONFIG_H
        #define included_tbox_Utilities
        #include "SAMRAI/tbox/SAMRAI_MPI.h"
        #undef included_tbox_Utilities
        #undef INCLUDED_SAMRAI_CONFIG_H
    #endif
    #undef HAVE_MPI
#elif defined( USE_SAMRAI )
    // SAMRAI serial builds define basic MPI types
    #include "SAMRAI/tbox/SAMRAI_MPI.h"
#elif defined( USE_PETSC )
    // petsc serial builds include mpi.h
    #include "petsc/mpiuni/mpi.h"
#elif defined( USE_TRILINOS )
    // trilinos serial builds include mpi.h
    #include "mpi.h" 
#elif defined(__has_include)
    // Check if another package defines mpi.h
    #if __has_include("mpi.h")
        #include "mpi.h"
    #else
        #define SET_MPI_TYPES
    #endif
#else
    #define SET_MPI_TYPES
#endif
#ifdef SET_MPI_TYPES
    typedef uint32_t MPI_Comm;
    typedef uint32_t MPI_Request;
    typedef uint32_t MPI_Status;
    typedef void *MPI_Errhandler;
    #define MPI_COMM_NULL ( (MPI_Comm) 0xF4000010 )
    #define MPI_COMM_WORLD ( (MPI_Comm) 0xF4000012 )
    #define MPI_COMM_SELF ( (MPI_Comm) 0xF4000011 )
#endif
// clang-format on


#define MPI_CLASS_ERROR_TYPE( MSG ) \
    MPI_CLASS_ERROR( std::string( MSG ) + " (" + typeid( TYPE ).name() + "(" );


namespace AMP {


// Function to test if a type can be passed by MPI
template<typename TYPE>
struct is_mpi_copyable {
    template<typename T>
    static constexpr bool test()
    {
        if constexpr ( std::is_trivially_copyable<T>::value )
            return true;
        if constexpr ( is_pair<T>::value )
            return test<typename T::first_type>() && test<typename T::second_type>();
        return false;
    }
    static constexpr bool value = test<TYPE>();
};


/************************************************************************
 *  getType                                                              *
 ************************************************************************/
#ifdef USE_MPI
template<class TYPE>
struct TypeIntStruct {
    TYPE x;
    int i;
};
// clang-format off
template<typename T> struct TypeIntSelector;
template<> struct TypeIntSelector<char>     { using type = TypeIntStruct<short>; };
template<> struct TypeIntSelector<uint8_t>  { using type = TypeIntStruct<short>; };
template<> struct TypeIntSelector<int8_t>   { using type = TypeIntStruct<short>; };
template<> struct TypeIntSelector<uint16_t> { using type = TypeIntStruct<int>;   };
template<> struct TypeIntSelector<int16_t>  { using type = TypeIntStruct<short>; };
template<> struct TypeIntSelector<uint32_t> { using type = TypeIntStruct<long>;  };
template<> struct TypeIntSelector<int32_t>  { using type = TypeIntStruct<int>;   };
template<> struct TypeIntSelector<uint64_t> { using type = TypeIntStruct<long>;  };
template<> struct TypeIntSelector<int64_t>  { using type = TypeIntStruct<long>;  };
template<> struct TypeIntSelector<float>    { using type = TypeIntStruct<float>; };
template<> struct TypeIntSelector<double>   { using type = TypeIntStruct<double>; };
template<typename T> struct TypeIntSelector { using type = TypeIntStruct<T>;     };
// clang-format on
template<class TYPE>
constexpr bool isMPItype()
{
    return std::is_same<TYPE, char>::value || std::is_same<TYPE, uint8_t>::value ||
           std::is_same<TYPE, int8_t>::value || std::is_same<TYPE, uint16_t>::value ||
           std::is_same<TYPE, int16_t>::value || std::is_same<TYPE, uint32_t>::value ||
           std::is_same<TYPE, int32_t>::value || std::is_same<TYPE, uint64_t>::value ||
           std::is_same<TYPE, int64_t>::value || std::is_same<TYPE, float>::value ||
           std::is_same<TYPE, double>::value || std::is_same<TYPE, std::complex<float>>::value ||
           std::is_same<TYPE, std::complex<double>>::value ||
           std::is_same<TYPE, TypeIntStruct<short>>::value ||
           std::is_same<TYPE, TypeIntStruct<int>>::value ||
           std::is_same<TYPE, TypeIntStruct<long>>::value ||
           std::is_same<TYPE, TypeIntStruct<float>>::value ||
           std::is_same<TYPE, TypeIntStruct<double>>::value;
}
template<class TYPE>
MPI_Datatype MPI_getType()
{
    if constexpr ( std::is_same<TYPE, char>::value )
        return MPI_CHAR;
    if constexpr ( std::is_same<TYPE, uint8_t>::value )
        return MPI_UINT8_T;
    if constexpr ( std::is_same<TYPE, int8_t>::value )
        return MPI_INT8_T;
    if constexpr ( std::is_same<TYPE, uint16_t>::value )
        return MPI_UINT16_T;
    if constexpr ( std::is_same<TYPE, int16_t>::value )
        return MPI_INT16_T;
    if constexpr ( std::is_same<TYPE, uint32_t>::value )
        return MPI_UINT32_T;
    if constexpr ( std::is_same<TYPE, int32_t>::value )
        return MPI_INT32_T;
    if constexpr ( std::is_same<TYPE, uint64_t>::value )
        return MPI_UINT64_T;
    if constexpr ( std::is_same<TYPE, int64_t>::value )
        return MPI_INT64_T;
    if constexpr ( std::is_same<TYPE, float>::value )
        return MPI_FLOAT;
    if constexpr ( std::is_same<TYPE, double>::value )
        return MPI_DOUBLE;
    if constexpr ( std::is_same<TYPE, std::complex<float>>::value )
        return MPI_C_FLOAT_COMPLEX;
    if constexpr ( std::is_same<TYPE, std::complex<double>>::value )
        return MPI_C_DOUBLE_COMPLEX;
    if constexpr ( std::is_same<TYPE, TypeIntStruct<short>>::value )
        return MPI_SHORT_INT;
    if constexpr ( std::is_same<TYPE, TypeIntStruct<int>>::value )
        return MPI_2INT;
    if constexpr ( std::is_same<TYPE, TypeIntStruct<long>>::value )
        return MPI_LONG_INT;
    if constexpr ( std::is_same<TYPE, TypeIntStruct<float>>::value )
        return MPI_FLOAT_INT;
    if constexpr ( std::is_same<TYPE, TypeIntStruct<double>>::value )
        return MPI_DOUBLE_INT;
    return MPI_DATATYPE_NULL;
}
#endif


/************************************************************************
 *  sumReduce                                                            *
 ************************************************************************/
template<class TYPE>
void MPI_CLASS::sumReduce( const TYPE *send, TYPE *recv, int n ) const
{
    if ( d_size == 1 ) {
        for ( int i = 0; i < n; i++ )
            recv[i] = send[i];
        return;
    }
#ifdef USE_MPI
    PROFILE_START( "sumReduce", profile_level );
    if constexpr ( isMPItype<TYPE>() ) {
        auto type = MPI_getType<TYPE>();
        MPI_CLASS_ASSERT( d_comm != MPI_COMM_NULL );
        MPI_Allreduce( send, recv, n, type, MPI_SUM, d_comm );
    } else {
        throw std::logic_error( "sumReduce is not valid for arbitrary types" );
    }
    PROFILE_STOP( "sumReduce", profile_level );
#endif
}
template<class TYPE>
TYPE MPI_CLASS::sumReduce( const TYPE value ) const
{
    if ( d_size > 1 ) {
        TYPE tmp = value;
        sumReduce( &value, &tmp, 1 );
        return tmp;
    } else {
        return value;
    }
}
template<class TYPE>
void MPI_CLASS::sumReduce( TYPE *x, int n ) const
{
    if ( d_size > 1 ) {
        auto tmp = new TYPE[n];
        for ( int i = 0; i < n; i++ )
            tmp[i] = x[i];
        sumReduce( tmp, x, n );
        delete[] tmp;
    }
}


/************************************************************************
 *  minReduce                                                            *
 ************************************************************************/
template<class TYPE>
void MPI_CLASS::minReduce( const TYPE *send, TYPE *recv, int n ) const
{
    if ( d_size == 1 ) {
        for ( int i = 0; i < n; i++ )
            recv[i] = send[i];
        return;
    }
#ifdef USE_MPI
    PROFILE_START( "minReduce", profile_level );
    if constexpr ( isMPItype<TYPE>() ) {
        auto type = MPI_getType<TYPE>();
        MPI_Allreduce( send, recv, n, type, MPI_MIN, d_comm );
    } else {
        throw std::logic_error( "minReduce is not valid for arbitrary types" );
    }
    PROFILE_STOP( "minReduce", profile_level );
#endif
}
template<class TYPE>
void MPI_CLASS::minReduce( const TYPE *send, TYPE *recv, int n, int *rank ) const
{
    if ( rank == nullptr )
        return minReduce( send, recv, n );
    if ( d_size == 1 ) {
        for ( int i = 0; i < n; i++ ) {
            recv[i] = send[i];
            rank[i] = 0;
        }
        return;
    }
#ifdef USE_MPI
    using TYPE2 = typename TypeIntSelector<TYPE>::type;
    if constexpr ( !isMPItype<TYPE>() )
        throw std::logic_error( "minReduce with rank is not supported for arbitrary types" );
    PROFILE_START( "minReduce2", profile_level );
    auto send2 = new TYPE2[n];
    auto recv2 = new TYPE2[n];
    for ( int i = 0; i < n; i++ ) {
        send2[i].x = send[i];
        send2[i].i = d_rank;
    }
    auto type2 = MPI_getType<TYPE2>();
    MPI_Allreduce( send2, recv2, n, type2, MPI_MINLOC, d_comm );
    for ( int i = 0; i < n; i++ ) {
        recv[i] = recv2[i].x;
        rank[i] = recv2[i].i;
    }
    delete[] send2;
    delete[] recv2;
    PROFILE_STOP( "minReduce2", profile_level );
#endif
}
template<class TYPE>
TYPE MPI_CLASS::minReduce( const TYPE value ) const
{
    if ( d_size > 1 ) {
        TYPE tmp = value;
        minReduce( &value, &tmp, 1 );
        return tmp;
    } else {
        return value;
    }
}
template<class TYPE>
void MPI_CLASS::minReduce( TYPE *x, int n ) const
{
    if ( d_size > 1 ) {
        auto tmp = new TYPE[n];
        for ( int i = 0; i < n; i++ )
            tmp[i] = x[i];
        minReduce( tmp, x, n );
        delete[] tmp;
    }
}
template<class TYPE>
void MPI_CLASS::minReduce( TYPE *x, int n, int *rank ) const
{
    if ( rank == nullptr )
        return minReduce( x, n );
    if ( d_size > 1 ) {
        auto tmp = new TYPE[n];
        for ( int i = 0; i < n; i++ )
            tmp[i] = x[i];
        minReduce( tmp, x, n, rank );
        delete[] tmp;
    } else {
        rank[0] = 0;
    }
}


/************************************************************************
 *  maxReduce                                                            *
 ************************************************************************/
template<class TYPE>
void MPI_CLASS::maxReduce( const TYPE *send, TYPE *recv, int n ) const
{
    if ( d_size == 1 ) {
        for ( int i = 0; i < n; i++ )
            recv[i] = send[i];
        return;
    }
#ifdef USE_MPI
    PROFILE_START( "maxReduce", profile_level );
    if constexpr ( isMPItype<TYPE>() ) {
        auto type = MPI_getType<TYPE>();
        MPI_Allreduce( send, recv, n, type, MPI_MAX, d_comm );
    } else {
        throw std::logic_error( "maxReduce is not valid for arbitrary types" );
    }
    PROFILE_STOP( "maxReduce", profile_level );
#endif
}
template<class TYPE>
void MPI_CLASS::maxReduce( const TYPE *send, TYPE *recv, int n, int *rank ) const
{
    if ( rank == nullptr )
        return maxReduce( send, recv, n );
    if ( d_size == 1 ) {
        for ( int i = 0; i < n; i++ ) {
            recv[i] = send[i];
            rank[i] = 0;
        }
        return;
    }
#ifdef USE_MPI
    using TYPE2 = typename TypeIntSelector<TYPE>::type;
    if constexpr ( !isMPItype<TYPE>() )
        throw std::logic_error( "maxReduce with rank is not supported for arbitrary types" );
    PROFILE_START( "maxReduce2", profile_level );
    auto send2 = new TYPE2[n];
    auto recv2 = new TYPE2[n];
    for ( int i = 0; i < n; i++ ) {
        send2[i].x = send[i];
        send2[i].i = d_rank;
    }
    auto type2 = MPI_getType<TYPE2>();
    MPI_Allreduce( send2, recv2, n, type2, MPI_MAXLOC, d_comm );
    for ( int i = 0; i < n; i++ ) {
        recv[i] = recv2[i].x;
        rank[i] = recv2[i].i;
    }
    delete[] send2;
    delete[] recv2;
    PROFILE_STOP( "maxReduce2", profile_level );
#endif
}
template<class TYPE>
TYPE MPI_CLASS::maxReduce( const TYPE value ) const
{
    if ( d_size > 1 ) {
        TYPE tmp = value;
        maxReduce( &value, &tmp, 1 );
        return tmp;
    } else {
        return value;
    }
}
template<class TYPE>
void MPI_CLASS::maxReduce( TYPE *x, int n ) const
{
    if ( d_size > 1 ) {
        auto tmp = new TYPE[n];
        for ( int i = 0; i < n; i++ )
            tmp[i] = x[i];
        maxReduce( tmp, x, n );
        delete[] tmp;
    }
}
template<class TYPE>
void MPI_CLASS::maxReduce( TYPE *x, int n, int *rank ) const
{
    if ( rank == nullptr )
        return minReduce( x, n );
    if ( d_size > 1 ) {
        auto tmp = new TYPE[n];
        for ( int i = 0; i < n; i++ )
            tmp[i] = x[i];
        maxReduce( tmp, x, n, rank );
        delete[] tmp;
    } else {
        rank[0] = 0;
    }
}


/************************************************************************
 *  bcast                                                                *
 ************************************************************************/
#ifdef USE_MPI
template<class TYPE>
void MPI_CLASS::bcast( TYPE *x, int n, int root ) const
{
    if ( d_size == 1 )
        return;
    if ( root >= d_size )
        MPI_CLASS_ERROR( "root cannot be >= size in bcast" );
    PROFILE_SCOPED( timer, "bcast", profile_level );
    if constexpr ( isMPItype<TYPE>() ) {
        auto type = MPI_getType<TYPE>();
        MPI_Bcast( x, n, type, root, d_comm );
    } else if constexpr ( is_mpi_copyable<TYPE>::value ) {
        MPI_Bcast( x, n * sizeof( TYPE ), MPI_BYTE, root, d_comm );
    } else if constexpr ( is_vector<TYPE>::value ) {
        for ( int i = 0; i < n; i++ ) {
            size_t N = bcast( x[i].size(), root );
            x[i].resize( N );
            bcast( x[i].data(), N, root );
        }
    } else {
        call_bcast<TYPE>( x, n, root );
    }
}
#else
template<class TYPE>
void MPI_CLASS::bcast( TYPE *, int, int ) const
{
}
#endif
// Default implimentation of bcast
template<class TYPE>
TYPE MPI_CLASS::bcast( TYPE x, int root ) const
{
    if ( d_size > 1 )
        bcast( &x, 1, root );
    return x;
}


/************************************************************************
 *  send                                                                 *
 ************************************************************************/
#ifdef USE_MPI
template<class TYPE>
void MPI_CLASS::send( const TYPE *buf, int length, int recv_proc, int tag ) const
{
    // Set the tag to 0 if it is < 0
    tag = ( tag >= 0 ) ? tag : 0;
    MPI_CLASS_INSIST( tag <= d_maxTag, "Maximum tag value exceeded" );
    // Send the data
    PROFILE_START( "send", profile_level );
    if constexpr ( isMPItype<TYPE>() ) {
        auto type = MPI_getType<TYPE>();
        MPI_Send( buf, length, type, recv_proc, tag, d_comm );
    } else if constexpr ( is_mpi_copyable<TYPE>::value ) {
        size_t bytes = length * sizeof( TYPE );
        MPI_Send( buf, bytes, MPI_BYTE, recv_proc, tag, d_comm );
    } else {
        MPI_CLASS_ERROR_TYPE( "Arbitrary types are not supported" );
    }
    PROFILE_STOP( "send", profile_level );
}
#else
template<class TYPE>
void MPI_CLASS::send( const TYPE *buf, int length, int recv_proc, int tag ) const
{
    int bytes = length * sizeof( TYPE );
    sendBytes( buf, bytes, recv_proc, tag );
}
#endif


/************************************************************************
 *  Isend                                                                *
 ************************************************************************/
#ifdef USE_MPI
template<class TYPE>
MPI_CLASS::Request MPI_CLASS::Isend( const TYPE *buf, int length, int recv_proc, int tag ) const
{
    // Set the tag to 0 if it is < 0
    tag = ( tag >= 0 ) ? tag : 0;
    MPI_CLASS_INSIST( tag <= d_maxTag, "Maximum tag value exceeded" );
    // Send the data
    PROFILE_START( "isend", profile_level );
    MPI_CLASS::Request request;
    if constexpr ( isMPItype<TYPE>() ) {
        auto type = MPI_getType<TYPE>();
        MPI_Isend( buf, length, type, recv_proc, tag, d_comm, &request );
    } else if constexpr ( is_mpi_copyable<TYPE>::value ) {
        size_t bytes = length * sizeof( TYPE );
        MPI_Isend( buf, bytes, MPI_BYTE, recv_proc, tag, d_comm, &request );
    } else {
        MPI_CLASS_ERROR_TYPE( "Arbitrary types are not supported" );
    }
    PROFILE_STOP( "isend", profile_level );
    return request;
}
#else
template<class TYPE>
MPI_CLASS::Request MPI_CLASS::Isend( const TYPE *buf, int length, int recv_proc, int tag ) const
{
    int bytes = length * sizeof( TYPE );
    return IsendBytes( buf, bytes, recv_proc, tag );
}
#endif


/************************************************************************
 *  recv                                                                 *
 ************************************************************************/
#ifdef USE_MPI
template<class TYPE>
void MPI_CLASS::recv( TYPE *buf, int &length, int send_proc, bool get_length, int tag ) const
{
    // Set the tag to 0 if it is < 0
    tag = ( tag >= 0 ) ? tag : 0;
    MPI_CLASS_INSIST( tag <= d_maxTag, "Maximum tag value exceeded" );
    PROFILE_START( "recv", profile_level );
    // Get the receive length if necessary
    if ( get_length ) {
        int bytes       = this->probe( send_proc, tag );
        int recv_length = bytes / sizeof( TYPE );
        MPI_CLASS_INSIST( length >= recv_length, "Recived length is larger than allocated array" );
        length = recv_length;
    }
    // Recieve the data
    MPI_Status status;
    if constexpr ( isMPItype<TYPE>() ) {
        auto type = MPI_getType<TYPE>();
        MPI_Recv( buf, length, type, send_proc, tag, d_comm, &status );
    } else if constexpr ( is_mpi_copyable<TYPE>::value ) {
        size_t bytes = length * sizeof( TYPE );
        MPI_Recv( buf, bytes, MPI_BYTE, send_proc, tag, d_comm, &status );
    } else {
        MPI_CLASS_ERROR_TYPE( "Arbitrary types are not supported" );
    }
    PROFILE_STOP( "recv", profile_level );
}
#else
template<class TYPE>
void MPI_CLASS::recv( TYPE *buf, int &length, int send_proc, bool get_length, int tag ) const
{
    // Get the receive length if necessary
    if ( get_length ) {
        int bytes       = this->probe( send_proc, tag );
        int recv_length = bytes / sizeof( TYPE );
        MPI_CLASS_INSIST( length >= recv_length, "Recived length is larger than allocated array" );
        length = recv_length;
    }
    // Recieve the data
    int bytes = length * sizeof( TYPE );
    recvBytes( buf, bytes, send_proc, tag );
}
#endif
template<class TYPE>
void MPI_CLASS::recv( TYPE *buf, int length, int send_proc, int tag ) const
{
    recv<TYPE>( buf, length, send_proc, false, tag );
}


/************************************************************************
 *  Irecv                                                                *
 ************************************************************************/
// Define specializations of recv(TYPE*, int&, int, bool, int)
#ifdef USE_MPI
template<class TYPE>
MPI_CLASS::Request MPI_CLASS::Irecv( TYPE *buf, int length, int send_proc, int tag ) const
{
    MPI_CLASS_INSIST( tag <= d_maxTag, "Maximum tag value exceeded" );
    MPI_CLASS_INSIST( tag >= 0, "tag must be >= 0" );
    MPI_CLASS::Request request;
    PROFILE_START( "Irecv", profile_level );
    if constexpr ( isMPItype<TYPE>() ) {
        auto type = MPI_getType<TYPE>();
        MPI_Irecv( buf, length, type, send_proc, tag, d_comm, &request );
    } else if constexpr ( is_mpi_copyable<TYPE>::value ) {
        size_t bytes = length * sizeof( TYPE );
        MPI_Irecv( buf, bytes, MPI_BYTE, send_proc, tag, d_comm, &request );
    } else {
        MPI_CLASS_ERROR_TYPE( "Arbitrary types are not supported" );
    }
    PROFILE_STOP( "Irecv", profile_level );
    return request;
}
#else
template<class TYPE>
MPI_CLASS::Request MPI_CLASS::Irecv( TYPE *buf, int length, int send_proc, int tag ) const
{
    int bytes = length * sizeof( TYPE );
    return IrecvBytes( buf, bytes, send_proc, tag );
}
#endif


/************************************************************************
 *  sendrecv                                                             *
 ************************************************************************/
template<class TYPE>
void MPI_CLASS::sendrecv( const TYPE *sendbuf,
                          int sendcount,
                          int dest,
                          int sendtag,
                          TYPE *recvbuf,
                          int recvcount,
                          int source,
                          int recvtag ) const
{
    if ( getSize() <= 1 ) {
        MPI_CLASS_ASSERT( dest == 0 );
        MPI_CLASS_ASSERT( source == 0 );
        MPI_CLASS_ASSERT( sendcount <= recvcount );
        MPI_CLASS_ASSERT( sendtag == recvtag );
        for ( int i = 0; i < sendcount; i++ )
            recvbuf[i] = sendbuf[i];
        return;
    }
#ifdef USE_MPI
    PROFILE_START( "sendrecv", profile_level );
    if constexpr ( isMPItype<TYPE>() ) {
        auto type = MPI_getType<TYPE>();
        MPI_Sendrecv( sendbuf,
                      sendcount,
                      type,
                      dest,
                      sendtag,
                      recvbuf,
                      recvcount,
                      type,
                      source,
                      recvtag,
                      d_comm,
                      MPI_STATUS_IGNORE );
    } else if constexpr ( is_mpi_copyable<TYPE>::value ) {
        size_t sendbytes = sendcount * sizeof( TYPE );
        size_t recvbytes = recvcount * sizeof( TYPE );
        MPI_Sendrecv( sendbuf,
                      sendbytes,
                      MPI_BYTE,
                      dest,
                      sendtag,
                      recvbuf,
                      recvbytes,
                      MPI_BYTE,
                      source,
                      recvtag,
                      d_comm,
                      MPI_STATUS_IGNORE );
    } else {
        MPI_CLASS_ERROR_TYPE( "Arbitrary types are not supported" );
    }
    PROFILE_STOP( "sendrecv", profile_level );
#endif
}


/************************************************************************
 *  gather                                                            *
 ************************************************************************/
template<class TYPE>
std::vector<TYPE> MPI_CLASS::gather( const TYPE &x, int root ) const
{
    if ( d_size == 1 )
        return { x };
    std::vector<int> recv_cnt, recv_disp;
    if ( d_rank == root ) {
        recv_cnt.resize( d_size, 1 );
        recv_disp.resize( d_size, 0 );
        for ( int i = 0; i < d_size; i++ )
            recv_disp[i] = i;
    }
    std::vector<TYPE> y( d_size );
    gather( &x, 1, y.data(), recv_cnt.data(), recv_disp.data(), root );
    return y;
}
template<class TYPE>
std::vector<TYPE> MPI_CLASS::gather( const std::vector<TYPE> &x, int root ) const
{
    if ( d_size == 1 )
        return x;
    auto recv_cnt = gather<int>( x.size(), root );
    std::vector<int> recv_disp( recv_cnt.size(), 0 );
    for ( size_t i = 1; i < recv_cnt.size(); i++ )
        recv_disp[i] = recv_disp[i - 1] + recv_cnt[i - 1];
    std::vector<TYPE> y;
    if ( d_rank == root ) {
        size_t N = 0;
        for ( auto s : recv_cnt )
            N += s;
        y.resize( N );
    }
    gather<TYPE>( x.data(), x.size(), y.data(), recv_cnt.data(), recv_disp.data(), root );
    return y;
}
template<class TYPE>
void MPI_CLASS::gather( const TYPE *sendbuf,
                        int sendcount,
                        TYPE *recvbuf,
                        const int *recv_cnt,
                        const int *recv_disp,
                        int root ) const
{
    NULL_USE( root );
    if ( d_size == 1 ) {
        MPI_CLASS_ASSERT( sendcount == recv_cnt[0] );
        for ( int i = 0; i < sendcount; i++ )
            recvbuf[i + recv_disp[0]] = sendbuf[i];
        return;
    }
#ifdef USE_MPI
    PROFILE_START( "gather", profile_level );
    if constexpr ( isMPItype<TYPE>() ) {
        auto type = MPI_getType<TYPE>();
        MPI_Gatherv( sendbuf, sendcount, type, recvbuf, recv_cnt, recv_disp, type, root, d_comm );
    } else if constexpr ( is_mpi_copyable<TYPE>::value ) {
        size_t N = sizeof( TYPE );
        std::vector<int> recv_cnt2( d_size );
        std::vector<int> recv_disp2( d_size );
        for ( int i = 0; i < d_size; i++ ) {
            recv_cnt2[i]  = N * recv_cnt[i];
            recv_disp2[i] = N * recv_disp[i];
        }
        MPI_Gatherv( sendbuf,
                     N * sendcount,
                     MPI_BYTE,
                     recvbuf,
                     recv_cnt2.data(),
                     recv_disp2.data(),
                     MPI_BYTE,
                     root,
                     d_comm );
    } else {
        call_gatherv( sendbuf, sendcount, recvbuf, recv_cnt, recv_disp, root );
    }
    PROFILE_STOP( "gather", profile_level );
#endif
}


/************************************************************************
 *  allGather                                                            *
 ************************************************************************/
template<class TYPE>
std::vector<TYPE> MPI_CLASS::allGather( const TYPE &x ) const
{
    std::vector<TYPE> data( getSize() );
    allGather( x, data.data() );
    return data;
}
template<class TYPE>
std::vector<TYPE> MPI_CLASS::allGather( const std::vector<TYPE> &x ) const
{
    static_assert( is_mpi_copyable<TYPE>::value, "Object is not trivially copyable" );
    if ( getSize() <= 1 )
        return x;
    std::vector<int> count = allGather<int>( x.size() );
    std::vector<int> disp( getSize(), 0 );
    size_t N = count[0];
    for ( size_t i = 1; i < count.size(); i++ ) {
        disp[i] = disp[i - 1] + count[i - 1];
        N += count[i];
    }
    std::vector<TYPE> data( N );
    allGather<TYPE>( x.data(), x.size(), data.data(), count.data(), disp.data(), true );
    return data;
}
// Default instantiation of MPI_CLASS::allGather
template<class TYPE>
void MPI_CLASS::allGather( const TYPE &x, TYPE *y ) const
{
    MPI_CLASS_ASSERT( y );
    if ( d_size == 1 ) {
        y[0] = x;
        return;
    }
#ifdef USE_MPI
    PROFILE_START( "allGather", profile_level );
    if constexpr ( isMPItype<TYPE>() ) {
        auto type = MPI_getType<TYPE>();
        MPI_Allgather( &x, 1, type, y, 1, type, d_comm );
    } else if constexpr ( is_mpi_copyable<TYPE>::value ) {
        size_t N = sizeof( TYPE );
        MPI_Allgather( &x, N, MPI_BYTE, y, N, MPI_BYTE, d_comm );
    } else {
        call_allGather( x, y );
    }
    PROFILE_STOP( "allGather", profile_level );
#endif
}
template<class TYPE>
int MPI_CLASS::allGather( const TYPE *send_data,
                          int send_cnt,
                          TYPE *recv_data,
                          int *recv_cnt,
                          int *recv_disp,
                          bool known_recv ) const
{
    if ( d_size == 1 ) {
        if ( !known_recv ) {
            for ( int i = 0; i < send_cnt; i++ )
                recv_data[i] = send_data[i];
            if ( recv_cnt )
                *recv_cnt = send_cnt;
            if ( recv_disp )
                *recv_disp = 0;
        } else {
            for ( int i = 0; i < send_cnt; i++ )
                recv_data[i + recv_disp[0]] = send_data[i];
        }
        return send_cnt;
    }
#ifdef USE_MPI
    PROFILE_START( "allGather2", profile_level );
    // Get the sizes of the received data (if necessary)
    int *recv_cnt2  = recv_cnt;
    int *recv_disp2 = recv_disp;
    if ( !known_recv ) {
        if ( !recv_cnt )
            recv_cnt2 = new int[d_size];
        if ( !recv_disp )
            recv_disp2 = new int[d_size];
        allGather( send_cnt, recv_cnt2 );
        recv_disp2[0] = 0;
        for ( int i = 1; i < d_size; i++ )
            recv_disp2[i] = recv_disp2[i - 1] + recv_cnt2[i - 1];
    }
    int N_recv = 0;
    for ( int i = 0; i < d_size; i++ )
        N_recv += recv_cnt2[i];
    // Send/recv the data
    if constexpr ( isMPItype<TYPE>() ) {
        auto type = MPI_getType<TYPE>();
        MPI_Allgatherv( send_data, send_cnt, type, recv_data, recv_cnt2, recv_disp2, type, d_comm );
    } else if constexpr ( is_mpi_copyable<TYPE>::value ) {
        auto cnt  = new int[d_size];
        auto disp = new int[d_size];
        for ( int i = 0; i < d_size; i++ ) {
            cnt[i]  = sizeof( TYPE ) * recv_cnt2[i];
            disp[i] = sizeof( TYPE ) * recv_disp2[i];
        }
        int N = send_cnt * sizeof( TYPE );
        MPI_Allgatherv( send_data, N, MPI_BYTE, recv_data, cnt, disp, MPI_BYTE, d_comm );
        delete[] cnt;
        delete[] disp;
    } else {
        call_allGather( send_data, send_cnt, recv_data, recv_cnt2, recv_disp2 );
    }
    // Delete any temporary memory
    if ( !recv_cnt )
        delete[] recv_cnt2;
    if ( !recv_disp )
        delete[] recv_disp2;
    PROFILE_STOP( "allGather2", profile_level );
    return N_recv;
#else
    return 0;
#endif
}


/************************************************************************
 *  setGather                                                            *
 ************************************************************************/
template<class TYPE>
void MPI_CLASS::setGather( std::set<TYPE> &set ) const
{
    std::vector<TYPE> send_buf( set.begin(), set.end() );
    std::vector<int> recv_cnt( this->d_size, 0 );
    this->allGather<int>( (int) send_buf.size(), &recv_cnt[0] );
    std::vector<int> recv_disp( this->d_size, 0 );
    for ( int i = 1; i < this->d_size; i++ )
        recv_disp[i] = recv_disp[i - 1] + recv_cnt[i - 1];
    size_t N_recv_tot = 0;
    for ( int i = 0; i < this->d_size; i++ )
        N_recv_tot += recv_cnt[i];
    if ( N_recv_tot == 0 )
        return;
    std::vector<TYPE> recv_buf( N_recv_tot );
    TYPE *send_data = nullptr;
    if ( send_buf.size() > 0 ) {
        send_data = &send_buf[0];
    }
    TYPE *recv_data = &recv_buf[0];
    static_assert( is_mpi_copyable<TYPE>::value, "Object is not trivially copyable" );
    this->allGather<TYPE>(
        send_data, (int) send_buf.size(), recv_data, &recv_cnt[0], &recv_disp[0], true );
    for ( size_t i = 0; i < recv_buf.size(); i++ )
        set.insert( recv_buf[i] );
}


/************************************************************************
 *  mapGather                                                            *
 ************************************************************************/
template<class KEY, class DATA>
void MPI_CLASS::mapGather( std::map<KEY, DATA> &map ) const
{
    std::vector<KEY> send_id;
    std::vector<DATA> send_data;
    send_id.reserve( map.size() );
    send_data.reserve( map.size() );
    for ( auto it = map.begin(); it != map.end(); ++it ) {
        send_id.push_back( it->first );
        send_data.push_back( it->second );
    }
    int send_size = (int) send_id.size();
    std::vector<int> recv_cnt( this->d_size, 0 );
    this->allGather<int>( send_size, &recv_cnt[0] );
    std::vector<int> recv_disp( this->d_size, 0 );
    for ( int i = 1; i < this->d_size; i++ )
        recv_disp[i] = recv_disp[i - 1] + recv_cnt[i - 1];
    size_t N_recv_tot = 0;
    for ( int i = 0; i < this->d_size; i++ )
        N_recv_tot += recv_cnt[i];
    if ( N_recv_tot == 0 )
        return;
    std::vector<KEY> recv_id( N_recv_tot );
    std::vector<DATA> recv_data( N_recv_tot );
    KEY *send_data1  = nullptr;
    DATA *send_data2 = nullptr;
    if ( send_id.size() > 0 ) {
        send_data1 = &send_id[0];
        send_data2 = &send_data[0];
    }
    static_assert( is_mpi_copyable<DATA>::value, "Object is not trivially copyable" );
    this->allGather<KEY>( send_data1, send_size, &recv_id[0], &recv_cnt[0], &recv_disp[0], true );
    this->allGather<DATA>(
        send_data2, send_size, &recv_data[0], &recv_cnt[0], &recv_disp[0], true );
    map = std::map<KEY, DATA>();
    for ( size_t i = 0; i < N_recv_tot; i++ )
        map.insert( std::pair<KEY, DATA>( recv_id[i], recv_data[i] ) );
}


/************************************************************************
 *  sumScan                                                              *
 ************************************************************************/
template<class TYPE>
TYPE MPI_CLASS::sumScan( const TYPE &x ) const
{
    TYPE y = x;
    sumScan<TYPE>( &x, &y, 1 );
    return y;
}
template<class TYPE>
void MPI_CLASS::sumScan( const TYPE *x, TYPE *y, int n ) const
{
    if ( d_size == 1 ) {
        for ( int i = 0; i < n; i++ )
            y[i] = x[i];
        return;
    }
#ifdef USE_MPI
    PROFILE_START( "sumScan", profile_level );
    if constexpr ( isMPItype<TYPE>() ) {
        auto type = MPI_getType<TYPE>();
        MPI_CLASS_ASSERT( d_comm != MPI_COMM_NULL );
        MPI_Scan( x, y, n, type, MPI_SUM, d_comm );
    } else {
        throw std::logic_error( "Default instantion of sumScan in parallel is not supported" );
    }
    PROFILE_STOP( "sumScan", profile_level );
#endif
}


/************************************************************************
 *  minScan                                                              *
 ************************************************************************/
template<class TYPE>
TYPE MPI_CLASS::minScan( const TYPE &x ) const
{
    TYPE y = x;
    minScan<TYPE>( &x, &y, 1 );
    return y;
}
template<class TYPE>
void MPI_CLASS::minScan( const TYPE *x, TYPE *y, int n ) const
{
    if ( d_size == 1 ) {
        for ( int i = 0; i < n; i++ )
            y[i] = x[i];
        return;
    }
#ifdef USE_MPI
    PROFILE_START( "minScan", profile_level );
    if constexpr ( isMPItype<TYPE>() ) {
        auto type = MPI_getType<TYPE>();
        MPI_CLASS_ASSERT( d_comm != MPI_COMM_NULL );
        MPI_Scan( x, y, n, type, MPI_MIN, d_comm );
    } else {
        throw std::logic_error( "Default instantion of minScan in parallel is not supported" );
    }
    PROFILE_STOP( "minScan", profile_level );
#endif
}


/************************************************************************
 *  maxScan                                                              *
 ************************************************************************/
template<class TYPE>
TYPE MPI_CLASS::maxScan( const TYPE &x ) const
{
    TYPE y = x;
    maxScan<TYPE>( &x, &y, 1 );
    return y;
}
template<class TYPE>
void MPI_CLASS::maxScan( const TYPE *x, TYPE *y, int n ) const
{
    if ( d_size == 1 ) {
        for ( int i = 0; i < n; i++ )
            y[i] = x[i];
        return;
    }
#ifdef USE_MPI
    PROFILE_START( "maxScan", profile_level );
    if constexpr ( isMPItype<TYPE>() ) {
        auto type = MPI_getType<TYPE>();
        MPI_CLASS_ASSERT( d_comm != MPI_COMM_NULL );
        MPI_Scan( x, y, n, type, MPI_MAX, d_comm );
    } else {
        throw std::logic_error( "Default instantion of maxScan in parallel is not supported" );
    }
    PROFILE_STOP( "maxScan", profile_level );
#endif
}


/************************************************************************
 *  allToAll                                                             *
 ************************************************************************/
template<class TYPE>
void MPI_CLASS::allToAll( int n, const TYPE *send, TYPE *recv ) const
{
    if ( d_size == 1 ) {
        for ( int i = 0; i < n; i++ )
            recv[i] = send[i];
        return;
    }
#ifdef USE_MPI
    if constexpr ( isMPItype<TYPE>() ) {
        auto type = MPI_getType<TYPE>();
        MPI_Alltoall( send, n, type, recv, n, type, d_comm );
    } else if constexpr ( is_mpi_copyable<TYPE>::value ) {
        int N = n * sizeof( TYPE );
        MPI_Alltoall( send, N, MPI_BYTE, recv, N, MPI_BYTE, d_comm );
    } else {
        MPI_CLASS_ERROR_TYPE( "Arbitrary types are not supported" );
    }
#endif
}


/************************************************************************
 *  allToAll                                                             *
 ************************************************************************/
template<class TYPE>
int MPI_CLASS::allToAll( const TYPE *send_data,
                         const int send_cnt[],
                         const int send_disp[],
                         TYPE *recv_data,
                         int *recv_cnt,
                         int *recv_disp,
                         bool known_recv ) const
{
    if ( d_size == 1 ) {
        // Special case for single-processor communicators
        if ( known_recv ) {
            if ( recv_cnt[0] != send_cnt[0] && send_cnt[0] > 0 )
                MPI_CLASS_ERROR( "Single processor send/recv are different sizes" );
        } else {
            if ( recv_cnt != nullptr )
                recv_cnt[0] = send_cnt[0];
            if ( recv_disp != nullptr )
                recv_disp[0] = send_disp[0];
        }
        for ( int i = 0; i < send_cnt[0]; i++ )
            recv_data[i + recv_disp[0]] = send_data[i + send_disp[0]];
        return send_cnt[0];
    }
#ifdef USE_MPI
    PROFILE_START( "allToAll", profile_level );
    // Get the recv sizes if needed
    int *recv_cnt2  = recv_cnt;
    int *recv_disp2 = recv_disp;
    if ( !known_recv ) {
        // Communicate the size we will be recieving from each processor
        if ( !recv_cnt )
            recv_cnt2 = new int[d_size];
        if ( !recv_disp )
            recv_disp2 = new int[d_size];
        allToAll<int>( 1, send_cnt, recv_cnt2 );
        recv_disp2[0] = 0;
        for ( int i = 1; i < d_size; i++ )
            recv_disp2[i] = recv_disp2[i - 1] + recv_cnt2[i - 1];
    }
    // Perform the communication
    MPI_CLASS_ASSERT( recv_cnt2 && recv_disp2 );
    if constexpr ( isMPItype<TYPE>() ) {
        auto type = MPI_getType<TYPE>();
        MPI_Alltoallv(
            send_data, send_cnt, send_disp, type, recv_data, recv_cnt2, recv_disp2, type, d_comm );
    } else if constexpr ( is_mpi_copyable<TYPE>::value ) {
        auto send_cnt3  = new int[d_size];
        auto send_disp3 = new int[d_size];
        auto recv_cnt3  = new int[d_size];
        auto recv_disp3 = new int[d_size];
        for ( int i = 0; i < d_size; i++ ) {
            send_cnt3[i]  = sizeof( TYPE ) * send_cnt[i];
            send_disp3[i] = sizeof( TYPE ) * send_disp[i];
            recv_cnt3[i]  = sizeof( TYPE ) * recv_cnt2[i];
            recv_disp3[i] = sizeof( TYPE ) * recv_disp2[i];
        }
        MPI_Alltoallv( send_data,
                       send_cnt3,
                       send_disp3,
                       MPI_BYTE,
                       recv_data,
                       recv_cnt3,
                       recv_disp3,
                       MPI_BYTE,
                       d_comm );
        delete[] send_cnt3;
        delete[] send_disp3;
        delete[] recv_cnt3;
        delete[] recv_disp3;
    } else {
        MPI_CLASS_ERROR( "Not finished" );
    }
    int N_received = 0;
    for ( int i = 0; i < d_size; i++ )
        N_received += recv_cnt2[i];
    // Delete the temporary memory if needed
    if ( !recv_cnt )
        delete[] recv_cnt2;
    if ( !recv_disp )
        delete[] recv_disp2;
    PROFILE_STOP( "allToAll", profile_level );
    return N_received;
#else
    return 0;
#endif
}


} // namespace AMP


#endif
