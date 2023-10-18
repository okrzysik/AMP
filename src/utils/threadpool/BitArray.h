#ifndef included_AMP_BitArray
#define included_AMP_BitArray

#include <atomic>
#include <vector>


namespace AMP {


/** \class BitArray
 *
 * \brief This class stores arrays of bits
 * \details This class stores arrays of bits.  It supports setting and unsetting bits
 *   in a thread-safe way using atomic operations.
 */
class BitArray final
{
public:
    BitArray() : d_N( 0 ), d_data( nullptr ) {}
    explicit BitArray( size_t N ) : d_N( N ), d_data( nullptr )
    {
        size_t N2 = ( d_N + 63 ) >> 6;
        d_data    = new std::atomic_uint64_t[N2];
        for ( size_t i = 0; i < N2; i++ )
            d_data[i] = 0;
    }
    BitArray( const BitArray &rhs ) : d_N( rhs.d_N ), d_data( nullptr )
    {
        size_t N2 = ( d_N + 63 ) >> 6;
        d_data    = new std::atomic_uint64_t[N2];
        for ( size_t i = 0; i < N2; i++ )
            d_data[i] = rhs.d_data[i].load();
    }
    BitArray( BitArray &&rhs ) : d_N( rhs.d_N ), d_data( rhs.d_data ) { rhs.d_data = nullptr; }
    BitArray &operator=( const BitArray &rhs )
    {
        if ( &rhs == this )
            return *this;
        d_N       = rhs.d_N;
        size_t N2 = ( d_N + 63 ) >> 6;
        d_data    = new std::atomic_uint64_t[N2];
        for ( size_t i = 0; i < N2; i++ )
            d_data[i] = rhs.d_data[i].load();
    }
    BitArray &operator=( BitArray &&rhs )
    {
        if ( &rhs == this )
            return *this;
        d_N        = rhs.d_N;
        d_data     = rhs.d_data;
        rhs.d_N    = 0;
        rhs.d_data = nullptr;
    }
    ~BitArray() { delete[] d_data; }
    inline void set( uint64_t index )
    {
        uint64_t mask = ( (uint64_t) 0x01 ) << ( index & 0x3F );
        d_data[index >> 6].fetch_or( mask );
    }
    inline void unset( uint64_t index )
    {
        uint64_t mask = ( (uint64_t) 0x01 ) << ( index & 0x3F );
        d_data[index >> 6].fetch_and( ~mask );
    }
    inline bool get( uint64_t index ) const
    {
        uint64_t mask = ( (uint64_t) 0x01 ) << ( index & 0x3F );
        return ( d_data[index >> 6] & mask ) != 0;
    }
    inline size_t sum() const
    {
        size_t count = 0;
        size_t N2    = ( d_N + 63 ) >> 6;
        for ( size_t i = 0; i < N2; i++ )
            count += popcount64( d_data[i] );
        return count;
    }
    inline std::vector<int> getIndicies() const
    {
        std::vector<int> index( sum() );
        size_t N2 = ( d_N + 63 ) >> 6;
        for ( size_t i = 0, j = 0, k = 0; i < N2; i++ ) {
            uint64_t mask = 0x01;
            for ( size_t m = 0; m < 64; m++, k++, mask <<= 1 ) {
                if ( ( d_data[i] & mask ) != 0 )
                    index[j++] = k;
            }
        }
        return index;
    }
    inline operator std::vector<bool>() const
    {
        size_t N2 = ( d_N + 63 ) >> 6;
        std::vector<bool> x( 64 * N2, false );
        for ( size_t i = 0; i < N2; i++ ) {
            uint64_t tmp = d_data[i];
            for ( size_t j = 0; j < 64; j++ ) {
                uint64_t mask = ( (uint64_t) 0x01 ) << j;
                if ( ( tmp & mask ) != 0 )
                    x[i * 64 + j] = true;
            }
        }
        x.resize( d_N );
        return x;
    }

private:
    static inline size_t popcount64( uint64_t x )
    {
        x = ( x & 0x5555555555555555LU ) + ( x >> 1 & 0x5555555555555555LU );
        x = ( x & 0x3333333333333333LU ) + ( x >> 2 & 0x3333333333333333LU );
        x = ( x + ( x >> 4 ) ) & 0x0F0F0F0F0F0F0F0FLU;
        x = ( x + ( x >> 8 ) );
        x = ( x + ( x >> 16 ) );
        x = ( x + ( x >> 32 ) ) & 0x0000007F;
        return x;
    }

private:
    size_t d_N;
    std::atomic_uint64_t *d_data;
};


} // namespace AMP

#endif
