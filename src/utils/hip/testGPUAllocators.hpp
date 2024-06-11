#include <stdio.h>

template<typename T>
void opDataW( T *d_A, size_t n );

template<typename T>
void setDataW( T *d_A, T v, size_t n );


template<typename T>
class KernelWrapper
{
public:
    void opData( T *d_A, size_t n ) { opDataW<T>( d_A, n ); }
    void setData( T *d_A, T v, size_t n ) { setDataW<T>( d_A, v, n ); }
};
