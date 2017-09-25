// this is just meant to experiment with how best to include
// templated kernels

void boo( void );

template <typename T>
class A
{
public:
    void foo( void ) { boo(); }
};
