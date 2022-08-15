#ifndef included_AMP_RGBA
#define included_AMP_RGBA

#include <cstdint>


namespace AMP {


//! Structure to store RGBA data
struct RGBA {
    uint8_t red;
    uint8_t green;
    uint8_t blue;
    uint8_t alpha;
    RGBA() : red( 0 ), green( 0 ), blue( 0 ), alpha( 255 ){};
    RGBA( uint8_t x ) : red( x ), green( x ), blue( x ), alpha( 255 ){};
    RGBA( uint8_t r, uint8_t g, uint8_t b, uint8_t a = 255 )
        : red( r ), green( g ), blue( b ), alpha( a ){};
    operator uint8_t() { return red * 0.299 + green * 0.587 + blue * 0.114; }
    operator uint16_t() { return red * 76.843 + green * 150.859 + blue * 29.298; }
    operator float()
    {
        return red * 0.001172549019608 + green * 0.002301960784314 + blue * 0.000447058823529;
    }
    operator double()
    {
        return red * 0.001172549019608 + green * 0.002301960784314 + blue * 0.000447058823529;
    }
    bool operator==( const RGBA &rhs ) const
    {
        return red == rhs.red && green == rhs.green && blue == rhs.blue && alpha == rhs.alpha;
    }
};


} // namespace AMP

#endif
