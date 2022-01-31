#ifndef included_AMP_typeid
#define included_AMP_typeid

#include <complex>
#include <type_traits>


namespace AMP {


//! Class to store type info
struct alignas( 8 ) typeID {
    char name[116] = { 0 };
    uint32_t hash  = { 0 };
    uint64_t bytes = { 0 };
    constexpr bool operator==( size_t rhs ) const { return hash == rhs; }
    constexpr bool operator!=( size_t rhs ) const { return hash != rhs; }
    constexpr bool operator==( const typeID &rhs ) const { return hash == rhs.hash; }
    constexpr bool operator!=( const typeID &rhs ) const { return hash != rhs.hash; }
};


// Helper function to copy a string
constexpr void copy( char *dst, const char *src )
{
    for ( size_t i = 0; src[i] != 0; i++ )
        dst[i] = src[i];
}


//! Get the type info (does not resolve dynamic types)
template<typename T0>
constexpr typeID getTypeID()
{
    typeID id = {};
    // Remove const/references
    using T1 = typename std::remove_reference<T0>::type;
    using T2 = typename std::remove_cv<T1>::type;
    using T  = typename std::remove_cv<T2>::type;
    // Get the name of the class
    if constexpr ( std::is_same<T, char>::value ) {
        copy( id.name, "char" );
    } else if constexpr ( std::is_same<T, unsigned char>::value ) {
        copy( id.name, "unsigned char" );
    } else if constexpr ( std::is_same<T, int8_t>::value ) {
        copy( id.name, "int8_t" );
    } else if constexpr ( std::is_same<T, uint8_t>::value ) {
        copy( id.name, "uint8_t" );
    } else if constexpr ( std::is_same<T, int>::value || std::is_same<T, int32_t>::value ) {
        copy( id.name, "int32_t" );
    } else if constexpr ( std::is_same<T, int64_t>::value ) {
        copy( id.name, "int64_t" );
    } else if constexpr ( std::is_same<T, unsigned>::value || std::is_same<T, uint32_t>::value ) {
        copy( id.name, "uint32_t" );
    } else if constexpr ( std::is_same<T, uint64_t>::value ) {
        copy( id.name, "unt64_t" );
    } else if constexpr ( std::is_same<T, float>::value ) {
        copy( id.name, "float" );
    } else if constexpr ( std::is_same<T, double>::value ) {
        copy( id.name, "double" );
    } else if constexpr ( std::is_same<T, std::complex<float>>::value ) {
        copy( id.name, "std::complex<float>" );
    } else if constexpr ( std::is_same<T, std::complex<double>>::value ) {
        copy( id.name, "std::complex<double>" );
    } else {
        // Get the name of the function to create the type name
        char name0[1024] = { 0 };
#ifdef __clang__
        copy( name0, __PRETTY_FUNCTION__ );
#elif defined( __GNUC__ )
        copy( name0, __PRETTY_FUNCTION__ );
#elif defined( _MSC_VER )
        copy( name0, __FUNCSIG__ );
#else
    #error "Not finished";
#endif
        // Get the type name from the function
        std::string_view name( name0 );
        if ( name.find( "]" ) != std::string::npos )
            name = name.substr( 0, name.find( "]" ) );
        if ( name.rfind( " = " ) != std::string::npos )
            name = name.substr( name.rfind( " = " ) + 3 );
        for ( size_t i = 0; i < name.size(); i++ )
            id.name[i] = name[i];
        id.name[name.size()] = 0;
    }
    // Create the hash
    if ( id.name[0] != 0 ) {
        id.hash = 5381;
        for ( unsigned char c : id.name ) {
            // hash = hash * 33 ^ c
            id.hash = ( ( id.hash << 5 ) + id.hash ) ^ c;
        }
    }
    // Set the size
    id.bytes = sizeof( T );
    return id;
}


} // namespace AMP

#endif
