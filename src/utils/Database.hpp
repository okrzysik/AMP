#ifndef included_AMP_Database_hpp
#define included_AMP_Database_hpp

#include "AMP/utils/Database.h"
#include "AMP/utils/TypeTraits.h"

#include <iomanip>
#include <limits>
#include <tuple>


#define DATABASE_ERROR( ... )          \
    do {                               \
        char msg[1000];                \
        sprintf( msg, __VA_ARGS__ );   \
        throw std::logic_error( msg ); \
    } while ( 0 )
#define DATABASE_INSIST( TEST, ... )       \
    do {                                   \
        if ( !( TEST ) ) {                 \
            char msg[1000];                \
            sprintf( msg, __VA_ARGS__ );   \
            throw std::logic_error( msg ); \
        }                                  \
    } while ( 0 )


namespace AMP {


/********************************************************************
 * Create database from arguments                                    *
 ********************************************************************/
template<class TYPE, class... Args>
inline void Database::addArgs( const std::string_view &key, TYPE value, Args... args )
{
    if constexpr ( is_vector<TYPE>::value ) {
        putVector( key, value );
    } else if constexpr ( is_Array<TYPE>::value ) {
        putArray( key, value );
    } else if constexpr ( std::is_same<TYPE, std::string>::value ||
                          std::is_same<TYPE, std::string_view>::value ) {
        putScalar( key, value );
    } else if constexpr ( has_size<TYPE>::value ) {
        typedef decltype( *value.begin() ) TYPE2;
        typedef typename std::remove_reference<TYPE2>::type TYPE3;
        typedef typename std::remove_cv<TYPE3>::type TYPE4;
        std::vector<TYPE4> data( value.begin(), value.end() );
        putVector( key, std::move( data ) );
    } else {
        putScalar( key, value );
    }
    if constexpr ( sizeof...( args ) > 0 )
        addArgs( args... );
}
template<class TYPE, class... Args>
inline void Database::addArgsWithUnits( const std::string_view &key,
                                        TYPE value,
                                        const Units &unit,
                                        Args... args )
{
    if constexpr ( is_vector<TYPE>::value ) {
        putVector( key, value, unit );
    } else if constexpr ( is_Array<TYPE>::value ) {
        putArray( key, value );
    } else if constexpr ( std::is_same<TYPE, std::string>::value ||
                          std::is_same<TYPE, std::string_view>::value ) {
        putScalar( key, value, unit );
    } else if constexpr ( has_size<TYPE>::value ) {
        typedef decltype( *value.begin() ) TYPE2;
        typedef typename std::remove_reference<TYPE2>::type TYPE3;
        typedef typename std::remove_cv<TYPE3>::type TYPE4;
        std::vector<TYPE4> data( value.begin(), value.end() );
        putVector( key, std::move( data ), unit );
    } else {
        putScalar( key, value, unit );
    }
    if constexpr ( sizeof...( args ) > 0 )
        addArgsWithUnits( args... );
}
template<class... Args>
inline std::unique_ptr<Database> Database::create( Args... args )
{
    constexpr size_t n = sizeof...( args );
    static_assert( n % 2 == 0 );
    auto db = std::make_unique<Database>();
    if ( n == 0 )
        return db;
    db->addArgs( args... );
    return db;
}
template<class... Args>
inline std::unique_ptr<Database> Database::createWithUnits( Args... args )
{
    constexpr size_t n = sizeof...( args );
    static_assert( n % 3 == 0 );
    auto db = std::make_unique<Database>();
    if ( n == 0 )
        return db;
    db->addArgsWithUnits( args... );
    return db;
}


/********************************************************************
 * Basic classes for primative data types                            *
 ********************************************************************/
template<class TYPE>
inline void printValue( std::ostream &os, const TYPE &value )
{
    if constexpr ( std::is_floating_point<TYPE>::value ) {
        if ( value != value ) {
            os << "nan";
        } else if ( value == std::numeric_limits<TYPE>::infinity() ) {
            os << "inf";
        } else if ( value == -std::numeric_limits<TYPE>::infinity() ) {
            os << "-inf";
        } else {
            os << std::setprecision( 14 ) << value;
        }
    } else if constexpr ( std::is_same<TYPE, bool>::value ) {
        if ( value )
            os << "true";
        else
            os << "false";
    } else if constexpr ( std::is_same<TYPE, char>::value ||
                          std::is_same<TYPE, std::string>::value ) {
        os << '"' << value << '"';
    } else if constexpr ( std::is_integral<TYPE>::value ) {
        os << static_cast<int64_t>( value );
    } else {
        os << value;
    }
}
class EmptyKeyData final : public KeyData
{
public:
    EmptyKeyData() {}
    virtual ~EmptyKeyData() {}
    std::unique_ptr<KeyData> clone() const override { return std::make_unique<EmptyKeyData>(); }
    void print( std::ostream &os, const std::string_view & = "" ) const override
    {
        os << std::endl;
    }
    std::string_view type() const override { return ""; }
    bool is_floating_point() const override { return true; }
    bool is_integral() const override { return true; }
    ArraySize arraySize() const override { return ArraySize(); }
    Array<double> convertToDouble() const override { return Array<double>(); }
    Array<int64_t> convertToInt64() const override { return Array<int64_t>(); }
    bool operator==( const KeyData &rhs ) const override { return rhs.convertToDouble().empty(); }
};
template<class TYPE>
class KeyDataScalar final : public KeyData
{
public:
    KeyDataScalar() = default;
    explicit KeyDataScalar( TYPE data, const Units &unit = Units() )
        : KeyData( unit ), d_data( std::move( data ) )
    {
    }
    virtual ~KeyDataScalar() {}
    std::unique_ptr<KeyData> clone() const override
    {
        return std::make_unique<KeyDataScalar>( d_data, d_unit );
    }
    void print( std::ostream &os, const std::string_view &indent = "" ) const override
    {
        os << indent;
        printValue( os, d_data );
        if ( !d_unit.isNull() )
            os << " " << d_unit.str();
        os << std::endl;
    }
    std::string_view type() const override { return typeid( TYPE ).name(); }
    bool is_floating_point() const override { return std::is_floating_point<TYPE>(); }
    bool is_integral() const override { return std::is_integral<TYPE>(); }
    ArraySize arraySize() const override { return ArraySize( 1 ); }
    Array<double> convertToDouble() const override;
    Array<int64_t> convertToInt64() const override;
    bool operator==( const KeyData &rhs ) const override;
    const TYPE &get() const { return d_data; }

private:
    TYPE d_data;
};
template<class TYPE>
class KeyDataArray final : public KeyData
{
public:
    KeyDataArray() = default;
    explicit KeyDataArray( Array<TYPE> data, const Units &unit = Units() )
        : KeyData( unit ), d_data( std::move( data ) )
    {
        data.clear(); // Suppress cppclean warning
    }
    virtual ~KeyDataArray() {}
    std::unique_ptr<KeyData> clone() const override
    {
        return std::make_unique<KeyDataArray>( d_data, d_unit );
    }
    void print( std::ostream &os, const std::string_view &indent = "" ) const override
    {
        os << indent;
        if ( d_data.ndim() == 1 ) {
            printResursive( os, d_data );
        } else {
            os << "[";
            printResursive( os, d_data );
            os << "]";
        }
        if ( !d_unit.isNull() )
            os << " " << d_unit.str();
        os << std::endl;
    }
    static void printResursive( std::ostream &os, const Array<TYPE> &x )
    {
        int ndim = x.ndim();
        if ( ndim > 1 ) {
            for ( size_t i = 0; i < x.size( ndim - 1 ); i++ ) {
                if ( i > 0 )
                    os << ",";
                os << "[";
                size_t dims[5] = { 1 };
                for ( int d = 0; d < ndim - 1; d++ )
                    dims[d] = x.size( d );
                ArraySize s( ndim - 1, dims );
                Array<TYPE> y;
                y.viewRaw( s, const_cast<TYPE *>( &x( i * s.length() ) ) );
                printResursive( os, y );
                os << ']';
            }
        } else {
            for ( size_t i = 0; i < x.length(); i++ ) {
                if ( i > 0 )
                    os << ", ";
                printValue( os, x( i ) );
            }
        }
    }
    std::string_view type() const override { return typeid( TYPE ).name(); }
    bool is_floating_point() const override { return std::is_floating_point<TYPE>(); }
    bool is_integral() const override { return std::is_integral<TYPE>(); }
    ArraySize arraySize() const override { return d_data.size(); }
    Array<double> convertToDouble() const override;
    Array<int64_t> convertToInt64() const override;
    bool operator==( const KeyData &rhs ) const override;
    const Array<TYPE> &get() const { return d_data; }

private:
    Array<TYPE> d_data;
};
class DatabaseVector final : public KeyData
{
public:
    DatabaseVector() = default;
    explicit DatabaseVector( std::vector<Database> &&data ) : d_data( std::move( data ) ) {}
    explicit DatabaseVector( const std::vector<Database> &data )
    {
        d_data.resize( data.size() );
        for ( size_t i = 0; i < data.size(); i++ )
            d_data[i] = std::move( *data[i].cloneDatabase() );
    }
    DatabaseVector( DatabaseVector &&rhs ) : d_data( std::move( rhs.d_data ) ) {}
    DatabaseVector &operator=( DatabaseVector &&rhs )
    {
        d_data = std::move( rhs.d_data );
        return *this;
    }
    virtual ~DatabaseVector() {}
    std::unique_ptr<KeyData> clone() const override
    {
        return std::make_unique<DatabaseVector>( d_data );
    }
    void print( std::ostream &os, const std::string_view &indent = "" ) const override
    {
        std::string indent2 = std::string( indent ) + "   ";
        for ( const auto &data : d_data ) {
            data.print( os, indent2 );
        }
    }
    std::string_view type() const override { return typeid( std::vector<Database> ).name(); }
    bool is_floating_point() const override { return false; }
    bool is_integral() const override { return false; }
    ArraySize arraySize() const override { return ArraySize( d_data.size() ); }
    Array<double> convertToDouble() const override
    {
        throw std::logic_error( "Cannot convert DatabaseVector to double" );
    }
    Array<int64_t> convertToInt64() const override
    {
        throw std::logic_error( "Cannot convert DatabaseVector to int64" );
    }
    bool operator==( const KeyData &rhs ) const override
    {
        auto rhs2 = dynamic_cast<const DatabaseVector *>( &rhs );
        if ( rhs2 == nullptr )
            throw std::logic_error( "rhs is not a DatabaseVector" );
        return d_data == rhs2->d_data;
    }
    const std::vector<Database> &get() const { return d_data; }

private:
    std::vector<Database> d_data;
};


/********************************************************************
 * Get data                                                          *
 ********************************************************************/
template<class TYPE>
void scaleData( Array<TYPE> &data, double factor );
template<class TYPE>
void scaleData( TYPE &data, double factor );
template<class TYPE>
Array<TYPE> convertFromDouble( const Array<double> &data );
template<class TYPE>
TYPE Database::getScalar( const std::string_view &key, Units unit ) const
{
    auto keyData = getData( key );
    DATABASE_INSIST( keyData, "Variable %s was not found in database", key.data() );
    TYPE data;
    auto scalarData = dynamic_cast<const KeyDataScalar<TYPE> *>( keyData );
    auto arrayData  = dynamic_cast<const KeyDataArray<TYPE> *>( keyData );
    if ( scalarData ) {
        // We are dealing with a scalar of the same type
        data = scalarData->get();
    } else if ( arrayData ) {
        // We are dealing with an Array of the same type
        const auto &data2 = arrayData->get();
        DATABASE_INSIST( data2.length() == 1, "Variable %s is not a scalar", key.data() );
        data = data2( 0 );
    } else {
        // We need to convert the data
        auto data2 = convertFromDouble<TYPE>( keyData->convertToDouble() );
        DATABASE_INSIST( data2.length() == 1, "Variable %s is not a scalar", key.data() );
        data = data2( 0 );
    }
    if ( !unit.isNull() ) {
        auto unit2 = keyData->unit();
        DATABASE_INSIST( !unit2.isNull(), "Field %s must have units", key.data() );
        double factor = unit2.convert( unit );
        DATABASE_INSIST( factor != 0, "Unit conversion failed" );
        scaleData( data, factor );
    }
    return data;
}
template<class TYPE>
TYPE Database::getWithDefault( const std::string_view &key, const TYPE &value, Units unit ) const
{
    auto keyData = getData( key );
    if ( !keyData )
        return value;
    TYPE data;
    auto scalarData = dynamic_cast<const KeyDataScalar<TYPE> *>( keyData );
    auto arrayData  = dynamic_cast<const KeyDataArray<TYPE> *>( keyData );
    if ( scalarData ) {
        // We are dealing with a scalar of the same type
        data = scalarData->get();
    } else if ( arrayData ) {
        // We are dealing with an Array of the same type
        const auto &data2 = arrayData->get();
        DATABASE_INSIST( data2.size() == 1, "Variable %s is not a scalar", key.data() );
        data = data2( 0 );
    } else if constexpr ( std::is_same<TYPE, bool>::value ) {
        DATABASE_ERROR(
            "Unable to convert key %s to bool from %s", key.data(), keyData->type().data() );
    } else {
        // We need to convert the data
        auto data2 = convertFromDouble<TYPE>( keyData->convertToDouble() );
        DATABASE_INSIST( data2.size() == 1, "Variable %s is not a scalar", key.data() );
        data = data2( 0 );
    }
    if ( !unit.isNull() ) {
        auto unit2 = keyData->unit();
        DATABASE_INSIST( !unit2.isNull(), "Field %s must have units", key.data() );
        double factor = unit2.convert( unit );
        DATABASE_INSIST( factor != 0, "Unit conversion failed" );
        scaleData( data, factor );
    }
    return data;
}
template<class TYPE>
std::vector<TYPE> Database::getVector( const std::string_view &key, Units unit ) const
{
    auto keyData = getData( key );
    DATABASE_INSIST( keyData, "Variable %s was not found in database", key.data() );
    Array<TYPE> data;
    auto scalarData = dynamic_cast<const KeyDataScalar<TYPE> *>( keyData );
    auto arrayData  = dynamic_cast<const KeyDataArray<TYPE> *>( keyData );
    if ( scalarData ) {
        // We are dealing with a scalar of the same type
        const auto &data2 = scalarData->get();
        data.resize( 1 );
        data( 0 ) = data2;
    } else if ( arrayData ) {
        // We are dealing with an array of the same type
        data = arrayData->get();
    } else {
        // We need to convert the data
        data = std::move( convertFromDouble<TYPE>( keyData->convertToDouble() ) );
    }
    if ( !unit.isNull() && !data.empty() ) {
        auto unit2 = keyData->unit();
        DATABASE_INSIST( !unit2.isNull(), "Field %s must have units", key.data() );
        double factor = unit2.convert( unit );
        DATABASE_INSIST( factor != 0, "Unit conversion failed" );
        scaleData( data, factor );
    }
    DATABASE_INSIST( data.ndim() <= 1, "Variable %s cannot e converted to a vector", key.data() );
    std::vector<TYPE> data2( data.length() );
    for ( size_t i = 0; i < data.length(); i++ )
        data2[i] = data( i );
    return data2;
}
template<class TYPE>
Array<TYPE> Database::getArray( const std::string_view &key, Units unit ) const
{
    auto keyData = getData( key );
    DATABASE_INSIST( keyData, "Variable %s was not found in database", key.data() );
    Array<TYPE> data;
    auto scalarData = dynamic_cast<const KeyDataScalar<TYPE> *>( keyData );
    auto arrayData  = dynamic_cast<const KeyDataArray<TYPE> *>( keyData );
    if ( scalarData ) {
        // We are dealing with a scalar of the same type
        const auto &data2 = scalarData->get();
        data.resize( 1 );
        data( 0 ) = data2;
    } else if ( arrayData ) {
        // We are dealing with an array of the same type
        data = arrayData->get();
    } else {
        // We need to convert the data
        data = std::move( convertFromDouble<TYPE>( keyData->convertToDouble() ) );
    }
    if ( !unit.isNull() && !data.empty() ) {
        auto unit2 = keyData->unit();
        DATABASE_INSIST( !unit2.isNull(), "Field %s must have units", key.data() );
        double factor = unit2.convert( unit );
        DATABASE_INSIST( factor != 0, "Unit conversion failed" );
        scaleData( data, factor );
    }
    return data;
}


/********************************************************************
 * Put data                                                          *
 ********************************************************************/
template<>
inline void
Database::putScalar<const char *>( const std::string_view &key, const char *value, Units unit )
{
    putScalar<std::string>( key, value, unit );
}
template<>
inline void Database::putScalar<std::string_view>( const std::string_view &key,
                                                   std::string_view value,
                                                   Units unit )
{
    putScalar<std::string>( key, std::string( value.data(), value.data() ), unit );
}
template<class TYPE>
inline void Database::putScalar( const std::string_view &key, TYPE value, Units unit )
{
    auto keyData = std::make_unique<KeyDataScalar<TYPE>>( std::move( value ), unit );
    putData( key, std::move( keyData ) );
}
template<class TYPE>
inline void
Database::putVector( const std::string_view &key, const std::vector<TYPE> &data, Units unit )
{
    Array<TYPE> x;
    x            = data;
    auto keyData = std::make_unique<KeyDataArray<TYPE>>( std::move( x ), unit );
    putData( key, std::move( keyData ) );
}
template<class TYPE>
inline void Database::putArray( const std::string_view &key, Array<TYPE> data, Units unit )
{
    auto keyData = std::make_unique<KeyDataArray<TYPE>>( std::move( data ), unit );
    putData( key, std::move( keyData ) );
}


} // namespace AMP

#endif
