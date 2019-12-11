#ifndef included_AMP_Database_hpp
#define included_AMP_Database_hpp

#include "AMP/utils/Database.h"

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
    virtual std::unique_ptr<KeyData> clone() const override
    {
        return std::make_unique<EmptyKeyData>();
    }
    virtual void print( std::ostream &os, const AMP::string_view & = "" ) const override
    {
        os << std::endl;
    }
    virtual AMP::string_view type() const override { return ""; }
    virtual bool is_floating_point() const override { return true; }
    virtual bool is_integral() const override { return true; }
    virtual std::vector<double> convertToDouble() const override { return std::vector<double>(); }
    virtual std::vector<int64_t> convertToInt64() const override { return std::vector<int64_t>(); }
    virtual bool operator==( const KeyData &rhs ) const override
    {
        return rhs.convertToDouble().empty();
    }
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
    virtual std::unique_ptr<KeyData> clone() const override
    {
        return std::make_unique<KeyDataScalar>( d_data, d_unit );
    }
    virtual void print( std::ostream &os, const AMP::string_view &indent = "" ) const override
    {
        os << indent;
        printValue( os, d_data );
        if ( !d_unit.isNull() )
            os << " " << d_unit.str();
        os << std::endl;
    }
    virtual AMP::string_view type() const override { return typeid( TYPE ).name(); }
    virtual bool is_floating_point() const override { return std::is_floating_point<TYPE>(); }
    virtual bool is_integral() const override { return std::is_integral<TYPE>(); }
    virtual std::vector<double> convertToDouble() const override;
    virtual std::vector<int64_t> convertToInt64() const override;
    virtual bool operator==( const KeyData &rhs ) const override;
    const TYPE &get() const { return d_data; }

private:
    TYPE d_data;
};
template<class TYPE>
class KeyDataVector final : public KeyData
{
public:
    KeyDataVector() = default;
    explicit KeyDataVector( std::vector<TYPE> data, const Units &unit = Units() )
        : KeyData( unit ), d_data( std::move( data ) )
    {
        data.clear(); // Suppress cppclean warning
    }
    virtual ~KeyDataVector() {}
    virtual std::unique_ptr<KeyData> clone() const override
    {
        return std::make_unique<KeyDataVector>( d_data, d_unit );
    }
    virtual void print( std::ostream &os, const AMP::string_view &indent = "" ) const override
    {
        os << indent;
        for ( size_t i = 0; i < d_data.size(); i++ ) {
            if ( i > 0 )
                os << ", ";
            printValue( os, d_data[i] );
        }
        if ( !d_unit.isNull() )
            os << " " << d_unit.str();
        os << std::endl;
    }
    virtual AMP::string_view type() const override { return typeid( TYPE ).name(); }
    virtual bool is_floating_point() const override { return std::is_floating_point<TYPE>(); }
    virtual bool is_integral() const override { return std::is_integral<TYPE>(); }
    virtual std::vector<double> convertToDouble() const override;
    virtual std::vector<int64_t> convertToInt64() const override;
    virtual bool operator==( const KeyData &rhs ) const override;
    const std::vector<TYPE> &get() const { return d_data; }

private:
    std::vector<TYPE> d_data;
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
    virtual std::unique_ptr<KeyData> clone() const override
    {
        return std::make_unique<DatabaseVector>( d_data );
    }
    virtual void print( std::ostream &os, const AMP::string_view &indent = "" ) const override
    {
        std::string indent2 = std::string( indent ) + "   ";
        for ( const auto &data : d_data ) {
            data.print( os, indent2 );
        }
    }
    virtual AMP::string_view type() const override
    {
        return typeid( std::vector<Database> ).name();
    }
    virtual bool is_floating_point() const override { return false; }
    virtual bool is_integral() const override { return false; }
    virtual std::vector<double> convertToDouble() const override
    {
        throw std::logic_error( "Cannot convert DatabaseVector to double" );
    }
    virtual std::vector<int64_t> convertToInt64() const override
    {
        throw std::logic_error( "Cannot convert DatabaseVector to int64" );
    }
    virtual bool operator==( const KeyData &rhs ) const override
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
void scaleData( std::vector<TYPE> &data, double factor );
template<class TYPE>
void scaleData( TYPE &data, double factor );
template<class TYPE>
std::vector<TYPE> convertFromDouble( const std::vector<double> &data );
template<class TYPE>
TYPE Database::getScalar( const AMP::string_view &key, Units unit ) const
{
    auto keyData = getData( key );
    DATABASE_INSIST( keyData, "Variable %s was not found in database", key.data() );
    TYPE data;
    auto scalarData = dynamic_cast<const KeyDataScalar<TYPE> *>( keyData );
    auto vectorData = dynamic_cast<const KeyDataVector<TYPE> *>( keyData );
    if ( scalarData ) {
        // We are dealing with a scalar of the same type
        data = scalarData->get();
    } else if ( vectorData ) {
        // We are dealing with a vector of the same type
        const auto &data2 = vectorData->get();
        DATABASE_INSIST( data2.size() == 1, "Variable %s is not a scalar", key.data() );
        data = data2[0];
    } else {
        // We need to convert the data
        auto data2 = convertFromDouble<TYPE>( keyData->convertToDouble() );
        DATABASE_INSIST( data2.size() == 1, "Variable %s is not a scalar", key.data() );
        data = data2[0];
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
TYPE Database::getWithDefault( const AMP::string_view &key, const TYPE &value, Units unit ) const
{
    auto keyData = getData( key );
    if ( !keyData )
        return value;
    TYPE data;
    auto scalarData = dynamic_cast<const KeyDataScalar<TYPE> *>( keyData );
    auto vectorData = dynamic_cast<const KeyDataVector<TYPE> *>( keyData );
    if ( scalarData ) {
        // We are dealing with a scalar of the same type
        data = scalarData->get();
    } else if ( vectorData ) {
        // We are dealing with a vector of the same type
        const auto &data2 = vectorData->get();
        DATABASE_INSIST( data2.size() == 1, "Variable %s is not a scalar", key.data() );
        data = data2[0];
    } else if constexpr ( std::is_same<TYPE, bool>::value ) {
        DATABASE_ERROR(
            "Unable to convert key %s to bool from %s", key.data(), keyData->type().data() );
    } else {
        // We need to convert the data
        auto data2 = convertFromDouble<TYPE>( keyData->convertToDouble() );
        DATABASE_INSIST( data2.size() == 1, "Variable %s is not a scalar", key.data() );
        data = data2[0];
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
std::vector<TYPE> Database::getVector( const AMP::string_view &key, Units unit ) const
{
    auto keyData = getData( key );
    DATABASE_INSIST( keyData, "Variable %s was not found in database", key.data() );
    std::vector<TYPE> data;
    auto scalarData = dynamic_cast<const KeyDataScalar<TYPE> *>( keyData );
    auto vectorData = dynamic_cast<const KeyDataVector<TYPE> *>( keyData );
    if ( scalarData ) {
        // We are dealing with a scalar of the same type
        const auto &data2 = scalarData->get();
        data.push_back( data2 );
    } else if ( vectorData ) {
        // We are dealing with a vector of the same type
        data = vectorData->get();
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
Database::putScalar<const char *>( const AMP::string_view &key, const char *value, Units unit )
{
    putScalar<std::string>( key, value, unit );
}
template<>
inline void Database::putScalar<AMP::string_view>( const AMP::string_view &key,
                                                   AMP::string_view value,
                                                   Units unit )
{
    putScalar<std::string>( key, std::string( value.data(), value.data() ), unit );
}
template<class TYPE>
inline void Database::putScalar( const AMP::string_view &key, TYPE value, Units unit )
{
    auto keyData = std::make_unique<KeyDataScalar<TYPE>>( std::move( value ), unit );
    putData( key, std::move( keyData ) );
}
template<class TYPE>
inline void Database::putVector( const AMP::string_view &key, std::vector<TYPE> data, Units unit )
{
    auto keyData = std::make_unique<KeyDataVector<TYPE>>( std::move( data ), unit );
    putData( key, std::move( keyData ) );
}


} // namespace AMP

#endif
