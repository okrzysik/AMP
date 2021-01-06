#include "AMP/utils/Database.h"
#include "AMP/utils/Array.h"
#include "AMP/utils/Utilities.h"

#include <algorithm>
#include <complex>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <string>
#include <tuple>


namespace AMP {


/********************************************************************
 * Helper functions                                                  *
 ********************************************************************/
static constexpr inline AMP::string_view deblank( const AMP::string_view &str )
{
    int i1 = 0, i2 = str.size() - 1;
    for ( ; i1 < (int) str.size() && ( str[i1] == ' ' || str[i1] == '\t' ); i1++ ) {}
    for ( ; i2 > 0 && ( str[i2] == ' ' || str[i2] == '\t' ); i2-- ) {}
    return i1 <= i2 ? str.substr( i1, i2 - i1 + 1 ) : AMP::string_view();
}
static inline bool strcmpi( const AMP::string_view &s1, const AMP::string_view &s2 )
{
    if ( s1.size() != s2.size() )
        return false;
    bool equal = true;
    for ( size_t i = 0; i < s1.size(); i++ ) {
        char a    = s1[i];
        char b    = s2[i];
        bool test = ( a & 0x1F ) == ( b & 0x1F ) && a > 64 && b > 64;
        equal     = equal && test;
    }
    return equal;
}
template<class TYPE>
static TYPE readValue( const AMP::string_view &str );
template<>
double readValue<double>( const AMP::string_view &str )
{
    double data = 0;
    if ( strcmpi( str, "inf" ) || strcmpi( str, "infinity" ) ) {
        data = std::numeric_limits<double>::infinity();
    } else if ( strcmpi( str, "inf" ) || strcmpi( str, "infinity" ) ) {
        data = -std::numeric_limits<double>::infinity();
    } else if ( strcmpi( str, "nan" ) ) {
        data = std::numeric_limits<double>::quiet_NaN();
    } else if ( str.find( '/' ) != std::string::npos ) {
        throw std::logic_error( "Error reading value" );
    } else {
        char *pos = nullptr;
        data      = strtod( str.data(), &pos );
        if ( static_cast<size_t>( pos - str.data() ) == str.size() + 1 )
            throw std::logic_error( "Error reading value" );
    }
    return data;
}
template<>
int readValue<int>( const AMP::string_view &str )
{
    char *pos = nullptr;
    int data  = strtol( str.data(), &pos, 10 );
    if ( static_cast<size_t>( pos - str.data() ) == str.size() + 1 )
        throw std::logic_error( "Error reading value" );
    return data;
}
template<>
std::complex<double> readValue<std::complex<double>>( const AMP::string_view &str )
{
    std::complex<double> data = 0;
    if ( str[0] != '(' ) {
        data = readValue<double>( str );
    } else {
        size_t pos = str.find( ',' );
        AMP::string_view s1( &str[1], pos - 1 );
        AMP::string_view s2( &str[pos + 1], str.size() - pos - 2 );
        data = std::complex<double>( readValue<double>( s1 ), readValue<double>( s2 ) );
    }
    return data;
}
template<class TYPE>
static std::tuple<TYPE, Units> readPair( const AMP::string_view &str )
{
    auto str0    = str;
    auto tmp     = deblank( std::move( str0 ) );
    size_t index = tmp.find( ' ' );
    if ( index != std::string::npos ) {
        return std::make_tuple( readValue<TYPE>( tmp.substr( 0, index ) ),
                                Units( tmp.substr( index + 1 ) ) );
    } else {
        return std::make_tuple( readValue<TYPE>( tmp ), Units() );
    }
}
static void strrep( std::string &str, const AMP::string_view &s, const AMP::string_view &r )
{
    size_t pos = str.find( s.data(), 0, s.size() );
    while ( pos != std::string::npos ) {
        str.replace( pos, s.size(), r.data(), r.size() );
        pos = str.find( s.data(), 0, s.size() );
    }
}


/********************************************************************
 * Constructors/destructor                                           *
 ********************************************************************/
Database::Database( Database &&rhs )
{
    std::swap( d_name, rhs.d_name );
    std::swap( d_hash, rhs.d_hash );
    std::swap( d_keys, rhs.d_keys );
    std::swap( d_data, rhs.d_data );
}
Database &Database::operator=( Database &&rhs )
{
    if ( this != &rhs ) {
        std::swap( d_name, rhs.d_name );
        std::swap( d_hash, rhs.d_hash );
        std::swap( d_keys, rhs.d_keys );
        std::swap( d_data, rhs.d_data );
    }
    return *this;
}


/********************************************************************
 * Clone the database                                                *
 ********************************************************************/
std::unique_ptr<KeyData> Database::clone() const
{
    auto db    = std::make_unique<Database>();
    db->d_name = d_name;
    db->d_hash = d_hash;
    db->d_keys = d_keys;
    db->d_data.resize( d_data.size() );
    for ( size_t i = 0; i < d_data.size(); i++ )
        db->d_data[i] = d_data[i]->clone();
    return db;
}
std::unique_ptr<Database> Database::cloneDatabase() const
{
    auto db    = std::make_unique<Database>();
    db->d_name = d_name;
    db->d_hash = d_hash;
    db->d_keys = d_keys;
    db->d_data.resize( d_data.size() );
    for ( size_t i = 0; i < d_data.size(); i++ )
        db->d_data[i] = d_data[i]->clone();
    return db;
}
void Database::copy( const Database &rhs )
{
    d_name = rhs.d_name;
    d_hash = rhs.d_hash;
    d_keys = rhs.d_keys;
    d_data.resize( rhs.d_data.size() );
    for ( size_t i = 0; i < d_data.size(); i++ )
        d_data[i] = rhs.d_data[i]->clone();
}


/********************************************************************
 * Check if the databases are equivalent                             *
 ********************************************************************/
bool Database::operator==( const Database &rhs ) const
{
    auto keys1 = getAllKeys();
    auto keys2 = rhs.getAllKeys();
    if ( keys1 != keys2 )
        return false;
    for ( const auto &key : keys1 ) {
        auto d1 = getData( key );
        auto d2 = rhs.getData( key );
        if ( *d1 != *d2 )
            return false;
    }
    return true;
}
bool Database::operator==( const KeyData &rhs ) const
{
    auto db = dynamic_cast<const Database *>( &rhs );
    if ( !db )
        return false;
    return operator==( *db );
}


/********************************************************************
 * Get the data object                                               *
 ********************************************************************/
bool Database::keyExists( const AMP::string_view &key ) const
{
    auto hash = hashString( key );
    int index = find( hash );
    return index != -1;
}
KeyData *Database::getData( const AMP::string_view &key )
{
    auto hash = hashString( key );
    int index = find( hash );
    return index == -1 ? nullptr : d_data[index].get();
}
const KeyData *Database::getData( const AMP::string_view &key ) const
{
    auto hash = hashString( key );
    int index = find( hash );
    return index == -1 ? nullptr : d_data[index].get();
}
bool Database::isDatabase( const AMP::string_view &key ) const
{
    auto hash = hashString( key );
    int index = find( hash );
    DATABASE_INSIST( index != -1, "Variable %s is not in database", key.data() );
    auto ptr2 = dynamic_cast<const Database *>( d_data[index].get() );
    return ptr2 != nullptr;
}
std::shared_ptr<Database> Database::getDatabase( const AMP::string_view &key )
{
    auto hash = hashString( key );
    int index = find( hash );
    DATABASE_INSIST( index != -1, "Variable %s is not in database", key.data() );
    auto ptr2 = std::dynamic_pointer_cast<Database>( d_data[index] );
    DATABASE_INSIST( ptr2, "Variable %s is not a database", key.data() );
    return ptr2;
}
std::shared_ptr<const Database> Database::getDatabase( const AMP::string_view &key ) const
{
    auto hash = hashString( key );
    int index = find( hash );
    DATABASE_INSIST( index != -1, "Variable %s is not in database", key.data() );
    auto ptr2 = std::dynamic_pointer_cast<const Database>( d_data[index] );
    DATABASE_INSIST( ptr2, "Variable %s is not a database", key.data() );
    return ptr2;
}
std::vector<std::string> Database::getAllKeys() const
{
    auto keys = d_keys;
    std::sort( keys.begin(), keys.end() );
    return keys;
}
void Database::putData( const AMP::string_view &key, std::unique_ptr<KeyData> data, bool check )
{
    auto hash = hashString( key );
    int index = find( hash );
    if ( index != -1 ) {
        if ( check )
            DATABASE_ERROR( "Variable %s already exists in database", key.data() );
        d_data[index] = std::move( data );
    } else {
        d_hash.emplace_back( hash );
        d_keys.emplace_back( key );
        d_data.emplace_back( std::move( data ) );
    }
}
void Database::erase( const AMP::string_view &key, bool check )
{
    auto hash = hashString( key );
    int index = find( hash );
    if ( index == -1 ) {
        if ( check )
            AMP_ERROR( std::string( key ) + " does not exist in database" );
        return;
    }
    std::swap( d_hash[index], d_hash.back() );
    std::swap( d_keys[index], d_keys.back() );
    std::swap( d_data[index], d_data.back() );
    d_hash.pop_back();
    d_keys.pop_back();
    d_data.pop_back();
}


/********************************************************************
 * Is the data of the given type                                     *
 ********************************************************************/
template<>
bool Database::isType<std::string>( const AMP::string_view &key ) const
{
    auto data = getData( key );
    DATABASE_INSIST( data, "Variable %s was not found in database", key.data() );
    auto type = data->type();
    return type == typeid( std::string ).name();
}
bool Database::isString( const AMP::string_view &key ) const { return isType<std::string>( key ); }
template<>
bool Database::isType<bool>( const AMP::string_view &key ) const
{
    auto data = getData( key );
    DATABASE_INSIST( data, "Variable %s was not found in database", key.data() );
    auto type  = data->type();
    auto type2 = typeid( bool ).name();
    return type == type2;
}
template<>
bool Database::isType<std::complex<float>>( const AMP::string_view &key ) const
{
    auto data = getData( key );
    DATABASE_INSIST( data, "Variable %s was not found in database", key.data() );
    auto type = data->type();
    return type == typeid( std::complex<float> ).name();
}
template<>
bool Database::isType<std::complex<double>>( const AMP::string_view &key ) const
{
    auto data = getData( key );
    DATABASE_INSIST( data, "Variable %s was not found in database", key.data() );
    auto type = data->type();
    return type == typeid( std::complex<double> ).name();
}
template<>
bool Database::isType<double>( const AMP::string_view &key ) const
{
    auto data = getData( key );
    DATABASE_INSIST( data, "Variable %s was not found in database", key.data() );
    auto type = data->type();
    if ( type == typeid( double ).name() )
        return true;
    bool is_floating = data->is_floating_point();
    bool is_integral = data->is_integral();
    return is_floating || is_integral;
}
template<>
bool Database::isType<DatabaseBox>( const AMP::string_view &key ) const
{
    auto data = getData( key );
    DATABASE_INSIST( data, "Variable %s was not found in database", key.data() );
    auto type = data->type();
    return type == typeid( DatabaseBox ).name();
}
template<class TYPE>
bool Database::isType( const AMP::string_view &key ) const
{
    auto data = getData( key );
    DATABASE_INSIST( data, "Variable %s was not found in database", key.data() );
    auto type = data->type();
    if ( type == typeid( TYPE ).name() )
        return true;
    if ( data->is_integral() ) {
        auto data2 = data->convertToInt64();
        bool pass  = true;
        for ( auto tmp : data2 )
            pass = pass && static_cast<int64_t>( static_cast<TYPE>( tmp ) ) == tmp;
        return pass;
    }
    if ( data->is_floating_point() ) {
        auto data2 = data->convertToDouble();
        bool pass  = true;
        for ( auto tmp : data2 )
            pass = pass && static_cast<double>( static_cast<TYPE>( tmp ) ) == tmp;
        return pass;
    }
    return false;
}
template bool Database::isType<char>( const AMP::string_view & ) const;
template bool Database::isType<uint8_t>( const AMP::string_view & ) const;
template bool Database::isType<uint16_t>( const AMP::string_view & ) const;
template bool Database::isType<uint32_t>( const AMP::string_view & ) const;
template bool Database::isType<uint64_t>( const AMP::string_view & ) const;
template bool Database::isType<int8_t>( const AMP::string_view & ) const;
template bool Database::isType<int16_t>( const AMP::string_view & ) const;
template bool Database::isType<int32_t>( const AMP::string_view & ) const;
template bool Database::isType<int64_t>( const AMP::string_view & ) const;
template bool Database::isType<float>( const AMP::string_view & ) const;
template bool Database::isType<long double>( const AMP::string_view & ) const;


/********************************************************************
 * Print the database                                                *
 ********************************************************************/
void Database::print( std::ostream &os, const AMP::string_view &indent ) const
{
    auto keys = getAllKeys(); //  We want the keys in sorted order
    for ( const auto &key : keys ) {
        os << indent << key;
        auto data  = getData( key );
        auto db    = dynamic_cast<const Database *>( data );
        auto dbVec = dynamic_cast<const DatabaseVector *>( data );
        if ( db ) {
            os << " {\n";
            db->print( os, std::string( indent ) + "   " );
            os << indent << "}\n";
        } else if ( dbVec ) {
            os << ":\n";
            dbVec->print( os, indent );
        } else {
            os << " = ";
            data->print( os, "" );
        }
    }
}
std::string Database::print( const AMP::string_view &indent ) const
{
    std::stringstream ss;
    print( ss, indent );
    return ss.str();
}


/********************************************************************
 * Read input database file                                          *
 ********************************************************************/
std::shared_ptr<Database> Database::parseInputFile( const std::string &filename )
{
    // Read the input file into memory
    FILE *fid = fopen( filename.data(), "rb" );
    DATABASE_INSIST( fid, "Error opening file %s", filename.data() );
    fseek( fid, 0, SEEK_END );
    size_t bytes = ftell( fid );
    rewind( fid );
    std::vector<char> buffer( bytes + 1 );
    size_t result = fread( buffer.data(), 1, bytes, fid );
    fclose( fid );
    DATABASE_INSIST( result == bytes, "Error reading file %s", filename.data() );
    buffer[bytes] = 0;
    // Create the database entries
    auto db = std::make_unique<Database>( filename );
    try {
        db->loadDatabase( buffer.data(), *db );
    } catch ( std::exception &err ) {
        throw std::logic_error( "Error loading database from file \"" + filename + "\"\n" +
                                err.what() );
    }
    return db;
}
std::unique_ptr<Database> Database::createFromString( const AMP::string_view &data )
{
    auto db = std::make_unique<Database>();
    loadDatabase( data.data(), *db );
    return db;
}
enum class token_type {
    newline,
    line_comment,
    block_start,
    block_stop,
    quote,
    comma,
    equal,
    bracket,
    end_bracket,
    end
};
static inline size_t length( token_type type )
{
    size_t len = 0;
    if ( type == token_type::newline || type == token_type::quote || type == token_type::equal ||
         type == token_type::bracket || type == token_type::end_bracket ||
         type == token_type::end ) {
        len = 1;
    } else if ( type == token_type::line_comment || type == token_type::block_start ||
                type == token_type::block_stop ) {
        len = 2;
    }
    return len;
}
static inline std::tuple<size_t, token_type> find_next_token( const char *buffer )
{
    size_t i = 0;
    while ( true ) {
        if ( buffer[i] == '\n' || buffer[i] == '\r' ) {
            return std::make_tuple( i + 1, token_type::newline );
        } else if ( buffer[i] == 0 ) {
            return std::make_tuple( i + 1, token_type::end );
        } else if ( buffer[i] == '"' ) {
            return std::make_tuple( i + 1, token_type::quote );
        } else if ( buffer[i] == ',' ) {
            return std::make_tuple( i + 1, token_type::comma );
        } else if ( buffer[i] == '=' ) {
            return std::make_tuple( i + 1, token_type::equal );
        } else if ( buffer[i] == '{' ) {
            return std::make_tuple( i + 1, token_type::bracket );
        } else if ( buffer[i] == '}' ) {
            return std::make_tuple( i + 1, token_type::end_bracket );
        } else if ( buffer[i] == '/' ) {
            if ( buffer[i + 1] == '/' ) {
                return std::make_tuple( i + 2, token_type::line_comment );
            } else if ( buffer[i + 1] == '*' ) {
                return std::make_tuple( i + 2, token_type::block_start );
            }
        } else if ( buffer[i] == '#' ) {
            return std::make_tuple( i + 1, token_type::line_comment );
        } else if ( buffer[i] == '*' ) {
            if ( buffer[i + 1] == '/' )
                return std::make_tuple( i + 2, token_type::block_stop );
        }
        i++;
    }
    return std::make_tuple<size_t, token_type>( 0, token_type::end );
}
static size_t skip_comment( const char *buffer )
{
    auto tmp          = find_next_token( buffer );
    auto comment_type = std::get<1>( tmp );
    auto end_comment =
        ( comment_type == token_type::line_comment ) ? token_type::newline : token_type::block_stop;
    size_t pos = 0;
    while ( std::get<1>( tmp ) != end_comment ) {
        if ( comment_type == token_type::block_start && std::get<1>( tmp ) == token_type::end )
            throw std::logic_error( "Encountered end of file before block comment end" );
        pos += std::get<0>( tmp );
        tmp = find_next_token( &buffer[pos] );
    }
    pos += std::get<0>( tmp );
    return pos;
}
enum class class_type { STRING, BOOL, INT, FLOAT, COMPLEX, BOX, UNKNOWN };
static std::tuple<size_t, std::unique_ptr<KeyData>> read_value( const char *buffer,
                                                                const AMP::string_view &key )
{
    // Split the value to an array of values
    size_t pos      = 0;
    token_type type = token_type::end;
    std::vector<AMP::string_view> values;
    class_type data_type = class_type::UNKNOWN;
    while ( type != token_type::newline ) {
        while ( buffer[pos] == ' ' || buffer[pos] == '\t' )
            pos++;
        size_t pos0 = pos;
        if ( buffer[pos0] == '(' ) {
            // We are dealing with a complex number
            data_type = class_type::COMPLEX;
            while ( buffer[pos] != ')' )
                pos++;
            size_t i;
            std::tie( i, type ) = find_next_token( &buffer[pos] );
            pos += i;
        } else if ( buffer[pos0] == '"' ) {
            // We are in a string
            data_type = class_type::STRING;
            pos++;
            while ( buffer[pos] != '"' )
                pos++;
            pos++;
            size_t i;
            std::tie( i, type ) = find_next_token( &buffer[pos] );
            pos += i;
        } else if ( buffer[pos0] == '[' && buffer[pos0 + 1] == '(' ) {
            // We are reading a SAMRAI box
            data_type = class_type::BOX;
            while ( buffer[pos] != ')' || buffer[pos + 1] != ']' )
                pos++;
            pos++;
            size_t i;
            std::tie( i, type ) = find_next_token( &buffer[pos] );
            pos += i;
        } else {
            std::tie( pos, type ) = find_next_token( &buffer[pos0] );
            pos += pos0;
        }
        AMP::string_view tmp( &buffer[pos0], pos - pos0 - length( type ) );
        if ( !tmp.empty() ) {
            if ( tmp.back() == ',' )
                tmp = AMP::string_view( tmp.data(), tmp.size() - 1 );
        }
        tmp = deblank( tmp );
        values.push_back( deblank( tmp ) );
        if ( type == token_type::comma ) {
            // We have multiple values
            continue;
        }
        if ( type == token_type::line_comment || type == token_type::block_start ) {
            // We encountered a comment
            pos += skip_comment( &buffer[pos - length( type )] ) - length( type );
            break;
        }
    }
    // Check if we are dealing with boolean values
    if ( strcmpi( values[0], "true" ) || strcmpi( values[0], "false" ) )
        data_type = class_type::BOOL;
    // Check if we are dealing with int
    if ( data_type == class_type::UNKNOWN ) {
        bool is_int = true;
        for ( size_t i = 0; i < values.size(); i++ ) {
            for ( size_t j = 0; j < values[i].size(); j++ ) {
                if ( values[i][j] < 42 || values[i][j] == 46 || values[i][j] >= 58 )
                    is_int = false;
            }
        }
        if ( is_int )
            data_type = class_type::INT;
    }
    // Default to an unknown type
    if ( data_type == class_type::UNKNOWN )
        data_type = class_type::FLOAT;
    // Convert the string value to the database value
    std::unique_ptr<KeyData> data;
    if ( values.empty() ) {
        data.reset( new EmptyKeyData() );
    } else if ( values.size() == 1 && values[0].empty() ) {
        data.reset( new EmptyKeyData() );
    } else if ( data_type == class_type::STRING ) {
        // We are dealing with strings
        for ( size_t i = 0; i < values.size(); i++ ) {
            if ( values[i][0] != '"' || values[i].back() != '"' )
                throw std::logic_error( "Error parsing string for key: " + std::string( key ) );
            values[i] = values[i].substr( 1, values[i].size() - 2 );
        }
        if ( values.size() == 1 ) {
            std::string str( values[0] );
            data = std::make_unique<KeyDataScalar<std::string>>( std::move( str ) );
        } else {
            std::vector<std::string> data2( values.size() );
            for ( size_t i = 0; i < values.size(); i++ )
                data2[i] = std::string( values[i].data(), values[i].size() );
            data = std::make_unique<KeyDataVector<std::string>>( std::move( data2 ) );
        }
    } else if ( data_type == class_type::BOOL ) {
        // We are dealing with logical values
        std::vector<bool> data2( values.size() );
        for ( size_t i = 0; i < values.size(); i++ ) {
            if ( !strcmpi( values[i], "true" ) && !strcmpi( values[i], "false" ) )
                throw std::logic_error( "Error converting " + std::string( key ) +
                                        " to logical array" );
            data2[i] = strcmpi( values[i], "true" );
        }
        if ( values.size() == 1 ) {
            data = std::make_unique<KeyDataScalar<bool>>( data2[0] );
        } else {
            data = std::make_unique<KeyDataVector<bool>>( std::move( data2 ) );
        }
    } else if ( data_type == class_type::INT ) {
        // We are dealing with floating point values
        std::vector<int> data2( values.size() );
        Units unit;
        for ( size_t i = 0; i < values.size(); i++ ) {
            Units unit2;
            std::tie( data2[i], unit2 ) = readPair<int>( values[i] );
            if ( !unit2.isNull() )
                unit = unit2;
        }
        if ( values.size() == 1 ) {
            data = std::make_unique<KeyDataScalar<int>>( data2[0], unit );
        } else {
            data = std::make_unique<KeyDataVector<int>>( std::move( data2 ), unit );
        }
    } else if ( data_type == class_type::FLOAT ) {
        // We are dealing with floating point values
        std::vector<double> data2( values.size() );
        Units unit;
        for ( size_t i = 0; i < values.size(); i++ ) {
            Units unit2;
            std::tie( data2[i], unit2 ) = readPair<double>( values[i] );
            if ( !unit2.isNull() )
                unit = unit2;
        }
        if ( values.size() == 1 ) {
            data = std::make_unique<KeyDataScalar<double>>( data2[0], unit );
        } else {
            data = std::make_unique<KeyDataVector<double>>( std::move( data2 ), unit );
        }
    } else if ( data_type == class_type::COMPLEX ) {
        // We are dealing with complex values
        std::vector<std::complex<double>> data2( values.size() );
        Units unit;
        for ( size_t i = 0; i < values.size(); i++ ) {
            Units unit2;
            std::tie( data2[i], unit2 ) = readPair<std::complex<double>>( values[i] );
            if ( !unit2.isNull() )
                unit = unit2;
        }
        if ( values.size() == 1 ) {
            data = std::make_unique<KeyDataScalar<std::complex<double>>>( data2[0], unit );
        } else {
            data =
                std::make_unique<KeyDataVector<std::complex<double>>>( std::move( data2 ), unit );
        }
    } else if ( data_type == class_type::BOX ) {
        std::vector<DatabaseBox> data2( values.size() );
        for ( size_t i = 0; i < values.size(); i++ )
            data2[i] = DatabaseBox( values[i] );
        if ( values.size() == 1 ) {
            data = std::make_unique<KeyDataScalar<DatabaseBox>>( data2[0] );
        } else {
            data = std::make_unique<KeyDataVector<DatabaseBox>>( std::move( data2 ) );
        }
    } else {
        // Treat unknown data as a string
        if ( values.size() == 1 ) {
            std::string str( values[0] );
            data = std::make_unique<KeyDataScalar<std::string>>( std::move( str ) );
        } else {
            std::vector<std::string> data2( values.size() );
            for ( size_t i = 0; i < values.size(); i++ )
                data2[i] = std::string( values[i].data(), values[i].size() );
            data = std::make_unique<KeyDataVector<std::string>>( std::move( data2 ) );
        }
    }
    return std::make_tuple( pos, std::move( data ) );
}
size_t Database::loadDatabase( const char *buffer, Database &db )
{
    size_t pos = 0;
    while ( true ) {
        size_t i;
        token_type type;
        std::tie( i, type ) = find_next_token( &buffer[pos] );
        AMP::string_view tmp( &buffer[pos], i - length( type ) );
        const auto key = deblank( tmp );
        if ( type == token_type::line_comment || type == token_type::block_start ) {
            // Comment
            DATABASE_INSIST( key.empty(), "Key should be empty: %s", key.data() );
            pos += skip_comment( &buffer[pos] );
        } else if ( type == token_type::newline ) {
            DATABASE_INSIST( key.empty(), "Key should be empty: %s", key.data() );
            pos += i;
        } else if ( type == token_type::equal ) {
            // Reading key/value pair
            DATABASE_INSIST( !key.empty(), "Empty key" );
            pos += i;
            std::unique_ptr<KeyData> data;
            std::tie( i, data ) = read_value( &buffer[pos], key );
            DATABASE_INSIST( data.get(), "null pointer" );
            db.putData( key, std::move( data ) );
            pos += i;
        } else if ( type == token_type::bracket ) {
            // Read database
            DATABASE_INSIST( !key.empty(), "Empty key" );
            pos += i;
            auto database = std::make_unique<Database>();
            pos += loadDatabase( &buffer[pos], *database );
            database->setName( std::string( key.data(), key.size() ) );
            db.putData( key, std::move( database ) );
        } else if ( type == token_type::end_bracket || type == token_type::end ) {
            // Finished with the database
            pos += i;
            break;
        } else {
            throw std::logic_error( "Error loading data" );
        }
    }
    return pos;
}
std::vector<double> Database::convertToDouble() const
{
    throw std::logic_error( "convertData on a database is not valid" );
}
std::vector<int64_t> Database::convertToInt64() const
{
    throw std::logic_error( "convertData on a database is not valid" );
}
bool Database::is_floating_point() const
{
    throw std::logic_error( "convertData on a database is not valid" );
}
bool Database::is_integral() const
{
    throw std::logic_error( "convertData on a database is not valid" );
}


/********************************************************************
 * Read YAML file                                                    *
 ********************************************************************/
static inline std::tuple<AMP::string_view, AMP::string_view>
splitYAML( const AMP::string_view &line )
{
    size_t pos = line.find_first_not_of( ' ' );
    if ( line[pos] == '-' )
        pos++;
    size_t pos2 = line.find( ':' );
    auto key    = deblank( line.substr( pos, pos2 - pos ) );
    auto value  = deblank( line.substr( pos2 + 1 ) );
    return std::tie( key, value );
}
static inline std::unique_ptr<KeyData> makeKeyData( std::vector<Database> &&data )
{
    if ( data.size() == 1 )
        return std::make_unique<Database>( std::move( data[0] ) );
    return std::make_unique<DatabaseVector>( std::move( data ) );
}
static inline void removeNewline( AMP::string_view &line )
{
    if ( line.find( '\n' ) != std::string::npos )
        line = line.substr( 0, line.find( '\n' ) );
    if ( line.find( '\r' ) != std::string::npos )
        line = line.substr( 0, line.find( '\r' ) );
}
static std::unique_ptr<KeyData> loadYAMLDatabase( FILE *fid, size_t indent = 0 )
{
    char tmp[4096];
    std::vector<Database> data;
    size_t indent2 = 10000;
    std::string lastKey;
    while ( true ) {
        tmp[0] = 0;
        if ( !fgets( tmp, sizeof( tmp ), fid ) )
            break;
        AMP::string_view line( tmp );
        // Remove newline
        removeNewline( line );
        // Remove the comments
        line = line.substr( 0, line.find( '#' ) );
        if ( line.empty() )
            continue;
        // Find the first non-whitespace character
        size_t pos = line.find_first_not_of( ' ' );
        if ( pos < indent ) {
            fseek( fid, -strlen( tmp ), SEEK_CUR );
            return makeKeyData( std::move( data ) );
        }
        if ( pos > indent2 ) {
            // We are dealing with a sub-item
            fseek( fid, -strlen( tmp ), SEEK_CUR );
            auto data2 = loadYAMLDatabase( fid, indent2 );
            data.back().putData( lastKey, std::move( data2 ) );
            continue;
        }
        if ( line[pos] == '-' ) {
            // We are dealing with a new item for the list
            data.resize( data.size() + 1 );
            indent2 = pos + 1 + line.substr( pos + 1 ).find_first_not_of( ' ' );
        }
        AMP::string_view key, value;
        std::tie( key, value ) = splitYAML( line );
        std::unique_ptr<KeyData> entry;
        if ( value == "|" ) {
            // Special case with block scalars
            Array<double> x;
            while ( true ) {
                char tmp2[4096] = { 0 };
                if ( !fgets( tmp2, sizeof( tmp2 ), fid ) )
                    break;
                AMP::string_view line2( tmp2 );
                removeNewline( line2 );
                size_t pos2 = line2.find_first_not_of( ' ' );
                if ( pos2 <= indent2 ) {
                    fseek( fid, -strlen( tmp2 ), SEEK_CUR );
                    break;
                }
                line2 = deblank( line2.substr( pos2 ) );
                std::string line3( line2 );
                strrep( line3, "  ", " " );
                strrep( line3, " ", "," );
                line3 += '\n';
                std::unique_ptr<KeyData> read_entry;
                std::tie( std::ignore, read_entry ) = read_value( line3.data(), key );
                auto y                              = read_entry->convertToDouble();
                size_t i                            = x.size( 0 );
                x.resize( i + 1, y.size() );
                for ( size_t j = 0; j < y.size(); j++ )
                    x( i, j ) = y[j];
            }
            // need to store a multi-dimentional array
            entry = std::make_unique<KeyDataScalar<std::string>>(
                std::string( "Not finished, need to store Array" ) );
        } else if ( !value.empty() ) {
            try {
                std::tie( std::ignore, entry ) = read_value( value.data(), key );
            } catch ( ... ) {
                entry = std::make_unique<KeyDataScalar<std::string>>( std::string( value ) );
            }
            try {
                if ( entry->convertToDouble() == std::vector<double>( 1, 0 ) )
                    entry = std::make_unique<KeyDataScalar<std::string>>( std::string( value ) );
            } catch ( ... ) {
            }
        }
        if ( entry ) {
            if ( data.empty() )
                data.resize( 1 );
            data.back().putData( key, std::move( entry ) );
        }
        lastKey = std::string( key.data(), key.size() );
    }
    return makeKeyData( std::move( data ) );
}
std::unique_ptr<KeyData> Database::readYAML( const AMP::string_view &filename )
{
    FILE *fid = fopen( filename.data(), "rb" );
    DATABASE_INSIST( fid, "Error opening file %s", filename.data() );
    auto data = loadYAMLDatabase( fid );
    fclose( fid );
    return data;
}


} // namespace AMP
