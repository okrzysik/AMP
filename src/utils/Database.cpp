#include "AMP/utils/Database.h"
#include "AMP/utils/Array.h"
#include "AMP/utils/Utilities.h"

#include <algorithm>
#include <complex>
#include <cstring>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <tuple>


namespace AMP {


// Forward declarations
static size_t
loadDatabase( const std::string &,
              const std::string_view &,
              size_t,
              Database &,
              std::map<std::string, const KeyData *> = std::map<std::string, const KeyData *>() );


/********************************************************************
 * Read the input file into memory                                   *
 ********************************************************************/
static std::vector<char> readFile( const std::string &filename )
{
    // Read the input file into memory
    FILE *fid = fopen( filename.data(), "rb" );
    DATABASE_INSIST( fid, "Error opening file %s", filename.data() );
    fseek( fid, 0, SEEK_END );
    size_t bytes = ftell( fid );
    rewind( fid );
    std::vector<char> data( bytes + 1, 0 );
    size_t result = fread( data.data(), 1, bytes, fid );
    fclose( fid );
    NULL_USE( result );
    return data;
}


/********************************************************************
 * Helper functions                                                  *
 ********************************************************************/
static constexpr inline std::string_view deblank( const std::string_view &str )
{
    if ( str.empty() )
        return std::string_view();
    int i1 = 0, i2 = str.size() - 1;
    for ( ; i1 < (int) str.size() && ( str[i1] == ' ' || str[i1] == '\t' ); i1++ ) {}
    for ( ; i2 > 0 && ( str[i2] == ' ' || str[i2] == '\t' || str[i2] == '\r' ); i2-- ) {}
    if ( i2 == 0 && ( str[i2] == ' ' || str[i2] == '\t' || str[i2] == '\r' ) )
        return std::string_view();
    return str.substr( i1, i2 - i1 + 1 );
}
static inline bool strcmpi( const std::string_view &s1, const std::string_view &s2 )
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
static TYPE readValue( const std::string_view &str );
template<>
double readValue<double>( const std::string_view &str )
{
    double data = 0;
    if ( strcmpi( str, "inf" ) || strcmpi( str, "infinity" ) ) {
        data = std::numeric_limits<double>::infinity();
    } else if ( strcmpi( str, "inf" ) || strcmpi( str, "infinity" ) ) {
        data = -std::numeric_limits<double>::infinity();
    } else if ( strcmpi( str, "nan" ) ) {
        data = std::numeric_limits<double>::quiet_NaN();
    } else if ( str.find( '/' ) != std::string::npos ) {
        AMP_ERROR( "Error reading value (double): " + std::string( str ) );
    } else {
        char *pos = nullptr;
        data      = strtod( str.data(), &pos );
        if ( static_cast<size_t>( pos - str.data() ) == str.size() + 1 )
            AMP_ERROR( "Error reading value (double): " + std::string( str ) );
    }
    return data;
}
template<>
int readValue<int>( const std::string_view &str )
{
    char *pos = nullptr;
    int data  = strtol( str.data(), &pos, 10 );
    if ( static_cast<size_t>( pos - str.data() ) == str.size() + 1 )
        AMP_ERROR( "Error reading value (int): " + std::string( str ) );
    return data;
}
template<class TYPE>
static std::tuple<TYPE, Units> readPair( const std::string_view &str )
{
    auto str0 = str;
    auto tmp  = deblank( std::move( str0 ) );
    if constexpr ( std::is_same<TYPE, std::complex<double>>::value ) {
        // We are trying to read a complex number
        if ( str[0] != '(' ) {
            // Read a double and convert to complex
            auto [value, unit] = readPair<double>( str );
            return std::make_tuple( std::complex<double>( value ), unit );
        }
        size_t p1 = str.find( ',' );
        size_t p2 = str.find( ')' );
        std::string_view s1( &str[1], p1 - 1 );
        std::string_view s2( &str[p1 + 1], p2 - p1 - 1 );
        auto value = std::complex<double>( readValue<double>( s1 ), readValue<double>( s2 ) );
        return std::make_tuple( value, Units( tmp.substr( p2 + 1 ) ) );
    } else {
        // Read a scalar type
        size_t index = tmp.find( ' ' );
        if ( index != std::string::npos ) {
            return std::make_tuple( readValue<TYPE>( tmp.substr( 0, index ) ),
                                    Units( tmp.substr( index + 1 ) ) );
        } else {
            return std::make_tuple( readValue<TYPE>( tmp ), Units() );
        }
    }
}
static void strrep( std::string &str, const std::string_view &s, const std::string_view &r )
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
Database::Database() : d_check( Check::WarnOverwrite ) {}
Database::Database( std::string name )
    : d_check( Check::WarnOverwrite ), d_name( std::move( name ) )
{
}
Database::Database( Database &&rhs )
{
    std::swap( d_check, rhs.d_check );
    std::swap( d_name, rhs.d_name );
    std::swap( d_hash, rhs.d_hash );
    std::swap( d_keys, rhs.d_keys );
    std::swap( d_data, rhs.d_data );
}
Database &Database::operator=( Database &&rhs )
{
    if ( this != &rhs ) {
        std::swap( d_check, rhs.d_check );
        std::swap( d_name, rhs.d_name );
        std::swap( d_hash, rhs.d_hash );
        std::swap( d_keys, rhs.d_keys );
        std::swap( d_data, rhs.d_data );
    }
    return *this;
}


/********************************************************************
 * Set default behavior for database                                 *
 ********************************************************************/
void Database::setDefaultAddKeyBehavior( Check check, bool setChildren )
{
    d_check = check;
    if ( setChildren ) {
        for ( auto tmp : d_data ) {
            auto db = std::dynamic_pointer_cast<Database>( tmp );
            if ( db )
                db->setDefaultAddKeyBehavior( check, true );
        }
    }
}


/********************************************************************
 * Clone the database                                                *
 ********************************************************************/
std::unique_ptr<KeyData> Database::clone() const
{
    auto db     = std::make_unique<Database>();
    db->d_check = d_check;
    db->d_name  = d_name;
    db->d_hash  = d_hash;
    db->d_keys  = d_keys;
    db->d_data.resize( d_data.size() );
    for ( size_t i = 0; i < d_data.size(); i++ )
        db->d_data[i] = d_data[i]->clone();
    return db;
}
std::unique_ptr<Database> Database::cloneDatabase() const
{
    auto db     = std::make_unique<Database>();
    db->d_check = d_check;
    db->d_name  = d_name;
    db->d_hash  = d_hash;
    db->d_keys  = d_keys;
    db->d_data.resize( d_data.size() );
    for ( size_t i = 0; i < d_data.size(); i++ )
        db->d_data[i] = d_data[i]->clone();
    return db;
}
void Database::copy( const Database &rhs )
{
    d_check = rhs.d_check;
    d_name  = rhs.d_name;
    d_hash  = rhs.d_hash;
    d_keys  = rhs.d_keys;
    d_data.resize( rhs.d_data.size() );
    for ( size_t i = 0; i < d_data.size(); i++ )
        d_data[i] = rhs.d_data[i]->clone();
}


/********************************************************************
 * Check if the databases are equivalent                             *
 ********************************************************************/
bool Database::operator==( const Database &rhs ) const
{
    auto keys1 = getAllKeys( true );
    auto keys2 = rhs.getAllKeys( true );
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
bool Database::keyExists( const std::string_view &key ) const
{
    auto hash = hashString( key );
    int index = find( hash );
    return index != -1;
}
KeyData *Database::getData( const std::string_view &key )
{
    auto hash = hashString( key );
    int index = find( hash );
    return index == -1 ? nullptr : d_data[index].get();
}
const KeyData *Database::getData( const std::string_view &key ) const
{
    auto hash = hashString( key );
    int index = find( hash );
    return index == -1 ? nullptr : d_data[index].get();
}
bool Database::isDatabase( const std::string_view &key ) const
{
    auto hash = hashString( key );
    int index = find( hash );
    DATABASE_INSIST( index != -1, "Variable %s is not in database", key.data() );
    auto ptr2 = dynamic_cast<const Database *>( d_data[index].get() );
    return ptr2 != nullptr;
}
std::shared_ptr<Database> Database::getDatabase( const std::string_view &key )
{
    auto hash = hashString( key );
    int index = find( hash );
    DATABASE_INSIST( index != -1, "Variable %s is not in database", key.data() );
    auto ptr2 = std::dynamic_pointer_cast<Database>( d_data[index] );
    DATABASE_INSIST( ptr2, "Variable %s is not a database", key.data() );
    return ptr2;
}
std::shared_ptr<const Database> Database::getDatabase( const std::string_view &key ) const
{
    auto hash = hashString( key );
    int index = find( hash );
    DATABASE_INSIST( index != -1, "Variable %s is not in database", key.data() );
    auto ptr2 = std::dynamic_pointer_cast<const Database>( d_data[index] );
    DATABASE_INSIST( ptr2, "Variable %s is not a database", key.data() );
    return ptr2;
}
const Database &Database::operator()( const std::string_view &key ) const
{
    auto hash = hashString( key );
    int index = find( hash );
    DATABASE_INSIST( index != -1, "Variable %s is not in database", key.data() );
    auto ptr2 = std::dynamic_pointer_cast<const Database>( d_data[index] );
    DATABASE_INSIST( ptr2, "Variable %s is not a database", key.data() );
    return *ptr2;
}
std::vector<std::string> Database::getAllKeys( bool sort ) const
{
    auto keys = d_keys;
    if ( sort )
        std::sort( keys.begin(), keys.end() );
    return keys;
}
void Database::putData( const std::string_view &key, std::unique_ptr<KeyData> data, Check check )
{
    if ( check == Check::GetDatabaseDefault )
        check = d_check;
    auto hash = hashString( key );
    int index = find( hash );
    if ( index != -1 ) {
        if ( check == Check::Error )
            DATABASE_ERROR( "Error: Variable '%s' already exists in database",
                            std::string( key ).data() );
        if ( check == Check::WarnOverwrite || check == Check::WarnKeep )
            DATABASE_WARNING( "Warning: variable '%s' already exists in database",
                              std::string( key ).data() );
        if ( check == Check::Overwrite || check == Check::WarnOverwrite )
            d_data[index] = std::move( data );
    } else {
        d_hash.emplace_back( hash );
        d_keys.emplace_back( key );
        d_data.emplace_back( std::move( data ) );
    }
}
void Database::erase( const std::string_view &key, bool check )
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
bool KeyData::isType<std::string>() const
{
    return type() == typeid( std::string ).name();
}
template<>
bool KeyData::isType<bool>() const
{
    return type() == typeid( bool ).name();
}
template<>
bool KeyData::isType<std::complex<float>>() const
{
    return type() == typeid( std::complex<float> ).name();
}
template<>
bool KeyData::isType<std::complex<double>>() const
{
    return type() == typeid( std::complex<double> ).name();
}
template<>
bool KeyData::isType<double>() const
{
    if ( type() == typeid( double ).name() )
        return true;
    return is_floating_point() || is_integral();
}
template<>
bool KeyData::isType<DatabaseBox>() const
{
    return type() == typeid( DatabaseBox ).name();
}
template<class TYPE>
bool KeyData::isType() const
{
    if ( type() == typeid( TYPE ).name() )
        return true;
    if ( is_integral() ) {
        auto data2 = convertToInt64();
        bool pass  = true;
        for ( auto tmp : data2 )
            pass = pass && static_cast<int64_t>( static_cast<TYPE>( tmp ) ) == tmp;
        return pass;
    }
    if ( is_floating_point() ) {
        auto data2 = convertToDouble();
        bool pass  = true;
        for ( auto tmp : data2 )
            pass = pass && static_cast<double>( static_cast<TYPE>( tmp ) ) == tmp;
        return pass;
    }
    return false;
}
template bool KeyData::isType<char>() const;
template bool KeyData::isType<uint8_t>() const;
template bool KeyData::isType<uint16_t>() const;
template bool KeyData::isType<uint32_t>() const;
template bool KeyData::isType<uint64_t>() const;
template bool KeyData::isType<int8_t>() const;
template bool KeyData::isType<int16_t>() const;
template bool KeyData::isType<int32_t>() const;
template bool KeyData::isType<int64_t>() const;
template bool KeyData::isType<float>() const;
template bool KeyData::isType<long double>() const;
bool Database::isString( const std::string_view &key ) const
{
    auto data = getData( key );
    DATABASE_INSIST( data, "Variable %s was not found in database", key.data() );
    return data->isType<std::string>();
}


/********************************************************************
 * Print the database                                                *
 ********************************************************************/
void Database::print( std::ostream &os, const std::string_view &indent, bool sort ) const
{
    auto keys = getAllKeys( sort ); //  We want the keys in sorted order
    for ( const auto &key : keys ) {
        os << indent << key;
        auto data  = getData( key );
        auto db    = dynamic_cast<const Database *>( data );
        auto dbVec = dynamic_cast<const DatabaseVector *>( data );
        if ( db ) {
            os << " {\n";
            db->print( os, std::string( indent ) + "   ", sort );
            os << indent << "}\n";
        } else if ( dbVec ) {
            os << ":\n";
            dbVec->print( os, indent, sort );
        } else {
            os << " = ";
            data->print( os, "", sort );
        }
    }
}
std::string Database::print( const std::string_view &indent, bool sort ) const
{
    std::stringstream ss;
    print( ss, indent, sort );
    return ss.str();
}


/********************************************************************
 * Read input database file                                          *
 ********************************************************************/
std::shared_ptr<Database> Database::parseInputFile( const std::string &filename )
{
    std::shared_ptr<Database> db;
    auto extension = filename.substr( filename.rfind( '.' ) + 1 );
    if ( extension == "yml" ) {
        std::shared_ptr<KeyData> data = readYAML( filename );
        auto db2                      = std::dynamic_pointer_cast<Database>( data );
        auto dbVec                    = std::dynamic_pointer_cast<DatabaseVector>( data );
        if ( db2 ) {
            db = db2;
        } else if ( dbVec ) {
            db = std::make_shared<Database>( filename );
            for ( const auto &db3 : dbVec->get() ) {
                auto key = db3.getName();
                AMP_ASSERT( !key.empty() );
                AMP_ASSERT( !db->keyExists( key ) );
                db->putDatabase( key, db3.cloneDatabase() );
            }
        } else {
            AMP_ERROR( "Unknown keyData" );
        }
        db->setName( filename );
    } else {
        db = std::make_shared<Database>( filename );
        db->readDatabase( filename );
    }
    return db;
}
void Database::readDatabase( const std::string &filename )
{
    // Read the input file into memory
    auto buffer = readFile( filename );
    // Create the database entries
    loadDatabase( "Error loading database from file \"" + filename + "\"\n",
                  buffer.data(),
                  buffer.size(),
                  *this );
}
std::unique_ptr<Database> Database::createFromString( const std::string_view &data )
{
    auto db = std::make_unique<Database>();
    loadDatabase( "Error creating database from file\n", data.data(), data.size(), *db );
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
    size_t pos        = 0;
    if ( comment_type == token_type::line_comment ) {
        // Line comment
        while ( std::get<1>( tmp ) != token_type::newline &&
                std::get<1>( tmp ) != token_type::end ) {
            pos += std::get<0>( tmp );
            tmp = find_next_token( &buffer[pos] );
        }
        pos += std::get<0>( tmp );
    } else {
        /* Block comment */
        while ( std::get<1>( tmp ) != token_type::block_stop ) {
            if ( comment_type == token_type::block_start && std::get<1>( tmp ) == token_type::end )
                throw std::logic_error( "Encountered end of file before block comment end" );
            pos += std::get<0>( tmp );
            tmp = find_next_token( &buffer[pos] );
        }
        pos += std::get<0>( tmp );
    }
    return pos;
}
enum class class_type { STRING, BOOL, INT, FLOAT, COMPLEX, BOX, ARRAY, DATABASE_ENTRY, UNKNOWN };
static class_type getType( std::string_view value0,
                           const std::map<std::string, const KeyData *> &databaseKeys )
{
    // Check for empty string
    if ( value0.empty() )
        return class_type::INT;
    // Check if we have units
    auto value = value0;
    if ( value.find( ' ' ) != std::string::npos ) {
        value = deblank( value.substr( 0, value.find( ' ' ) ) );
    }
    // Check if we are dealing with a simple bool
    if ( strcmpi( value, "true" ) || strcmpi( value, "false" ) )
        return class_type::BOOL;
    // Check if we could be an int or float
    bool is_int   = true;
    bool is_float = true;
    for ( size_t j = 0; j < value.size(); j++ ) {
        char c = value[j];
        if ( c < 42 || c == 46 || c >= 58 )
            is_int = false;
        if ( ( c < 42 || c >= 58 ) && ( c != 69 && c != 101 ) )
            is_float = false;
    }
    if ( is_int )
        return class_type::INT;
    if ( is_float )
        return class_type::FLOAT;
    // Check for special floating point types
    if ( strcmpi( value, "inf" ) || strcmpi( value, "nan" ) )
        return class_type::FLOAT;
    // Check if we are a database entry
    if ( databaseKeys.find( std::string( value0 ) ) != databaseKeys.end() )
        return class_type::DATABASE_ENTRY;
    // Return unknown
    return class_type::UNKNOWN;
}
template<class T>
static bool isType( const std::vector<const KeyData *> data )
{
    bool test = true;
    for ( auto tmp : data )
        test = test && tmp->isType<T>();
    return test;
}
// Convert the string value to the database value
static std::unique_ptr<KeyData>
createKeyData( const std::string_view &key,
               class_type data_type,
               std::vector<std::string_view> &values,
               const std::map<std::string, const KeyData *> &databaseKeys )
{
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
            Array<std::string> data2( values.size() );
            for ( size_t i = 0; i < values.size(); i++ )
                data2( i ) = std::string( values[i].data(), values[i].size() );
            data = std::make_unique<KeyDataArray<std::string>>( std::move( data2 ) );
        }
    } else if ( data_type == class_type::BOOL ) {
        // We are dealing with logical values
        Array<bool> data2( values.size() );
        for ( size_t i = 0; i < values.size(); i++ ) {
            if ( !strcmpi( values[i], "true" ) && !strcmpi( values[i], "false" ) )
                throw std::logic_error( "Error converting " + std::string( key ) +
                                        " to logical array" );
            data2( i ) = strcmpi( values[i], "true" );
        }
        if ( values.size() == 1 ) {
            data = std::make_unique<KeyDataScalar<bool>>( data2( 0 ) );
        } else {
            data = std::make_unique<KeyDataArray<bool>>( std::move( data2 ) );
        }
    } else if ( data_type == class_type::INT ) {
        // We are dealing with integer values
        Array<int> data2( values.size() );
        Units unit;
        for ( size_t i = 0; i < values.size(); i++ ) {
            Units unit2;
            std::tie( data2( i ), unit2 ) = readPair<int>( values[i] );
            if ( !unit2.isNull() )
                unit = unit2;
        }
        if ( values.size() == 1 ) {
            data = std::make_unique<KeyDataScalar<int>>( data2( 0 ), unit );
        } else {
            data = std::make_unique<KeyDataArray<int>>( std::move( data2 ), unit );
        }
    } else if ( data_type == class_type::FLOAT ) {
        // We are dealing with floating point values
        Array<double> data2( values.size() );
        Units unit;
        for ( size_t i = 0; i < values.size(); i++ ) {
            Units unit2;
            std::tie( data2( i ), unit2 ) = readPair<double>( values[i] );
            if ( !unit2.isNull() )
                unit = unit2;
        }
        if ( values.size() == 1 ) {
            data = std::make_unique<KeyDataScalar<double>>( data2( 0 ), unit );
        } else {
            data = std::make_unique<KeyDataArray<double>>( std::move( data2 ), unit );
        }
    } else if ( data_type == class_type::COMPLEX ) {
        // We are dealing with complex values
        Array<std::complex<double>> data2( values.size() );
        Units unit;
        for ( size_t i = 0; i < values.size(); i++ ) {
            Units unit2;
            std::tie( data2( i ), unit2 ) = readPair<std::complex<double>>( values[i] );
            if ( !unit2.isNull() )
                unit = unit2;
        }
        if ( values.size() == 1 ) {
            data = std::make_unique<KeyDataScalar<std::complex<double>>>( data2( 0 ), unit );
        } else {
            data = std::make_unique<KeyDataArray<std::complex<double>>>( std::move( data2 ), unit );
        }
    } else if ( data_type == class_type::BOX ) {
        Array<DatabaseBox> data2( values.size() );
        for ( size_t i = 0; i < values.size(); i++ )
            data2( i ) = DatabaseBox( values[i] );
        if ( values.size() == 1 ) {
            data = std::make_unique<KeyDataScalar<DatabaseBox>>( data2( 0 ) );
        } else {
            data = std::make_unique<KeyDataArray<DatabaseBox>>( std::move( data2 ) );
        }
    } else if ( data_type == class_type::ARRAY ) {
        // We are dealing with an Array
        size_t k = 0;
        for ( size_t i = 0; i < values[0].size(); i++ ) {
            if ( values[0][i] == ']' )
                k = i;
        }
        auto array_str = values[0].substr( 0, k + 1 );
        auto unit_str  = values[0].substr( k + 1 );
        Units unit( unit_str );
        if ( deblank( array_str.substr( 1, array_str.size() - 2 ) ).empty() ) {
            // We are dealing with an empty array
            data = std::make_unique<KeyDataArray<double>>( Array<double>( 0 ), unit );
        } else {
            // Get the array size
            size_t ndim = 1, dims[10] = { 1 };
            for ( size_t i = 1, d = 0; i < array_str.size() - 1; i++ ) {
                if ( array_str[i] == '[' ) {
                    d++;
                    ndim    = std::max( ndim, d + 1 );
                    dims[d] = 1;
                } else if ( array_str[i] == ']' ) {
                    d--;
                } else if ( array_str[i] == ',' ) {
                    dims[d]++;
                }
            }
            size_t dims2[10] = { 0 };
            for ( size_t d = 0; d < ndim; d++ )
                dims2[d] = dims[ndim - d - 1];
            ArraySize size( ndim, dims2 );
            // Get the Array values
            values.clear();
            values.resize( size.length() );
            for ( size_t i1 = 0, i2 = 0, j = 0; j < values.size(); j++ ) {
                while ( array_str[i1] == '[' || array_str[i1] == ']' || array_str[i1] == ',' ||
                        array_str[i1] == ' ' || array_str[i1] == '\n' )
                    i1++;
                i2 = i1 + 1;
                while ( array_str[i2] != '[' && array_str[i2] != ']' && array_str[i2] != ',' &&
                        array_str[i2] != ' ' )
                    i2++;
                values[j] = deblank( array_str.substr( i1, i2 - i1 ) );
                i1        = i2;
            }
            // Create the array
            if ( array_str.find( '"' ) != std::string::npos ) {
                // String array
                Array<std::string> A( size );
                for ( size_t i = 0; i < size.length(); i++ )
                    A( i ) = std::string( values[i].substr( 1, values[i].size() - 2 ) );
                data = std::make_unique<KeyDataArray<std::string>>( std::move( A ), unit );
            } else if ( array_str.find( '(' ) != std::string::npos ) {
                // Complex array
                Array<std::complex<double>> A( size );
                A.fill( 0.0 );
                for ( size_t i = 0; i < size.length(); i++ )
                    std::tie( A( i ), std::ignore ) = readPair<std::complex<double>>( values[i] );
                data = std::make_unique<KeyDataArray<std::complex<double>>>( std::move( A ), unit );
            } else if ( array_str.find( "true" ) != std::string::npos ||
                        array_str.find( "false" ) != std::string::npos ) {
                // Bool array
                Array<bool> A( size );
                A.fill( false );
                for ( size_t i = 0; i < values.size(); i++ ) {
                    if ( !strcmpi( values[i], "true" ) && !strcmpi( values[i], "false" ) )
                        throw std::logic_error( "Error converting " + std::string( key ) +
                                                " to logical array" );
                    A( i ) = strcmpi( values[i], "true" );
                }
                data = std::make_unique<KeyDataArray<bool>>( std::move( A ), unit );
            } else if ( array_str.find( '.' ) != std::string::npos ||
                        array_str.find( 'e' ) != std::string::npos ) {
                // Floating point array
                Array<double> A( size );
                A.fill( 0.0 );
                for ( size_t i = 0; i < size.length(); i++ )
                    std::tie( A( i ), std::ignore ) = readPair<double>( values[i] );
                data = std::make_unique<KeyDataArray<double>>( std::move( A ), unit );
            } else {
                // Integer point array
                Array<int> A( size );
                A.fill( 0 );
                for ( size_t i = 0; i < size.length(); i++ )
                    std::tie( A( i ), std::ignore ) = readPair<int>( values[i] );
                data = std::make_unique<KeyDataArray<int>>( std::move( A ), unit );
            }
        }
    } else if ( data_type == class_type::DATABASE_ENTRY ) {
        std::vector<const KeyData *> data2( values.size(), nullptr );
        for ( size_t i = 0; i < values.size(); i++ )
            data2[i] = databaseKeys.find( std::string( values[i] ) )->second;
        if ( values.size() == 1 ) {
            data = data2[0]->clone();
        } else {
            for ( size_t i = 0; i < values.size(); i++ ) {
                if ( data2[i]->arraySize().length() != 1 )
                    throw std::logic_error(
                        "Using multiple values from database only works for scalars: " +
                        std::string( key ) );
            }
            if ( isType<bool>( data2 ) ) {
                AMP::Array<bool> x( values.size() );
                for ( size_t i = 0; i < values.size(); i++ )
                    x( i ) = dynamic_cast<const KeyDataScalar<bool> *>( data2[i] )->get();
                data = std::make_unique<KeyDataArray<bool>>( std::move( x ) );
            } else if ( isType<int>( data2 ) ) {
                AMP::Array<int> x( values.size() );
                for ( size_t i = 0; i < values.size(); i++ )
                    x( i ) = dynamic_cast<const KeyDataScalar<int> *>( data2[i] )->get();
                data = std::make_unique<KeyDataArray<int>>( std::move( x ) );
            } else if ( isType<double>( data2 ) ) {
                AMP::Array<double> x( values.size() );
                for ( size_t i = 0; i < values.size(); i++ )
                    x( i ) = dynamic_cast<const KeyDataScalar<double> *>( data2[i] )->get();
                data = std::make_unique<KeyDataArray<double>>( std::move( x ) );
            } else if ( isType<std::string>( data2 ) ) {
                AMP::Array<std::string> x( values.size() );
                for ( size_t i = 0; i < values.size(); i++ )
                    x( i ) = dynamic_cast<const KeyDataScalar<std::string> *>( data2[i] )->get();
                data = std::make_unique<KeyDataArray<std::string>>( std::move( x ) );
            } else {
                throw std::logic_error(
                    "Using multiple values from database - unable to convert data: " +
                    std::string( key ) );
            }
        }
    } else if ( data_type == class_type::UNKNOWN ) {
        // Treat unknown data as a string
        if ( values.size() == 1 ) {
            std::string str( values[0] );
            data = std::make_unique<KeyDataScalar<std::string>>( std::move( str ) );
        } else {
            Array<std::string> data2( values.size() );
            for ( size_t i = 0; i < values.size(); i++ )
                data2( i ) = std::string( values[i].data(), values[i].size() );
            data = std::make_unique<KeyDataArray<std::string>>( std::move( data2 ) );
        }
    } else {
        throw std::logic_error( "Internal error" );
    }
    return data;
}
static std::tuple<size_t, std::unique_ptr<KeyData>>
read_value( const std::string_view &buffer,
            const std::string_view &key,
            const std::map<std::string, const KeyData *> &databaseKeys =
                std::map<std::string, const KeyData *>() )
{
    AMP_ASSERT( !buffer.empty() );
    // Split the value to an array of values
    size_t pos      = 0;
    token_type type = token_type::end;
    std::vector<std::string_view> values;
    class_type data_type = class_type::UNKNOWN;
    auto nextChar        = []( const char *c ) {
        while ( *c == ' ' )
            c++;
        return *c;
    };
    while ( type != token_type::newline ) {
        AMP_ASSERT( pos < buffer.size() );
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
        } else if ( buffer[pos0] == '[' && nextChar( &buffer[pos0 + 1] ) == '(' ) {
            // We are (probably) reading a SAMRAI box
            data_type = class_type::BOX;
            while ( buffer[pos] != ']' )
                pos++;
            size_t i;
            std::tie( i, type ) = find_next_token( &buffer[pos] );
            pos += i;
            // Check that it is a box and not an array of complex numbers
        } else if ( buffer[pos0] == '[' ) {
            // We are reading a multi-dimensional array
            data_type = class_type::ARRAY;
            int count = 1;
            pos       = pos0 + 1;
            while ( count != 0 && pos < buffer.size() ) {
                if ( buffer[pos] == '[' )
                    count++;
                if ( buffer[pos] == ']' )
                    count--;
                pos++;
            }
            size_t i;
            std::tie( i, type ) = find_next_token( &buffer[pos] );
            pos += i;
        } else {
            std::tie( pos, type ) = find_next_token( &buffer[pos0] );
            pos += pos0;
            if ( buffer[pos - 1] == '"' ) {
                while ( buffer[pos] != '"' )
                    pos++;
                size_t pos2           = pos + 1;
                std::tie( pos, type ) = find_next_token( &buffer[pos2] );
                pos += pos2;
            }
        }
        std::string_view tmp( &buffer[pos0], pos - pos0 - length( type ) );
        if ( !tmp.empty() ) {
            if ( tmp.back() == ',' )
                tmp = std::string_view( tmp.data(), tmp.size() - 1 );
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
    // Check the type of each entry
    if ( data_type == class_type::UNKNOWN ) {
        data_type = getType( values[0], databaseKeys );
        for ( size_t i = 1; i < values.size(); i++ ) {
            auto type = getType( values[i], databaseKeys );
            if ( type == class_type::UNKNOWN ) {
                data_type = class_type::UNKNOWN;
                break;
            } else if ( ( type == class_type::INT || type == class_type::FLOAT ) &&
                        ( data_type == class_type::INT || data_type == class_type::FLOAT ) ) {
                data_type = class_type::FLOAT;
            } else if ( type != data_type ) {
                throw std::logic_error( "Mismatched types in '" + std::string( key ) + "'" );
            }
        }
    }
    if ( data_type == class_type::UNKNOWN ) {
        throw std::logic_error( "Unkown type in '" + std::string( key ) + "':\n   " +
                                std::string( buffer.data(), pos ) );
    }
    // Convert the string value to the database value
    auto data = createKeyData( key, data_type, values, databaseKeys );
    return std::make_tuple( pos, std::move( data ) );
}
static size_t loadDatabase( const std::string &errMsgPrefix,
                            const std::string_view &buffer,
                            size_t N,
                            Database &db,
                            std::map<std::string, const KeyData *> databaseKeys )
{
    size_t pos = 0;
    while ( pos < N ) {
        size_t i;
        token_type type;
        std::tie( i, type ) = find_next_token( &buffer[pos] );
        std::string_view tmp( &buffer[pos], i - length( type ) );
        const auto key = deblank( tmp );
        if ( type == token_type::line_comment || type == token_type::block_start ) {
            // Comment
            AMP_INSIST( key.empty(), errMsgPrefix + "  Error reading comment" );
            pos += skip_comment( &buffer[pos] );
        } else if ( type == token_type::newline ) {
            if ( !key.empty() ) {
                auto msg = errMsgPrefix + "Key is not assigned a value: " + std::string( key );
                AMP_ERROR( msg );
            }
            pos += i;
        } else if ( type == token_type::equal ) {
            // Reading key/value pair
            DATABASE_INSIST( !key.empty(), "Empty key" );
            pos += i;
            std::unique_ptr<KeyData> data;
            try {
                std::tie( i, data ) = read_value( &buffer[pos], key, databaseKeys );
            } catch ( StackTrace::abort_error &err ) {
                AMP_ERROR( errMsgPrefix + "   Error loading key '" + std::string( key ) + "'\n" +
                           "      in file " + err.filename + " at line " +
                           std::to_string( err.line ) + "\n" + "   " + err.message );
            } catch ( std::exception &err ) {
                AMP_ERROR( errMsgPrefix + "    Error loading key '" + std::string( key ) +
                           "', unhandled exception:\n" + err.what() );
            }
            if ( !data )
                AMP_ERROR( errMsgPrefix + "    Error loading key '" + std::string( key ) + "'\n" +
                           " - empty data\n" );
            databaseKeys[std::string( key )] = data.get();
            db.putData( key, std::move( data ) );
            pos += i;
        } else if ( type == token_type::bracket ) {
            // Read database
            DATABASE_INSIST( !key.empty(), "Empty key" );
            pos += i;
            auto database = std::make_unique<Database>();
            pos += loadDatabase( errMsgPrefix, &buffer[pos], N - pos, *database, databaseKeys );
            database->setName( std::string( key ) );
            databaseKeys[std::string( key )] = database.get();
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
Array<double> Database::convertToDouble() const
{
    throw std::logic_error( "convertData on a database is not valid" );
}
Array<int64_t> Database::convertToInt64() const
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
static inline std::tuple<std::string_view, std::string_view>
splitYAML( const std::string_view &line )
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
static inline size_t getLine( const char *buffer, size_t pos )
{
    size_t i = pos;
    while ( buffer[i] != 0 && buffer[i] != '\n' ) {
        i++;
    }
    return i;
}
size_t loadYAMLDatabase( const char *buffer, Database &db, size_t pos = 0, size_t indent = 0 )
{
    std::string lastKey;
    while ( buffer[pos] != 0 ) {
        // Get the next line
        auto pos2 = getLine( buffer, pos );
        std::string_view line( &buffer[pos], pos2 - pos );
        // Remove the comments
        line = line.substr( 0, line.find( '#' ) );
        // Find the first non-whitespace character
        size_t p = line.find_first_not_of( ' ' );
        if ( p < indent )
            return pos; // End of list
        // Remove empty space
        line = deblank( line );
        if ( line.empty() ) {
            pos = pos2 + 1;
            continue;
        }
        if ( line[0] == '-' ) {
            // We are dealing with a new item (database)
            auto p2 = line.find( ':' );
            std::string name;
            if ( p2 != std::string::npos )
                name = deblank( line.substr( p2 + 1 ) );
            else
                name = deblank( line.substr( 1 ) );
            auto db2 = std::make_unique<Database>( name );
            pos      = loadYAMLDatabase( buffer, *db2, pos2 + 1, p + 1 );
            db.putDatabase( name, std::move( db2 ) );
            continue;
        }
        auto [key, value] = splitYAML( line );
        AMP_ASSERT( !key.empty() );
        if ( value.empty() ) {
            // Treat the key as a new database to load
            Database tmp;
            pos = loadYAMLDatabase( buffer, tmp, pos2 + 1, p + 1 );
            for ( auto key : tmp.getAllKeys() )
                db.putData( key, tmp.getData( key )->clone() );
            continue;
        } else if ( value == "|" ) {
            // Special case with block scalars
            pos2++;
            size_t pos3 = getLine( buffer, pos2 );
            std::string_view line2( &buffer[pos2 + 1], pos3 - pos2 );
            size_t p0 = line2.find_first_not_of( ' ' );
            Array<double> x;
            while ( true ) {
                size_t pos3 = getLine( buffer, pos2 );
                std::string_view line2( &buffer[pos2], pos3 - pos2 );
                size_t p2 = line2.find_first_not_of( ' ' );
                if ( p2 < p0 )
                    break;
                pos2  = pos3 + 1;
                line2 = deblank( line2 );
                std::string line3( line2 );
                strrep( line3, "  ", " " );
                strrep( line3, " ", "," );
                line3 += '\n';
                auto [pos4, read_entry] = read_value( line3.data(), key );
                NULL_USE( pos4 );
                auto y   = read_entry->convertToDouble();
                size_t i = x.size( 0 );
                x.resize( i + 1, y.length() );
                for ( size_t j = 0; j < y.length(); j++ )
                    x( i, j ) = y( j );
            }
            auto data = std::make_unique<KeyDataArray<double>>( std::move( x ) );
            db.putData( key, std::move( data ) );
        } else if ( !value.empty() ) {
            std::unique_ptr<KeyData> entry;
            try {
                std::tie( std::ignore, entry ) = read_value( value.data(), key );
            } catch ( ... ) {
                entry = std::make_unique<KeyDataScalar<std::string>>( std::string( value ) );
            }
            try {
                if ( entry->convertToDouble() == Array<double>( 1, 0 ) )
                    entry = std::make_unique<KeyDataScalar<std::string>>( std::string( value ) );
            } catch ( ... ) {
            }
            if ( !entry )
                AMP_ERROR( "Unable to parse value: " + std::string( value ) );
            db.putData( key, std::move( entry ) );
        } else {
            AMP_ERROR( "Not finished" );
        }
        pos = pos2 + 1;
    }
    return pos;
}
std::unique_ptr<KeyData> Database::readYAML( const std::string_view &filename )
{
    // Read the file into memory
    auto buffer = readFile( std::string( filename ) );
    // Read the file
    std::vector<Database> data;
    for ( size_t i = 0; i < buffer.size(); ) {
        data.resize( data.size() + 1 );
        i = loadYAMLDatabase( &buffer[i], data.back() ) + 1;
    }
    // Return the result
    return makeKeyData( std::move( data ) );
}


} // namespace AMP
