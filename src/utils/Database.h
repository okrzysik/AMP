#ifndef included_AMP_Database
#define included_AMP_Database

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "AMP/utils/Units.h"
#include "AMP/utils/string_view.h"
#include <memory>


// Forward declare SAMRAI's database
namespace SAMRAI {
namespace tbox {
class Database;
class DatabaseBox;
} // namespace tbox
} // namespace SAMRAI


// AMP namespace
namespace AMP {


//! Base class to hold data of a given type
class KeyData
{
public:
    //! Destructor
    virtual ~KeyData() {}
    //! Copy the data
    virtual std::unique_ptr<KeyData> clone() const = 0;
    //! Print the data to a stream
    virtual void print( std::ostream &os, const AMP::string_view &indent = "" ) const = 0;
    //! Return the native data type
    virtual AMP::string_view type() const = 0;
    //! Return true if the type is a floating point type
    virtual bool is_floating_point() const = 0;
    //! Return true if the type is a integer point type
    virtual bool is_integral() const = 0;
    //! Return the data as a std::vector<double> (throw error if this is not valid)
    virtual std::vector<double> convertToDouble() const = 0;
    //! Return the data as a std::vector<int64_t> (throw error if this is not valid)
    virtual std::vector<int64_t> convertToInt64() const = 0;
    //! Check if two sets of data are equal
    virtual bool operator==( const KeyData &rhs ) const = 0;
    //! Check if two sets of data are not equal
    inline bool operator!=( const KeyData &rhs ) const { return !operator==( rhs ); }
    //! Return the units
    const Units &unit() const { return d_unit; }

protected:
    KeyData() {}
    KeyData( const Units &unit ) : d_unit( unit ) {}
    KeyData( KeyData && )      = delete;
    KeyData( const KeyData & ) = delete;
    KeyData &operator=( KeyData && ) = delete;
    KeyData &operator=( const KeyData & ) = delete;

protected:
    Units d_unit;
};


//! Class to a database
class Database final : public KeyData
{
public:
    //! Empty constructor
    Database() = default;

    //! Basic constructor
    Database( std::string name ) : d_name( std::move( name ) ) {}

    /**
     * Open an database file.
     * @param filename       Name of input file to open
     */
    static std::shared_ptr<Database> parseInputFile( const std::string &filename );

    // Read a YAML database
    static std::unique_ptr<KeyData> readYAML( const AMP::string_view &filename );

    /** \brief Create a database from key/value pairs
     * \details  This function will create a database from a set of key/value pairs
     *    of the form: create( "key1", value1, "key2", value2, ... ).
     *    Note that for simplicity each value must either be a scalar value
     *       (int, double, string, etc) or a std::vector of scalar values
     * \param[in]  args         The input arguments
     */
    template<class... Args>
    static std::unique_ptr<Database> create( Args... args );

    /** \brief Create a database from key/value/unit triplets
     * \details  This function will create a database from a set of key/value/unit triplets
     *    of the form: create( "key1", value1, unit1, "key2", value2, unit2, ... ).
     *    Note that for simplicity each value must either be a scalar value
     *       (int, double, string, etc) or a std::vector of scalar values
     * \param[in]  args         The input arguments
     */
    template<class... Args>
    static std::unique_ptr<Database> createWithUnits( Args... args );

    /**
     * Create database from string
     * @param data       String containing the database data
     */
    static std::unique_ptr<Database> createFromString( const AMP::string_view &data );

    //! Copy constructor
    Database( const Database & ) = delete;

    //! Assignment operator
    Database &operator=( const Database & ) = delete;

    //! Move constructor
    Database( Database &&rhs );

    //! Move assignment operator
    Database &operator=( Database &&rhs );

    //! Destructor
    virtual ~Database() = default;

    //! Copy the data
    std::unique_ptr<KeyData> clone() const override;

    //! Copy the data
    void copy( const Database &rhs );

    //! Copy the data
    std::unique_ptr<Database> cloneDatabase() const;

    //! Get the name of the database
    inline const std::string &getName() const { return d_name; }

    //! Get the name of the database
    inline void setName( std::string name ) { d_name = std::move( name ); }

    /**
     * Return true if the specified key exists in the database and false
     *     otherwise.
     * @param[in] key           Key name to lookup.
     */
    bool keyExists( const AMP::string_view &key ) const;


    /**
     * Return all keys in the database.
     */
    std::vector<std::string> getAllKeys() const;


    //! Return true if the database is empty
    inline size_t empty() const { return d_data.size() == 0; }


    //! Return the number of entries in the database
    inline size_t size() const { return d_data.size(); }


    //! Return true if the databases are equivalent
    bool operator==( const Database &rhs ) const;


    //! Return true if the databases are equivalent
    bool operator==( const KeyData &rhs ) const override;


    //! Return the number of entries in the database
    inline bool operator!=( const Database &rhs ) const { return !operator==( rhs ); }


    /**
     * Get the key as a string
     *
     * @param[in] key           Key name in database.
     */
    inline std::string getString( const AMP::string_view &key ) const
    {
        return getScalar<std::string>( key );
    }


    /**
     * Get the scalar entry from the database with the specified key
     * name.  If the specified key does not exist in the database or
     * is not a scalar of the given type, then an error message is printed and
     * the program exits.
     *
     * @param[in] key           Key name in database.
     * @param[in] unit          Desired units
     */
    template<class TYPE>
    TYPE getScalar( const AMP::string_view &key, Units unit = Units() ) const;


    /**
     * Get the scalar entry from the database with the specified key
     * name.  If the specified key does not exist in the database the
     * the default value will be printed
     *
     * @param[in] key           Key name in database
     * @param[in] value         Default value
     * @param[in] unit          Desired units
     */
    template<class TYPE>
    TYPE
    getWithDefault( const AMP::string_view &key, const TYPE &value, Units unit = Units() ) const;


    /**
     * Get the vector entries from the database with the specified key
     * name.  If the specified key does not exist in the database or
     * is not of the given type, then an error message is printed and
     * the program exits.
     *
     * @param key           Key name in database.
     * @param unit          Desired units
     */
    template<class TYPE>
    std::vector<TYPE> getVector( const AMP::string_view &key, Units unit = Units() ) const;


    /**
     * Put the scalar entry into the database with the specified key name.
     * @param key           Key name in database.
     * @param value         Value to store
     * @param unit          Desired units
     */
    template<class TYPE>
    inline void putScalar( const AMP::string_view &key, TYPE value, Units unit = Units() );


    /**
     * Put the vector entries into the database with the specified key
     * name.  If the specified key does not exist in the database or
     * is not of the given type, then an error message is printed and
     * the program exits.
     *
     * @param key           Key name in database.
     * @param data          Data to store
     * @param unit          Desired units
     */
    template<class TYPE>
    inline void
    putVector( const AMP::string_view &key, std::vector<TYPE> data, Units unit = Units() );


    /**
     * Get a raw pointer to the data for a key in the database.
     * If the specified key does not exist, a null pointer is returned.
     *
     * @param key Key name in database.
     */
    KeyData *getData( const AMP::string_view &key );

    /**
     * Get a raw pointer to the data for a key in the database.
     * If the specified key does not exist, a null pointer is returned.
     *
     * @param key Key name in database.
     */
    const KeyData *getData( const AMP::string_view &key ) const;


    /**
     * Put the data for a key in the database.
     *
     * @param key       Key name in database.
     * @param data      Data to store
     * @param check     Check if the key exists and throw an error if does
     */
    void putData( const AMP::string_view &key, std::unique_ptr<KeyData> data, bool check = false );


    // Check if the key is a database object
    bool isDatabase( const AMP::string_view &key ) const;


    // Check if the key is a database object
    bool isString( const AMP::string_view &key ) const;


    // Check if the entry can be stored as the given type
    template<class TYPE>
    bool isType( const AMP::string_view &key ) const;


    /**
     * Get a raw pointer to the database for a key in the database.
     * If the specified key does not exist, a null pointer is returned.
     *
     * @param key Key name in database.
     */
    std::shared_ptr<Database> getDatabase( const AMP::string_view &key );

    /**
     * Get a raw pointer to the database for a key in the database.
     * If the specified key does not exist, a null pointer is returned.
     *
     * @param key Key name in database.
     */
    std::shared_ptr<const Database> getDatabase( const AMP::string_view &key ) const;


    /**
     * Put the database for a key in the database.
     * If the specified key already exists in the database an error is thrown.
     *
     * @param key       Key name in database.
     * @param db        Database to store
     */
    inline void putDatabase( const AMP::string_view &key, std::unique_ptr<Database> db )
    {
        putData( key, std::move( db ) );
    }


    /**
     * Create a database and return it
     * If the specified key already exists in the database an error is thrown.
     *
     * @param key       Key name in database.
     */
    inline std::shared_ptr<Database> putDatabase( const AMP::string_view &key )
    {
        putData( key, std::make_unique<Database>( std::string( key.data(), key.size() ) ) );
        return getDatabase( key );
    }


    /**
     * Erase the given key
     * If the specified key does not exists and check is true, an error is thrown.
     *
     * @param key       Key name in database
     * @param check     Check if the key exists
     */
    void erase( const AMP::string_view &key, bool check = true );


    /**
     * Print the data to a stream
     * @param os        Output stream
     * @param indent    Indenting to use before each line
     */
    void print( std::ostream &os, const AMP::string_view &indent = "" ) const override;


    //! Print the type
    AMP::string_view type() const override { return "database"; }


    /**
     * Print the data to a string
     * @return          Output string
     */
    std::string print( const AMP::string_view &indent = "" ) const;


#ifdef USE_SAMRAI
public: // SAMRAI interfaces
    //! Construct a database from a SAMRAI database
    Database( SAMRAI::tbox::Database & );

    //! Construct a database from a SAMRAI database
    inline Database( std::shared_ptr<SAMRAI::tbox::Database> ptr ) : Database( *ptr ) {}

    //! Create a SAMRAI database from this
    std::shared_ptr<SAMRAI::tbox::Database> cloneToSAMRAI() const;

#endif


protected: // Internal data and functions
    std::string d_name;
    std::vector<uint32_t> d_hash;
    std::vector<std::string> d_keys;
    std::vector<std::shared_ptr<KeyData>> d_data;

    // Function to load a database from a buffer
    static size_t loadDatabase( const char *buffer, Database &db );

    // Function to add arguments to the database
    template<class TYPE, class... Args>
    void addArgs( const AMP::string_view &key, TYPE value, Args... args );

    // Function to add arguments to the database
    template<class TYPE, class... Args>
    void
    addArgsWithUnits( const AMP::string_view &key, TYPE value, const Units &unit, Args... args );

    // Hash a string
    static constexpr uint32_t hashString( const AMP::string_view &s )
    {
        uint32_t hash = 5381;
        for ( size_t i = 0; i < s.size(); i++ )
            hash = ( ( hash << 5 ) + hash ) ^ s[i];
        return hash;
    }

    // Find an entry
    inline int find( uint32_t hash ) const
    {
        int index = -1;
        for ( size_t i = 0; i < d_hash.size(); i++ )
            if ( hash == d_hash[i] )
                index = i;
        return index;
    }

    // Functions inherited from KeyData that really aren't valid
    std::vector<double> convertToDouble() const override;
    std::vector<int64_t> convertToInt64() const override;
    bool is_floating_point() const override;
    bool is_integral() const override;
};


//! Class to store a box
class DatabaseBox final
{
public:
    //! Empty constructor
    DatabaseBox();

    //! Default constructor
    DatabaseBox( int dim, const int *lower, const int *upper );

    //! Construct from a string of the format "[(0,0,0), (7,7,7)]"
    DatabaseBox( const AMP::string_view &str );

#ifdef USE_SAMRAI
    //! Construct a DatabaseBox from a SAMRAI DatabaseBox
    DatabaseBox( const SAMRAI::tbox::DatabaseBox & );

    //! Create a SAMRAI DatabaseBox from this
    SAMRAI::tbox::DatabaseBox cloneToSAMRAI() const;
#endif

    //! Return a non-const reference to the number of dimensions
    uint8_t &dim();

    //! Return the number of dimensions
    uint8_t dim() const;

    //! Return whether the box is empty.
    //! A box is empty if it has dimension zero or if any part of the
    //! upper index is less than its corresponding part of the lower index.
    bool empty() const;

    //! Return non-const reference to the specified component of the lower index
    int &lower( uint8_t i );

    //! Return the specified component of the lower index
    int lower( uint8_t i ) const;

    //! Return non-const reference to the specified component of the upper index
    int &upper( uint8_t i );

    //! Return the specified component of the upper index
    int upper( uint8_t i ) const;

    bool operator==( const DatabaseBox &box ) const;

private:
    uint8_t d_dim;
    std::array<int, 5> d_lower, d_upper;
};

std::ostream &operator<<( std::ostream &out, const DatabaseBox & );


} // namespace AMP


#include "AMP/utils/Database.hpp"

#endif
