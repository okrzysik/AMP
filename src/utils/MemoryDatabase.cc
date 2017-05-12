//
// File:	$URL:
// file:///usr/casc/samrai/repository/AMP/tags/v-2-4-4/source/toolbox/database/MemoryDatabase.C $
// Package:	AMP toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2620 $
// Modified:	$LastChangedDate: 2008-11-19 14:24:28 -0800 (Wed, 19 Nov 2008) $
// Description:	An memory database structure that stores (key,value) pairs in memory
//
// NOTE: This is a modification of the MemoryDatabase class from SAMRAI
//       We have simply used it with modifications

#include "MemoryDatabase.h"

#include "Utilities.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdio.h>

#include <complex>
#include <stdlib.h>


#define PRINT_DEFAULT ( 1 )
#define PRINT_INPUT ( 2 )
#define PRINT_UNUSED ( 4 )

#define SSTREAM_BUFFER ( 4096 )

#define MEMORY_DB_ERROR( X )                                          \
    do {                                                              \
        pout << "MemoryDatabase: " << X << std::endl << std::flush;   \
        printClassData( pout );                                       \
        pout << "Program abort called..." << std::endl << std::flush; \
        comm.abort();                                                 \
    } while ( 0 )

namespace AMP {


/************************************************************************
*  Constructors/destructors                                             *
************************************************************************/
MemoryDatabase::MemoryDatabase( const std::string &name )
    : d_database_name( name ), comm( AMP_COMM_WORLD )
{
}


/************************************************************************
 *									                                    *
 * The virtual destructor deallocates database data.			        *
 *									                                    *
************************************************************************/

MemoryDatabase::~MemoryDatabase() {}


/************************************************************************
 *                                                                       *
 * Create memory data file specified by name.                            *
 *                                                                       *
************************************************************************/

bool MemoryDatabase::create( const std::string &name )
{
    d_database_name = name;
    d_keyvalues.clear();

    return true;
}


/************************************************************************
 *                                                                      *
 * Open memory data file specified by name                              *
 *                                                                      *
************************************************************************/

bool MemoryDatabase::open( const std::string &name )
{
    d_database_name = name;
    d_keyvalues.clear();

    return true;
}

/************************************************************************
 *                                                                       *
 * Close the open data file.                                             *
 *                                                                       *
************************************************************************/

bool MemoryDatabase::close()
{
    d_database_name = "";
    d_keyvalues.clear();

    return true;
}

/************************************************************************
 *									*
 * Return whether the key exists in the database.			*
 *									*
************************************************************************/

bool MemoryDatabase::keyExists( const std::string &key )
{
    return ( findKeyData( key ) ? true : false );
}

/************************************************************************
 *									*
 * Return all of the keys in the database.				*
 *									*
************************************************************************/

std::vector<std::string> MemoryDatabase::getAllKeys()
{
    const int n = d_keyvalues.size();
    std::vector<std::string> keys( n );

    int k = 0;
    for ( auto &elem : d_keyvalues ) {
        keys[k++] = ( elem ).d_key;
    }

    return ( keys );
}

/************************************************************************
 *									*
 * Get the type of the array entry associated with the specified key	*
 *									*
************************************************************************/
enum Database::DataType MemoryDatabase::getArrayType( const std::string &key )
{
    KeyData *keydata = findKeyData( key );

    if ( keydata ) {
        return keydata->d_type;
    } else {
        return Database::AMP_INVALID;
    }
}

/************************************************************************
 *									*
 * Get the size of the array entry associated with the specified key;	*
 * return 0 if the key does not exist.					*
 *									*
************************************************************************/

int MemoryDatabase::getArraySize( const std::string &key )
{
    KeyData *keydata = findKeyData( key );
    if ( keydata && keydata->d_type != Database::AMP_DATABASE ) {
        return keydata->d_array_size;
    } else {
        return 0;
    }
}

/************************************************************************
 *									*
 * Member functions that manage the database values within the database.	*
 *									*
************************************************************************/

bool MemoryDatabase::isDatabase( const std::string &key )
{
    KeyData *keydata = findKeyData( key );
    return ( keydata ? keydata->d_type == Database::AMP_DATABASE : false );
}

AMP::shared_ptr<Database> MemoryDatabase::putDatabase( const std::string &key )
{
    deleteKeyIfFound( key );
    KeyData keydata;
    keydata.d_key          = key;
    keydata.d_type         = Database::AMP_DATABASE;
    keydata.d_array_size   = 1;
    keydata.d_accessed     = false;
    keydata.d_from_default = false;
    keydata.d_database.reset( new MemoryDatabase( key ) );
    d_keyvalues.push_back( keydata );
    return ( keydata.d_database );
}

AMP::shared_ptr<Database> MemoryDatabase::getDatabase( const std::string &key )
{
    KeyData *keydata = findKeyDataOrExit( key );
    if ( keydata->d_type != Database::AMP_DATABASE ) {
        MEMORY_DB_ERROR( "Key=" << key << " is not a database..." );
    }
    keydata->d_accessed = true;
    return ( keydata->d_database );
}

/************************************************************************
 *									*
 * Member functions that manage boolean values within the database.	*
 *									*
************************************************************************/

bool MemoryDatabase::isBool( const std::string &key )
{
    KeyData *keydata = findKeyData( key );
    return ( keydata ? keydata->d_type == Database::AMP_BOOL : false );
}

void MemoryDatabase::putBool( const std::string &key, const bool &data )
{
    unsigned char uchar_data = (unsigned char) data;
    putBoolArray( key, &uchar_data, 1 );
}

void MemoryDatabase::putBoolArray( const std::string &key, const std::vector<unsigned char> &data )
{
    this->putBoolArray( key, &( data[0] ), static_cast<int>( data.size() ) );
}

void MemoryDatabase::putBoolArray( const std::string &key,
                                   const unsigned char *const data,
                                   const int nelements )
{
    deleteKeyIfFound( key );
    KeyData keydata;
    keydata.d_key          = key;
    keydata.d_type         = Database::AMP_BOOL;
    keydata.d_array_size   = nelements;
    keydata.d_accessed     = false;
    keydata.d_from_default = false;
    keydata.d_boolean      = std::vector<unsigned char>( nelements );

    for ( int i = 0; i < nelements; i++ ) {
        keydata.d_boolean[i] = data[i];
    }

    d_keyvalues.push_back( keydata );
}

bool MemoryDatabase::getBool( const std::string &key )
{
    KeyData *keydata = findKeyDataOrExit( key );
    if ( ( keydata->d_type != Database::AMP_BOOL ) || ( keydata->d_array_size != 1 ) ) {
        MEMORY_DB_ERROR( "Key=" << key << " is not a boolean scalar..." );
    }
    keydata->d_accessed = true;
    return ( keydata->d_boolean[0] );
}

bool MemoryDatabase::getBoolWithDefault( const std::string &key, const bool &defaultvalue )
{
    KeyData *keydata = findKeyData( key );
    if ( keydata )
        return ( this->getBool( key ) );
    putBool( key, defaultvalue );
    d_keyvalues.back().d_from_default = true;
    return ( defaultvalue );
}

std::vector<unsigned char> MemoryDatabase::getBoolArray( const std::string &key )
{
    KeyData *keydata = findKeyDataOrExit( key );
    if ( keydata->d_type != Database::AMP_BOOL ) {
        MEMORY_DB_ERROR( "Key=" << key << " is not a boolean..." );
    }
    keydata->d_accessed = true;
    return ( keydata->d_boolean );
}

void MemoryDatabase::getBoolArray( const std::string &key, bool *data, const int nelements )
{
    std::vector<unsigned char> tmp = this->getBoolArray( key );
    const int tsize                = tmp.size();

    if ( nelements != (int) tmp.size() ) {
        MEMORY_DB_ERROR( "Incorrect array size=" << nelements << " specified for key=" << key
                                                 << " with array size="
                                                 << tsize
                                                 << "..." );
    }

    for ( int i = 0; i < tsize; i++ ) {
        data[i] = tmp[i];
    }
}

/************************************************************************
 *									*
 * Member functions that manage box values within the database.		*
 *									*
************************************************************************/

bool MemoryDatabase::isDatabaseBox( const std::string &key )
{
    KeyData *keydata = findKeyData( key );
    return ( keydata ? keydata->d_type == Database::AMP_BOX : false );
}

void MemoryDatabase::putDatabaseBox( const std::string &key, const DatabaseBox &data )
{
    putDatabaseBoxArray( key, &data, 1 );
}

void MemoryDatabase::putDatabaseBoxArray( const std::string &key,
                                          const std::vector<DatabaseBox> &data )
{
    this->putDatabaseBoxArray( key, &( data[0] ), static_cast<int>( data.size() ) );
}

void MemoryDatabase::putDatabaseBoxArray( const std::string &key,
                                          const DatabaseBox *const data,
                                          const int nelements )
{
    deleteKeyIfFound( key );
    KeyData keydata;
    keydata.d_key          = key;
    keydata.d_type         = Database::AMP_BOX;
    keydata.d_array_size   = nelements;
    keydata.d_accessed     = false;
    keydata.d_from_default = false;
    keydata.d_box          = std::vector<DatabaseBox>( nelements );

    for ( int i = 0; i < nelements; i++ ) {
        keydata.d_box[i] = data[i];
    }

    d_keyvalues.push_back( keydata );
}

DatabaseBox MemoryDatabase::getDatabaseBox( const std::string &key )
{
    KeyData *keydata = findKeyDataOrExit( key );
    if ( ( keydata->d_type != Database::AMP_BOX ) || ( keydata->d_array_size != 1 ) ) {
        MEMORY_DB_ERROR( "Key=" << key << " is not a single box..." );
    }
    keydata->d_accessed = true;
    return ( keydata->d_box[0] );
}

DatabaseBox MemoryDatabase::getDatabaseBoxWithDefault( const std::string &key,
                                                       const DatabaseBox &defaultvalue )
{
    KeyData *keydata = findKeyData( key );
    if ( keydata )
        return ( this->getDatabaseBox( key ) );
    putDatabaseBox( key, defaultvalue );
    d_keyvalues.back().d_from_default = true;
    return ( defaultvalue );
}

std::vector<DatabaseBox> MemoryDatabase::getDatabaseBoxArray( const std::string &key )
{
    KeyData *keydata = findKeyDataOrExit( key );
    if ( keydata->d_type != Database::AMP_BOX ) {
        MEMORY_DB_ERROR( "Key=" << key << " is not a box..." );
    }
    keydata->d_accessed = true;
    return ( keydata->d_box );
}

void MemoryDatabase::getDatabaseBoxArray( const std::string &key,
                                          DatabaseBox *data,
                                          const int nelements )
{
    std::vector<DatabaseBox> tmp = this->getDatabaseBoxArray( key );
    const int tsize              = tmp.size();

    if ( nelements != (int) tmp.size() ) {
        MEMORY_DB_ERROR( "Incorrect array size=" << nelements << " specified for key=" << key
                                                 << " with array size="
                                                 << tsize
                                                 << "..." );
    }

    for ( int i = 0; i < tsize; i++ ) {
        data[i] = tmp[i];
    }
}

/************************************************************************
 *									*
 * Member functions that manage character values within the database.	*
 *									*
************************************************************************/

bool MemoryDatabase::isChar( const std::string &key )
{
    KeyData *keydata = findKeyData( key );
    return ( keydata ? keydata->d_type == Database::AMP_CHAR : false );
}

void MemoryDatabase::putChar( const std::string &key, const char &data )
{
    putCharArray( key, &data, 1 );
}

void MemoryDatabase::putCharArray( const std::string &key, const std::vector<char> &data )
{
    this->putCharArray( key, &( data[0] ), static_cast<int>( data.size() ) );
}

void MemoryDatabase::putCharArray( const std::string &key,
                                   const char *const data,
                                   const int nelements )
{
    deleteKeyIfFound( key );
    KeyData keydata;
    keydata.d_key          = key;
    keydata.d_type         = Database::AMP_CHAR;
    keydata.d_array_size   = nelements;
    keydata.d_accessed     = false;
    keydata.d_from_default = false;
    keydata.d_char         = std::vector<char>( nelements );

    for ( int i = 0; i < nelements; i++ ) {
        keydata.d_char[i] = data[i];
    }

    d_keyvalues.push_back( keydata );
}

char MemoryDatabase::getChar( const std::string &key )
{
    KeyData *keydata = findKeyDataOrExit( key );
    if ( ( keydata->d_type != Database::AMP_CHAR ) || ( keydata->d_array_size != 1 ) ) {
        MEMORY_DB_ERROR( "Key=" << key << " is not a single character..." );
    }
    keydata->d_accessed = true;
    return ( keydata->d_char[0] );
}

char MemoryDatabase::getCharWithDefault( const std::string &key, const char &defaultvalue )
{
    KeyData *keydata = findKeyData( key );
    if ( keydata )
        return ( this->getChar( key ) );
    putChar( key, defaultvalue );
    d_keyvalues.back().d_from_default = true;
    return ( defaultvalue );
}

std::vector<char> MemoryDatabase::getCharArray( const std::string &key )
{
    KeyData *keydata = findKeyDataOrExit( key );
    if ( keydata->d_type != Database::AMP_CHAR ) {
        MEMORY_DB_ERROR( "Key=" << key << " is not a character..." );
    }
    keydata->d_accessed = true;
    return ( keydata->d_char );
}

void MemoryDatabase::getCharArray( const std::string &key, char *data, const int nelements )
{
    std::vector<char> tmp = this->getCharArray( key );
    const int tsize       = tmp.size();

    if ( nelements != (int) tmp.size() ) {
        MEMORY_DB_ERROR( "Incorrect array size=" << nelements << " specified for key=" << key
                                                 << " with array size="
                                                 << tsize
                                                 << "..." );
    }

    for ( int i = 0; i < tsize; i++ ) {
        data[i] = tmp[i];
    }
}

/************************************************************************
 *									*
 * Member functions that manage complex values within the database.	*
 * Note that complex numbers may be promoted from integers, floats,	*
 * and doubles.								*
 *									*
************************************************************************/

bool MemoryDatabase::isComplex( const std::string &key )
{
    KeyData *keydata = findKeyData( key );
    return ( !keydata ? false : ( keydata->d_type == Database::AMP_COMPLEX ||
                                  keydata->d_type == Database::AMP_INT ||
                                  keydata->d_type == Database::AMP_FLOAT ||
                                  keydata->d_type == Database::AMP_DOUBLE ) );
}

void MemoryDatabase::putComplex( const std::string &key, const std::complex<double> &data )
{
    putComplexArray( key, &data, 1 );
}

void MemoryDatabase::putComplexArray( const std::string &key,
                                      const std::vector<std::complex<double>> &data )
{
    this->putComplexArray( key, &( data[0] ), static_cast<int>( data.size() ) );
}

void MemoryDatabase::putComplexArray( const std::string &key,
                                      const std::complex<double> *const data,
                                      const int nelements )
{
    deleteKeyIfFound( key );
    KeyData keydata;
    keydata.d_key          = key;
    keydata.d_type         = Database::AMP_COMPLEX;
    keydata.d_array_size   = nelements;
    keydata.d_accessed     = false;
    keydata.d_from_default = false;
    keydata.d_complex      = std::vector<std::complex<double>>( nelements );

    for ( int i = 0; i < nelements; i++ ) {
        keydata.d_complex[i] = data[i];
    }

    d_keyvalues.push_back( keydata );
}

std::complex<double> MemoryDatabase::getComplex( const std::string &key )
{
    std::complex<double> value( 0.0, 0.0 );
    KeyData *keydata = findKeyDataOrExit( key );

    if ( keydata->d_array_size != 1 ) {
        MEMORY_DB_ERROR( "Key=" << key << " is not a single complex..." );
    }

    switch ( keydata->d_type ) {
    case Database::AMP_INT:
        value = std::complex<double>( (double) keydata->d_integer[0], 0.0 );
        break;
    case Database::AMP_FLOAT:
        value = std::complex<double>( (double) keydata->d_float[0], 0.0 );
        break;
    case Database::AMP_DOUBLE:
        value = std::complex<double>( keydata->d_double[0], 0.0 );
        break;
    case Database::AMP_COMPLEX:
        value = keydata->d_complex[0];
        break;
    default:
        MEMORY_DB_ERROR( "Key=" << key << " is not a single complex..." );
    }

    keydata->d_accessed = true;
    return ( value );
}

std::complex<double>
MemoryDatabase::getComplexWithDefault( const std::string &key,
                                       const std::complex<double> &defaultvalue )
{
    KeyData *keydata = findKeyData( key );
    if ( keydata )
        return ( this->getComplex( key ) );
    putComplex( key, defaultvalue );
    d_keyvalues.back().d_from_default = true;
    return ( defaultvalue );
}

std::vector<std::complex<double>> MemoryDatabase::getComplexArray( const std::string &key )
{
    KeyData *keydata = findKeyDataOrExit( key );
    std::vector<std::complex<double>> array;
    switch ( keydata->d_type ) {
    case Database::AMP_INT: {
        array = std::vector<std::complex<double>>( keydata->d_integer.size() );
        for ( unsigned int i = 0; i < keydata->d_integer.size(); i++ ) {
            array[i] = std::complex<double>( (double) keydata->d_integer[i], 0.0 );
        }
        break;
    }
    case Database::AMP_FLOAT: {
        array = std::vector<std::complex<double>>( keydata->d_float.size() );
        for ( unsigned int i = 0; i < keydata->d_float.size(); i++ ) {
            array[i] = std::complex<double>( (double) keydata->d_float[i], 0.0 );
        }
        break;
    }
    case Database::AMP_DOUBLE: {
        array = std::vector<std::complex<double>>( keydata->d_double.size() );
        for ( unsigned int i = 0; i < keydata->d_float.size(); i++ ) {
            array[i] = std::complex<double>( keydata->d_double[i], 0.0 );
        }
        break;
    }
    case Database::AMP_COMPLEX:
        array = keydata->d_complex;
        break;
    default:
        MEMORY_DB_ERROR( "Key=" << key << " is not a complex..." );
    }
    keydata->d_accessed = true;
    return ( array );
}

void MemoryDatabase::getComplexArray( const std::string &key,
                                      std::complex<double> *data,
                                      const int nelements )
{
    std::vector<std::complex<double>> tmp = this->getComplexArray( key );
    const int tsize                       = tmp.size();

    if ( nelements != (int) tmp.size() ) {
        MEMORY_DB_ERROR( "Incorrect array size=" << nelements << " specified for key=" << key
                                                 << " with array size="
                                                 << tsize
                                                 << "..." );
    }

    for ( int i = 0; i < tsize; i++ ) {
        data[i] = tmp[i];
    }
}

/************************************************************************
 *									*
 * Member functions that manage double values within the database.	*
 * Note that doubles may be promoted from integers or floats.		*
 *									*
************************************************************************/

bool MemoryDatabase::isDouble( const std::string &key )
{
    KeyData *keydata = findKeyData( key );
    return ( !keydata ? false : ( keydata->d_type == Database::AMP_DOUBLE ||
                                  keydata->d_type == Database::AMP_INT ||
                                  keydata->d_type == Database::AMP_FLOAT ) );
}

void MemoryDatabase::putDouble( const std::string &key, const double &data )
{
    putDoubleArray( key, &data, 1 );
}

void MemoryDatabase::putDoubleArray( const std::string &key, const std::vector<double> &data )
{
    this->putDoubleArray( key, &( data[0] ), static_cast<int>( data.size() ) );
}

void MemoryDatabase::putDoubleArray( const std::string &key,
                                     const double *const data,
                                     const int nelements )
{
    deleteKeyIfFound( key );
    KeyData keydata;
    keydata.d_key          = key;
    keydata.d_type         = Database::AMP_DOUBLE;
    keydata.d_array_size   = nelements;
    keydata.d_accessed     = false;
    keydata.d_from_default = false;
    keydata.d_double       = std::vector<double>( nelements );

    for ( int i = 0; i < nelements; i++ ) {
        keydata.d_double[i] = data[i];
    }

    d_keyvalues.push_back( keydata );
}

double MemoryDatabase::getDouble( const std::string &key )
{
    double value     = 0.0;
    KeyData *keydata = findKeyDataOrExit( key );

    if ( keydata->d_array_size != 1 ) {
        MEMORY_DB_ERROR( "Key=" << key << " is not a single double..." );
    }

    switch ( keydata->d_type ) {
    case Database::AMP_INT:
        value = (double) keydata->d_integer[0];
        break;
    case Database::AMP_FLOAT:
        value = (double) keydata->d_float[0];
        break;
    case Database::AMP_DOUBLE:
        value = keydata->d_double[0];
        break;
    default:
        MEMORY_DB_ERROR( "Key=" << key << " is not a single double..." );
    }

    keydata->d_accessed = true;
    return ( value );
}

double MemoryDatabase::getDoubleWithDefault( const std::string &key, const double &defaultvalue )
{
    KeyData *keydata = findKeyData( key );
    if ( keydata )
        return ( this->getDouble( key ) );
    putDouble( key, defaultvalue );
    d_keyvalues.back().d_from_default = true;
    return ( defaultvalue );
}

std::vector<double> MemoryDatabase::getDoubleArray( const std::string &key )
{
    KeyData *keydata = findKeyDataOrExit( key );
    std::vector<double> array;
    switch ( keydata->d_type ) {
    case Database::AMP_INT: {
        array = std::vector<double>( keydata->d_integer.size() );
        for ( unsigned int i = 0; i < keydata->d_integer.size(); i++ ) {
            array[i] = (double) keydata->d_integer[i];
        }
        break;
    }
    case Database::AMP_FLOAT: {
        array = std::vector<double>( keydata->d_float.size() );
        for ( unsigned int i = 0; i < keydata->d_float.size(); i++ ) {
            array[i] = (double) keydata->d_float[i];
        }
        break;
    }
    case Database::AMP_DOUBLE: {
        array = keydata->d_double;
        break;
    }
    default:
        MEMORY_DB_ERROR( "Key=" << key << " is not a double..." );
    }
    keydata->d_accessed = true;
    return ( array );
}

void MemoryDatabase::getDoubleArray( const std::string &key, double *data, const int nelements )
{
    std::vector<double> tmp = this->getDoubleArray( key );
    const int tsize         = tmp.size();

    if ( nelements != (int) tmp.size() ) {
        MEMORY_DB_ERROR( "Incorrect array size=" << nelements << " specified for key=" << key
                                                 << " with array size="
                                                 << tsize
                                                 << "..." );
    }

    for ( int i = 0; i < tsize; i++ ) {
        data[i] = tmp[i];
    }
}

/************************************************************************
 *									*
 * Member functions that manage float values within the database.	*
 * Note that floats may be promoted from integers or truncated from	*
 * doubles (without a warning).						*
 *									*
************************************************************************/

bool MemoryDatabase::isFloat( const std::string &key )
{
    KeyData *keydata = findKeyData( key );
    return ( !keydata ? false : ( keydata->d_type == Database::AMP_DOUBLE ||
                                  keydata->d_type == Database::AMP_INT ||
                                  keydata->d_type == Database::AMP_FLOAT ) );
}

void MemoryDatabase::putFloat( const std::string &key, const float &data )
{
    putFloatArray( key, &data, 1 );
}

void MemoryDatabase::putFloatArray( const std::string &key, const std::vector<float> &data )
{
    this->putFloatArray( key, &( data[0] ), static_cast<int>( data.size() ) );
}

void MemoryDatabase::putFloatArray( const std::string &key,
                                    const float *const data,
                                    const int nelements )
{
    deleteKeyIfFound( key );
    KeyData keydata;
    keydata.d_key          = key;
    keydata.d_type         = Database::AMP_FLOAT;
    keydata.d_array_size   = nelements;
    keydata.d_accessed     = false;
    keydata.d_from_default = false;
    keydata.d_float        = std::vector<float>( nelements );

    for ( int i = 0; i < nelements; i++ ) {
        keydata.d_float[i] = data[i];
    }

    d_keyvalues.push_back( keydata );
}

float MemoryDatabase::getFloat( const std::string &key )
{

// Disable Intel warning about conversions
#ifdef __INTEL_COMPILER
#pragma warning( disable : 810 )
#endif

    float value      = 0.0;
    KeyData *keydata = findKeyDataOrExit( key );

    if ( keydata->d_array_size != 1 ) {
        MEMORY_DB_ERROR( "Key=" << key << " is not a single float..." );
    }

    switch ( keydata->d_type ) {
    case Database::AMP_INT:
        value = static_cast<float>( keydata->d_integer[0] );
        break;
    case Database::AMP_FLOAT:
        value = keydata->d_float[0];
        break;
    case Database::AMP_DOUBLE:
        value = static_cast<float>( keydata->d_double[0] );
        break;
    default:
        MEMORY_DB_ERROR( "Key=" << key << " is not a single float..." );
    }

    keydata->d_accessed = true;
    return ( value );
}

float MemoryDatabase::getFloatWithDefault( const std::string &key, const float &defaultvalue )
{
    KeyData *keydata = findKeyData( key );
    if ( keydata )
        return ( this->getFloat( key ) );
    putFloat( key, defaultvalue );
    d_keyvalues.back().d_from_default = true;
    return ( defaultvalue );
}

std::vector<float> MemoryDatabase::getFloatArray( const std::string &key )
{
    KeyData *keydata = findKeyDataOrExit( key );
    std::vector<float> array;
    switch ( keydata->d_type ) {
    case Database::AMP_INT: {
        array = std::vector<float>( keydata->d_integer.size() );
        for ( unsigned int i = 0; i < keydata->d_integer.size(); i++ ) {
            array[i] = static_cast<float>( keydata->d_integer[i] );
        }
        break;
    }
    case Database::AMP_FLOAT:
        array = keydata->d_float;
        break;
    case Database::AMP_DOUBLE: {
        array = std::vector<float>( keydata->d_double.size() );
        for ( unsigned int i = 0; i < keydata->d_double.size(); i++ ) {
            array[i] = static_cast<float>( keydata->d_double[i] );
        }
        break;
    }
    default:
        MEMORY_DB_ERROR( "Key=" << key << " is not a float..." );
    }
    keydata->d_accessed = true;
    return ( array );
}

void MemoryDatabase::getFloatArray( const std::string &key, float *data, const int nelements )
{
    std::vector<float> tmp = this->getFloatArray( key );
    const int tsize        = tmp.size();

    if ( nelements != (int) tmp.size() ) {
        MEMORY_DB_ERROR( "Incorrect array size=" << nelements << " specified for key=" << key
                                                 << " with array size="
                                                 << tsize
                                                 << "..." );
    }

    for ( int i = 0; i < tsize; i++ ) {
        data[i] = tmp[i];
    }
}

/************************************************************************
 *									*
 * Member functions that manage integer values within the database.	*
 *									*
************************************************************************/

bool MemoryDatabase::isInteger( const std::string &key )
{
    KeyData *keydata = findKeyData( key );
    return ( !keydata ? false : keydata->d_type == Database::AMP_INT );
}

void MemoryDatabase::putInteger( const std::string &key, const int &data )
{
    putIntegerArray( key, &data, 1 );
}

void MemoryDatabase::putIntegerArray( const std::string &key, const std::vector<int> &data )
{
    this->putIntegerArray( key, &( data[0] ), static_cast<int>( data.size() ) );
}

void MemoryDatabase::putIntegerArray( const std::string &key,
                                      const int *const data,
                                      const int nelements )
{
    deleteKeyIfFound( key );
    KeyData keydata;
    keydata.d_key          = key;
    keydata.d_type         = Database::AMP_INT;
    keydata.d_array_size   = nelements;
    keydata.d_accessed     = false;
    keydata.d_from_default = false;
    keydata.d_integer      = std::vector<int>( nelements );

    for ( int i = 0; i < nelements; i++ ) {
        keydata.d_integer[i] = data[i];
    }

    d_keyvalues.push_back( keydata );
}

int MemoryDatabase::getInteger( const std::string &key )
{
    KeyData *keydata = findKeyDataOrExit( key );
    if ( ( keydata->d_type != Database::AMP_INT ) || ( keydata->d_array_size != 1 ) ) {
        MEMORY_DB_ERROR( "Key=" << key << " is not an integer scalar..." );
    }
    keydata->d_accessed = true;
    return ( keydata->d_integer[0] );
}

int MemoryDatabase::getIntegerWithDefault( const std::string &key, const int &defaultvalue )
{
    KeyData *keydata = findKeyData( key );
    if ( keydata )
        return ( this->getInteger( key ) );
    putInteger( key, defaultvalue );
    d_keyvalues.back().d_from_default = true;
    return ( defaultvalue );
}

std::vector<int> MemoryDatabase::getIntegerArray( const std::string &key )
{
    KeyData *keydata = findKeyDataOrExit( key );
    if ( keydata->d_type != Database::AMP_INT ) {
        MEMORY_DB_ERROR( "Key=" << key << " is not an integer..." );
    }
    keydata->d_accessed = true;
    return ( keydata->d_integer );
}

void MemoryDatabase::getIntegerArray( const std::string &key, int *data, const int nelements )
{
    std::vector<int> tmp = this->getIntegerArray( key );
    const int tsize      = tmp.size();

    if ( nelements != (int) tmp.size() ) {
        MEMORY_DB_ERROR( "Incorrect array size=" << nelements << " specified for key=" << key
                                                 << " with array size="
                                                 << tsize
                                                 << "..." );
    }

    for ( int i = 0; i < tsize; i++ ) {
        data[i] = tmp[i];
    }
}

/************************************************************************
 *									*
 * Member functions that manage string values within the database.	*
 *									*
************************************************************************/

bool MemoryDatabase::isString( const std::string &key )
{
    KeyData *keydata = findKeyData( key );
    return ( !keydata ? false : keydata->d_type == Database::AMP_STRING );
}

void MemoryDatabase::putString( const std::string &key, const std::string &data )
{
    putStringArray( key, &data, 1 );
}

void MemoryDatabase::putStringArray( const std::string &key, const std::vector<std::string> &data )
{
    this->putStringArray( key, &( data[0] ), static_cast<int>( data.size() ) );
}

void MemoryDatabase::putStringArray( const std::string &key,
                                     const std::string *const data,
                                     const int nelements )
{
    deleteKeyIfFound( key );
    KeyData keydata;
    keydata.d_key          = key;
    keydata.d_type         = Database::AMP_STRING;
    keydata.d_array_size   = nelements;
    keydata.d_accessed     = false;
    keydata.d_from_default = false;
    keydata.d_string       = std::vector<std::string>( nelements );

    for ( int i = 0; i < nelements; i++ ) {
        keydata.d_string[i] = data[i];
    }

    d_keyvalues.push_back( keydata );
}

std::string MemoryDatabase::getString( const std::string &key )
{
    KeyData *keydata = findKeyDataOrExit( key );
    if ( ( keydata->d_type != Database::AMP_STRING ) || ( keydata->d_array_size != 1 ) ) {
        MEMORY_DB_ERROR( "Key=" << key << " is not a single string..." );
    }
    keydata->d_accessed = true;
    return ( keydata->d_string[0] );
}

std::string MemoryDatabase::getStringWithDefault( const std::string &key,
                                                  const std::string &defaultvalue )
{
    KeyData *keydata = findKeyData( key );
    if ( keydata )
        return ( this->getString( key ) );
    putString( key, defaultvalue );
    d_keyvalues.back().d_from_default = true;
    return ( defaultvalue );
}

std::vector<std::string> MemoryDatabase::getStringArray( const std::string &key )
{
    KeyData *keydata = findKeyDataOrExit( key );
    if ( keydata->d_type != Database::AMP_STRING ) {
        MEMORY_DB_ERROR( "Key=" << key << " is not a string..." );
    }
    keydata->d_accessed = true;
    return ( keydata->d_string );
}

void MemoryDatabase::getStringArray( const std::string &key,
                                     std::string *data,
                                     const int nelements )
{
    std::vector<std::string> tmp = this->getStringArray( key );
    const int tsize              = tmp.size();

    if ( nelements != (int) tmp.size() ) {
        MEMORY_DB_ERROR( "Incorrect array size=" << nelements << " specified for key=" << key
                                                 << " with array size="
                                                 << tsize
                                                 << "..." );
    }

    for ( int i = 0; i < tsize; i++ ) {
        data[i] = tmp[i];
    }
}

std::string MemoryDatabase::getName( void ) { return d_database_name; }

/************************************************************************
 *									*
 * Search the current database for a matching key.  If found, delete	*
 * that key and return true.  If the key does not exist, then return	*
 * false.								*
 *									*
************************************************************************/

bool MemoryDatabase::deleteKeyIfFound( const std::string &key )
{
    for ( std::list<KeyData>::iterator i = d_keyvalues.begin(); i != d_keyvalues.end(); i++ ) {
        if ( ( *i ).d_key == key ) {
            d_keyvalues.erase( i );
            return ( true );
        }
    }
    return ( false );
}

/************************************************************************
 *									*
 * Find the key data associated with the specified key and return a	*
 * pointer to the record.  If no such key data exists, then return NULL.	*
 *									*
************************************************************************/

MemoryDatabase::KeyData *MemoryDatabase::findKeyData( const std::string &key )
{
    for ( auto &elem : d_keyvalues ) {
        if ( key == ( elem ).d_key )
            return ( &( elem ) );
    }
    return ( nullptr );
}

/************************************************************************
 *									*
 * Find the key data associated with the specified key and return a	*
 * pointer to the record.  If no such key data exists, then exit with	*
 * an error message.							*
 *									*
************************************************************************/

MemoryDatabase::KeyData *MemoryDatabase::findKeyDataOrExit( const std::string &key )
{
    for ( auto &elem : d_keyvalues ) {
        if ( key == ( elem ).d_key )
            return ( &( elem ) );
    }
    MEMORY_DB_ERROR( "Key ``" << key << "'' does not exist in the database..." );
    return ( nullptr );
}

/************************************************************************
 *									*
 * Print the entire database to the specified output stream.	        *
 *									*
************************************************************************/

void MemoryDatabase::printClassData( std::ostream &os )
{
    printDatabase( os, 0, PRINT_DEFAULT | PRINT_INPUT | PRINT_UNUSED );
}

/************************************************************************
 *									*
 * Print unused database keys to the specified output stream.	        *
 *									*
************************************************************************/

void MemoryDatabase::printUnusedKeys( std::ostream &os ) const
{
    printDatabase( os, 0, PRINT_UNUSED );
}

/************************************************************************
 *									*
 * Print default database keys to the specified output stream.     	*
 *									*
************************************************************************/

void MemoryDatabase::printDefaultKeys( std::ostream &os ) const
{
    printDatabase( os, 0, PRINT_DEFAULT );
}

/************************************************************************
 *									*
 * Indent the output stream by the specified indentation factor.		*
 *									*
************************************************************************/

void MemoryDatabase::indentStream( std::ostream &os, const int indent )
{
    for ( int i = 0; i < indent; i++ ) {
        os << " ";
    }
}

/************************************************************************
 *									*
 * Print database data to the specified output stream.			*
 *									*
************************************************************************/

void MemoryDatabase::printDatabase( std::ostream &, const int, const int ) const
{

#if 0
    /*
     * Get the maximum key width in the output (excluding databases)
     */

    int width = 0;
    for (std::list<KeyData>::iterator k = d_keyvalues.begin(); k!=d_keyvalues.end(); k++) {
      if ( ( ((*k).d_from_default) && (toprint & PRINT_DEFAULT))
          || ( ((*k).d_accessed)     && (toprint & PRINT_INPUT  ))
          || (!((*k).d_accessed)     && (toprint & PRINT_UNUSED ))) {
        if ((*k).d_type != Database::AMP_DATABASE) {
          const int keywidth = (*k).d_key.length();
          if (keywidth > width) width = keywidth;
        }
      }
    }

    /*
     * Iterate over all non-database keys in the database and output key values
     */

    indentStream(os, indent);
    os << d_database_name << " {\n";
    for (std::list<KeyData>::iterator i = d_keyvalues.begin(); i!=d_keyvalues.end(); i++) {

      if ( ( ((*i).d_from_default) && (toprint & PRINT_DEFAULT))
          || ( ((*i).d_accessed)     && (toprint & PRINT_INPUT  ))
          || (!((*i).d_accessed)     && (toprint & PRINT_UNUSED ))) {

#ifndef LACKS_SSTREAM
        std::ostringstream sstream;
#else
        char sstream_buffer[SSTREAM_BUFFER];
        std::ostrstream sstream(sstream_buffer, SSTREAM_BUFFER);
#endif

        switch((*i).d_type) {

          case Database::AMP_INVALID: {
                                        break;
                                      }

          case Database::AMP_DATABASE: {
                                         break;
                                       }

          case Database::AMP_BOOL: {
                                     indentStream(sstream, indent+3);
                                     sstream << (*i).d_key;
                                     indentStream(sstream, width-(*i).d_key.length());
                                     sstream << " = ";
                                     const int n = (*i).d_boolean.size();
                                     for (int j = 0; j < n; j++) {
                                       sstream << ((*i).d_boolean[j] ? "TRUE" : "FALSE");
                                       if (j < n-1) sstream << ", ";
                                     }
                                     break;
                                   }

          case Database::AMP_BOX: {
                                    indentStream(sstream, indent+3);
                                    sstream << (*i).d_key;
                                    indentStream(sstream, width-(*i).d_key.length());
                                    sstream << " = ";
                                    const int n = (*i).d_box.size();
                                    for (int j = 0; j < n; j++) {
                                      const int m = (*i).d_box[j].getDimension();
                                      sstream << "[(";
                                      for (int k = 0; k < m; k++) {
                                        sstream << (*i).d_box[j].lower(k);
                                        if (k < m-1) sstream << ",";
                                      }
                                      sstream << "),(";
                                      for (int l = 0; l < m; l++) {
                                        sstream << (*i).d_box[j].upper(l);
                                        if (l < m-1) sstream << ",";
                                      }
                                      sstream << ")]";
                                      if (j < n-1) sstream << ", ";
                                    }
                                    break;
                                  }

          case Database::AMP_CHAR: {
                                     indentStream(sstream, indent+3);
                                     sstream << (*i).d_key;
                                     indentStream(sstream, width-(*i).d_key.length());
                                     sstream << " = ";
                                     const int n = (*i).d_char.size();
                                     for (int j = 0; j < n; j++) {
                                       sstream << "'" << (*i).d_char[j] << "'";
                                       if (j < n-1) sstream << ", ";
                                     }
                                     break;
                                   }

          case Database::AMP_COMPLEX: {
                                        indentStream(sstream, indent+3);
                                        sstream << (*i).d_key;
                                        indentStream(sstream, width-(*i).d_key.length());
                                        sstream << " = ";
                                        const int n = (*i).d_complex.size();
                                        for (int j = 0; j < n; j++) {
                                          sstream << (*i).d_complex[j];
                                          if (j < n-1) sstream << ", ";
                                        }
                                        break;
                                      }

          case Database::AMP_DOUBLE: {
                                       indentStream(sstream, indent+3);
                                       sstream << (*i).d_key;
                                       indentStream(sstream, width-(*i).d_key.length());
                                       sstream << " = ";
                                       const int n = (*i).d_double.size();
                                       for (int j = 0; j < n; j++) {
                                         sstream << (*i).d_double[j];
                                         if (j < n-1) sstream << ", ";
                                       }
                                       break;
                                     }

          case Database::AMP_FLOAT: {
                                      indentStream(sstream, indent+3);
                                      sstream << (*i).d_key;
                                      indentStream(sstream, width-(*i).d_key.length());
                                      sstream << " = ";
                                      const int n = (*i).d_float.size();
                                      for (int j = 0; j < n; j++) {
                                        sstream << (*i).d_float[j];
                                        if (j < n-1) sstream << ", ";
                                      }
                                      break;
                                    }

          case Database::AMP_INT: {
                                    indentStream(sstream, indent+3);
                                    sstream << (*i).d_key;
                                    indentStream(sstream, width-(*i).d_key.length());
                                    sstream << " = ";
                                    const int n = (*i).d_integer.size();
                                    for (int j = 0; j < n; j++) {
                                      sstream << (*i).d_integer[j];
                                      if (j < n-1) sstream << ", ";
                                    }
                                    break;
                                  }

          case Database::AMP_STRING: {
                                       indentStream(sstream, indent+3);
                                       sstream << (*i).d_key;
                                       indentStream(sstream, width-(*i).d_key.length());
                                       sstream << " = ";
                                       const int n = (*i).d_string.size();
                                       for (int j = 0; j < n; j++) {
                                         sstream << "\"" << (*i).d_string[j] << "\"";
                                         if (j < n-1) sstream << ", ";
                                       }
                                       break;
                                     }
        }

        /*
         * Output whether the key was used or default in column 60
         */

        if ((*i).d_type != Database::AMP_DATABASE) {
#ifndef LACKS_SSTREAM
          const int tab = 59 - sstream.str().length();
#else
          const int tab = 59 - sstream.pcount();
#endif
          if (tab > 0) indentStream(sstream, tab);
          if ((*i).d_from_default) {
            sstream << " // from default";
          } else if ((*i).d_accessed) {
            sstream << " // input used";
          } else {
            sstream << " // input not used";
          }

          //            sstream << std::endl << ends;
          sstream << std::endl;
          os << sstream.str();
        }
      }
    }

    /*
     * Finally, output all databases in the current key list
     */

    for (std::list<KeyData>::iterator j = d_keyvalues.begin(); j!=d_keyvalues.end(); j++) {
      if ((*j).d_type == Database::AMP_DATABASE) {
        AMP::shared_ptr<MemoryDatabase> db = (*j).d_database;
        db->printDatabase(os, indent+3, toprint);
      }
    }

    indentStream(os, indent);
    os << "}\n";
#endif
}
}
