#include "AMP/AMP_TPLs.h"

#ifdef AMP_USE_SAMRAI

    #include "AMP/utils/Database.h"
    #include "AMP/utils/Utilities.h"

    #include "SAMRAI/tbox/InputManager.h"
    #include "SAMRAI/tbox/MemoryDatabase.h"

    #include <complex>


/********************************************************************
 * Construct database from a SAMRAI database                         *
 ********************************************************************/
AMP::Database::Database( SAMRAI::tbox::Database &db ) : d_name( db.getName() )
{
    auto keys = db.getAllKeys();
    for ( const auto &key : keys ) {
        auto type = db.getArrayType( key );
        if ( type == SAMRAI::tbox::Database::SAMRAI_DATABASE ) {
            auto db2 = db.getDatabase( key );
            putDatabase( key, std::make_unique<AMP::Database>( db2 ) );
        } else if ( type == SAMRAI::tbox::Database::SAMRAI_BOOL ) {
            auto data = db.getBoolVector( key );
            putVector( key, data );
        } else if ( type == SAMRAI::tbox::Database::SAMRAI_CHAR ) {
            auto data = db.getCharVector( key );
            putVector( key, data );
        } else if ( type == SAMRAI::tbox::Database::SAMRAI_INT ) {
            auto data = db.getIntegerVector( key );
            putVector( key, data );
        } else if ( type == SAMRAI::tbox::Database::SAMRAI_COMPLEX ) {
            auto data = db.getComplexVector( key );
            putVector( key, data );
        } else if ( type == SAMRAI::tbox::Database::SAMRAI_FLOAT ) {
            auto data = db.getFloatVector( key );
            putVector( key, data );
        } else if ( type == SAMRAI::tbox::Database::SAMRAI_DOUBLE ) {
            auto data = db.getDoubleVector( key );
            putVector( key, data );
        } else if ( type == SAMRAI::tbox::Database::SAMRAI_STRING ) {
            auto data = db.getStringVector( key );
            putVector( key, data );
        } else if ( type == SAMRAI::tbox::Database::SAMRAI_BOX ) {
            auto data = db.getDatabaseBoxVector( key );
            std::vector<AMP::DatabaseBox> data2( data.size() );
            for ( size_t i = 0; i < data.size(); i++ )
                data2[i] = AMP::DatabaseBox( data[i] );
            putVector( key, data2 );
        } else {
            AMP_ERROR( "Unknown type" );
        }
    }
}


/********************************************************************
 * Construct a SAMRAI database                                       *
 ********************************************************************/
static void
addScalar( const AMP::Database &db0, SAMRAI::tbox::Database &db, const std::string &key )
{
    using AMP::getTypeID;
    auto data = db0.getData( key );
    auto id   = data->getDataType();
    if ( id == getTypeID<std::string>() ) {
        db.putString( key, db0.getScalar<std::string>( key ) );
    } else if ( id == getTypeID<bool>() ) {
        db.putBool( key, db0.getScalar<bool>( key ) );
    } else if ( id == getTypeID<char>() ) {
        db.putChar( key, db0.getScalar<char>( key ) );
    } else if ( id == getTypeID<int>() ) {
        db.putInteger( key, db0.getScalar<int>( key ) );
    } else if ( id == getTypeID<float>() ) {
        db.putFloat( key, db0.getScalar<float>( key ) );
    } else if ( id == getTypeID<double>() ) {
        db.putDouble( key, db0.getScalar<double>( key ) );
    } else if ( id == getTypeID<std::complex<double>>() ||
                id == getTypeID<std::complex<float>>() ) {
        db.putComplex( key, db0.getScalar<std::complex<double>>( key ) );
    } else if ( id == getTypeID<AMP::DatabaseBox>() ) {
        auto box = db0.getScalar<AMP::DatabaseBox>( key );
        db.putDatabaseBox( key, box.cloneToSAMRAI() );
    } else if ( data->is_integral() ) {
        db.putInteger( key, db0.getScalar<int>( key ) );
    } else if ( data->is_floating_point() ) {
        db.putDouble( key, db0.getScalar<double>( key ) );
    } else {
        AMP_ERROR( "Unknown type: " + std::string( id.name ) );
    }
}
static void
addVector( const AMP::Database &db0, SAMRAI::tbox::Database &db, const std::string &key )
{
    using AMP::getTypeID;
    auto data = db0.getData( key );
    auto id   = data->getDataType();
    if ( id == getTypeID<std::string>() ) {
        db.putStringVector( key, db0.getVector<std::string>( key ) );
    } else if ( id == getTypeID<bool>() ) {
        db.putBoolVector( key, db0.getVector<bool>( key ) );
    } else if ( id == getTypeID<char>() ) {
        db.putCharVector( key, db0.getVector<char>( key ) );
    } else if ( id == getTypeID<int>() ) {
        db.putIntegerVector( key, db0.getVector<int>( key ) );
    } else if ( id == getTypeID<float>() ) {
        db.putFloatVector( key, db0.getVector<float>( key ) );
    } else if ( id == getTypeID<double>() ) {
        db.putDoubleVector( key, db0.getVector<double>( key ) );
    } else if ( id == getTypeID<std::complex<double>>() ||
                id == getTypeID<std::complex<float>>() ) {
        db.putComplexVector( key, db0.getVector<std::complex<double>>( key ) );
    } else if ( id == getTypeID<AMP::DatabaseBox>() ) {
        auto boxes = db0.getVector<AMP::DatabaseBox>( key );
        std::vector<SAMRAI::tbox::DatabaseBox> boxes2;
        for ( const auto &box : boxes )
            boxes2.push_back( box.cloneToSAMRAI() );
        db.putDatabaseBoxVector( key, boxes2 );
    } else if ( data->is_integral() ) {
        db.putIntegerVector( key, db0.getVector<int>( key ) );
    } else if ( data->is_floating_point() ) {
        db.putDoubleVector( key, db0.getVector<double>( key ) );
    } else {
        AMP_ERROR( "Unknown type: " + std::string( id.name ) );
    }
}
std::shared_ptr<SAMRAI::tbox::Database> AMP::Database::cloneToSAMRAI() const
{
    auto db   = std::make_shared<SAMRAI::tbox::MemoryDatabase>( getName() );
    auto keys = getAllKeys();
    for ( const auto &key : keys ) {
        if ( isDatabase( key ) ) {
            auto db2 = getDatabase( key );
            db->putDatabase( key )->copyDatabase( db2->cloneToSAMRAI() );
        } else {
            auto data = getData( key );
            AMP_ASSERT( data );
            auto size = data->arraySize();
            if ( size.length() == 1 ) {
                addScalar( *this, *db, key );
            } else {
                addVector( *this, *db, key );
            }
        }
    }
    return db;
}


/********************************************************************
 * Construct an AMP database using SAMRAI                            *
 ********************************************************************/
std::shared_ptr<AMP::Database> AMP::Database::readThroughSAMRAI( const std::string &filename )
{
    auto db = std::make_shared<SAMRAI::tbox::MemoryDatabase>( filename );
    SAMRAI::tbox::InputManager::getManager()->parseInputFile( filename, db );
    return std::make_shared<AMP::Database>( db );
}


/********************************************************************
 * Convert DatabaseBox                                               *
 ********************************************************************/
AMP::DatabaseBox::DatabaseBox( const SAMRAI::tbox::DatabaseBox &box ) : d_dim( box.getDimVal() )
{
    d_lower.fill( 0 );
    d_upper.fill( 0 );
    for ( int d = 0; d < d_dim; d++ ) {
        d_lower[d] = box.lower( d );
        d_upper[d] = box.upper( d );
    }
}
SAMRAI::tbox::DatabaseBox AMP::DatabaseBox::cloneToSAMRAI() const
{
    return SAMRAI::tbox::DatabaseBox(
        SAMRAI::tbox::Dimension( d_dim ), d_lower.data(), d_upper.data() );
}


#endif
