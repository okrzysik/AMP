#include "AMP/TPLs.h"

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
        } else if ( type == SAMRAI::tbox::Database::SAMRAI_DOUBLE ||
                    type == SAMRAI::tbox::Database::SAMRAI_FLOAT ) {
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
std::shared_ptr<SAMRAI::tbox::Database> AMP::Database::cloneToSAMRAI() const
{
    auto db   = std::make_shared<SAMRAI::tbox::MemoryDatabase>( getName() );
    auto keys = getAllKeys();
    for ( const auto &key : keys ) {
        if ( isDatabase( key ) ) {
            auto db2 = getDatabase( key );
            db->putDatabase( key )->copyDatabase( db2->cloneToSAMRAI() );
        } else if ( isType<std::string>( key ) ) {
            db->putStringVector( key, getVector<std::string>( key ) );
        } else if ( isType<bool>( key ) ) {
            db->putBoolVector( key, getVector<bool>( key ) );
        } else if ( isType<int>( key ) ) {
            db->putIntegerVector( key, getVector<int>( key ) );
        } else if ( isType<std::complex<double>>( key ) ) {
            db->putComplexVector( key, getVector<std::complex<double>>( key ) );
        } else if ( isType<double>( key ) ) {
            db->putDoubleVector( key, getVector<double>( key ) );
        } else if ( isType<AMP::DatabaseBox>( key ) ) {
            auto boxes = getVector<DatabaseBox>( key );
            std::vector<SAMRAI::tbox::DatabaseBox> boxes2;
            for ( const auto &box : boxes )
                boxes2.push_back( box.cloneToSAMRAI() );
            db->putDatabaseBoxVector( key, boxes2 );
        } else {
            AMP_ERROR( "Unknown type" );
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
