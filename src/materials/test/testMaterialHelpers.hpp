#ifndef included_AMP_materials_testMaterialHelper_hpp
#define included_AMP_materials_testMaterialHelper_hpp


#include "testMaterialHelpers.h"


inline std::string xlate( bool val )
{
    std::string y( "yes" ), n( "no " ), r;
    r = val ? y : n;
    return r;
}


void MatTestResult::print() const
{
    std::cout << "for material " << name << ": ";
    std::cout << "creation=" << xlate( creationGood ) << " ";
    std::cout << "undefined=" << xlate( undefined ) << " ";
    std::cout << "unknown=" << xlate( unknown ) << " ";
    std::cout << std::endl;
    std::cout << "    property name                           range params nevalv nargsize unknown"
              << std::endl;
    for ( auto &_j : propResults ) {
        std::cout << "    ";
        unsigned int osize = std::cout.width();
        std::cout.width( 29 );
        std::cout << _j.name << " ";
        std::cout.width( osize );
        std::cout << "          ";
        std::cout << xlate( _j.range ) << "   ";
        std::cout << xlate( _j.params ) << "    ";
        unsigned int nsuccess = 0, nargeval = 0;
        for ( bool succes : _j.success )
            if ( succes )
                nsuccess++;
        std::cout << nsuccess << "/" << NSUCCESS << "    ";
        for ( bool k : _j.nargeval )
            if ( k )
                nargeval++;
        std::cout << nargeval << "/" << NARGEVAL << "      ";
        if ( _j.isVector ) {
            unsigned int nvector = 0;
            for ( bool k : _j.vector )
                if ( k )
                    nvector++;
            std::cout << nvector << "/" << NVECTOR << "      ";
        }
        if ( _j.isTensor ) {
            unsigned int ntensor = 0;
            for ( bool k : _j.tensor )
                if ( k )
                    ntensor++;
            std::cout << ntensor << "/" << NTENSOR << "      ";
        }
        std::cout << xlate( _j.unknown ) << "     ";
        std::cout << std::endl;
    }
}
void MatTestResult::record( AMP::UnitTest &ut ) const
{
    std::string msg = "material " + name + " ";
    if ( creationGood )
        ut.passes( msg + "created" );
    else
        ut.failure( msg + "created" );
    for ( const auto &result : propResults ) {
        msg = "material " + name + " property" + " " + result.name + " ";
        if ( result.params )
            ut.passes( msg + "get/set parameters" );
        else
            ut.failure( msg + "get/set parameters" );

        if ( !result.unknown )
            ut.passes( msg + "unknown error" );
        else
            ut.failure( msg + "unknown error" );

        if ( result.range )
            ut.passes( msg + "in_range std::vector" );
        else
            ut.failure( msg + "in_range std::vector" );

        if ( result.success[0] )
            ut.passes( msg + "evalv std::vector" );
        else
            ut.failure( msg + "evalv std::vector" );

        if ( result.success[1] )
            ut.passes( msg + "evalv std::vector out of range lo 1" );
        else
            ut.failure( msg + "evalv std::vector out of range lo 1" );

        if ( result.success[2] )
            ut.passes( msg + "evalv std::vector out of range hi 1" );
        else
            ut.failure( msg + "evalv std::vector out of range hi 1" );

        if ( result.success[4] )
            ut.passes( msg + "evalv AMP::Vector" );
        else
            ut.failure( msg + "evalv AMP::Vector" );

        if ( result.success[5] )
            ut.passes( msg + "evalv AMP::Vector out of range lo 1" );
        else
            ut.failure( msg + "evalv AMP::Vector out of range lo 1" );

        if ( result.success[6] )
            ut.passes( msg + "evalv AMP::Vector out of range hi 1" );
        else
            ut.failure( msg + "evalv AMP::Vector out of range hi 1" );

        if ( result.success[7] )
            ut.passes( msg + "AMP::Multivector translator" );
        else
            ut.failure( msg + "AMP::Multivector translator" );

        if ( result.success[8] )
            ut.passes( msg + "make_map" );
        else
            ut.failure( msg + "make_map" );

        if ( result.success[9] )
            ut.passes( msg + "evalv AMP::MultiVector" );
        else
            ut.failure( msg + "evalv AMP::MultiVector" );

        if ( result.success[10] )
            ut.passes( msg + "evalv AMP::MultiVector out of range lo 1" );
        else
            ut.failure( msg + "evalv AMP::MultiVector out of range lo 1" );

        if ( result.success[11] )
            ut.passes( msg + "evalv AMP::MultiVector out of range hi 1" );
        else
            ut.failure( msg + "evalv AMP::MultiVector out of range hi 1" );

        if ( result.success[3] )
            ut.passes( msg + "evalv agrees std::vector, AMP::Vector, AMP::MultiVector" );
        else
            ut.failure( msg + "evalv agrees std::vector, AMP::Vector, AMP::MultiVector" );

        if ( result.nargeval[0] )
            ut.passes( msg + "set/get defaults" );
        else
            ut.failure( msg + "set/get defaults" );

        if ( result.nargeval[1] )
            ut.passes( msg + "evalv with missing arguments" );
        else
            ut.failure( msg + "evalv with missing arguments" );
        if ( result.isVector ) {
            if ( result.vector[0] )
                ut.passes( msg + "get_dimension() ok" );
            else
                ut.failure( msg + "get_dimension() ok" );

            if ( result.vector[1] )
                ut.passes( msg + "number of components positive" );
            else
                ut.failure( msg + "number of components positive" );

            if ( result.vector[2] )
                ut.passes( msg + "not a scalar" );
            else
                ut.failure( msg + "not a scalar" );

            if ( result.vector[3] )
                ut.passes( msg + "scalar evaluator std::vector disabled" );
            else
                ut.failure( msg + "scalar evaluator std::vector disabled" );

            if ( result.vector[4] )
                ut.passes( msg + "scalar evaluator AMP::Vector disabled" );
            else
                ut.failure( msg + "scalar evaluator AMP::Vector disabled" );

            if ( result.vector[5] )
                ut.passes( msg + "scalar evaluator AMP::Multivector disabled" );
            else
                ut.failure( msg + "scalar evaluator AMP::Multivector disabled" );
        }
    }
}


#endif
