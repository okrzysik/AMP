/*
 * test_FixedClad.cc
 *
 *  Created on: Mar 11, 2010
 *      Author: bm, gad
 */

#include "AMP/materials/Material.h"
#include "AMP/materials/TensorProperty.h"
#include "AMP/materials/VectorProperty.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"

#include <iostream>
#include <string>
#include <valarray>
using namespace std;

int main( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    using namespace AMP::Materials;

    bool good = true;

    // get material pointer
    Material::shared_ptr mat =
        AMP::voodoo::Factory<AMP::Materials::Material>::instance().create( "FixedClad" );
    auto prop = mat->property( "ThermalConductivity" );

    // test property accessors
    string tcname = prop->get_name();
    string tcsorc = prop->get_source();
    good          = good && tcname == string( "FixedClad" ) + string( "_ThermalConductivity" );
    good          = good && tcsorc == prop->get_source();
    std::cout << "thermal conductivity name is " << tcname << "\n";
    std::cout << "thermal conductivity source is " << tcsorc << "\n";

    // test material accessors
    size_t n                 = 10;
    std::vector<double> &tv  = *new std::vector<double>( n ); // temperature vector
    std::vector<double> &uv  = *new std::vector<double>( n ); // concentration
    std::vector<double> &bv  = *new std::vector<double>( n ); // burnup
    std::vector<double> &prv = *new std::vector<double>( n ); // poisson ratio (result)
    std::vector<double> &tcv = *new std::vector<double>( n ); // thermal conductivity (result)
    for ( size_t i = 0; i < n; i++ ) {
        tv[i] = 563.4 + i / 10.;
        uv[i] = .05 + i / 100., bv[i] = 0. + i * 100;
    } // set input arguments
    std::map<std::string, AMP::shared_ptr<std::vector<double>>> argMap;
    argMap.insert( std::make_pair( "temperance", AMP::shared_ptr<std::vector<double>>( &tv ) ) );
    argMap.insert( std::make_pair( "burningman", AMP::shared_ptr<std::vector<double>>( &bv ) ) );
    argMap.insert( std::make_pair( "upstart", AMP::shared_ptr<std::vector<double>>( &uv ) ) );

    mat->property( "PoissonRatio" )->evalv( prv, argMap );
    mat->property( "ThermalConductivity" )->evalv( tcv, argMap );

    std::vector<AMP::shared_ptr<std::vector<double>>> vfcv( 3 );
    for ( size_t i = 0; i < 3; i++ )
        vfcv[i] = AMP::make_shared<std::vector<double>>( n );

    auto vectorProperty = AMP::dynamic_pointer_cast<AMP::Materials::VectorProperty<double>>(
        mat->property( "VectorFickCoefficient" ) );
    vectorProperty->set_dimension( 3 );
    double vparams[] = { 1.1, 2.2, 3.3 };
    vectorProperty->set_parameters_and_number( vparams, 3 );
    vectorProperty->evalv( vfcv, argMap );

    std::vector<std::vector<AMP::shared_ptr<std::vector<double>>>> tfcv(
        3, std::vector<AMP::shared_ptr<std::vector<double>>>( 3 ) );
    for ( size_t i = 0; i < 3; i++ )
        for ( size_t j = 0; j < 3; j++ )
            tfcv[i][j] = AMP::make_shared<std::vector<double>>( n );

    AMP::shared_ptr<AMP::Materials::TensorProperty<double>> tensorProperty =
        AMP::dynamic_pointer_cast<AMP::Materials::TensorProperty<double>>(
            mat->property( "TensorFickCoefficient" ) );
    tensorProperty->set_dimensions( std::vector<size_t>( 2, 3U ) );
    double tparams[9] = { 1.1, 2.2, 3.3, 11., 22., 33., 111., 222., 333. };
    tensorProperty->set_parameters_and_number( tparams, 9 );
    tensorProperty->evalv( tfcv, argMap );

    std::vector<double> arg( 3 );
    arg[0] = tv[1];
    arg[1] = bv[1];
    arg[2] = uv[1];
    tcv[1] = prop->eval( arg );

    prop->evalv( tcv, argMap );

    good                    = good && AMP::Utilities::approx_equal( tcv[1], tcv[n - 1] );
    good                    = good && AMP::Utilities::approx_equal( tcv[2], tcv[n - 1] );
    valarray<double> params = prop->get_parameters();
    good                    = good && AMP::Utilities::approx_equal( tcv[1], params[0] );

    std::valarray<double> sparams = mat->property( "PoissonRatio" )->get_parameters();
    for ( size_t i = 0; i < n; i++ ) {
        good = good && prv[i] == sparams[0];
    }

    for ( size_t i = 0; i < 3; i++ ) {
        for ( size_t j = 0; j < n; j++ ) {
            double val = ( *vfcv[i] )[j];
            good       = good && val == vparams[i];
        }
    }

    for ( size_t i = 0; i < 3; i++ ) {
        for ( size_t j = 0; j < 3; j++ ) {
            for ( size_t k = 0; k < n; k++ ) {
                double val = ( *tfcv[i][j] )[k];
                good       = good && val == tparams[i * 3 + j];
            }
        }
    }

    if ( good )
        ut.passes( "basic tests of FixedClad" );
    else
        ut.failure( "basic tests of FixedClad" );

    // test parameter change
    double param[1] = { 1.2345 };

    prop->set_parameters( param, 1U );
    double tcs = prop->eval( arg );
    good       = good && AMP::Utilities::approx_equal( tcs, param[0] );

    mat->property( "ThermalConductivity" )->set_parameters( param, 1U );
    mat->property( "ThermalConductivity" )->evalv( tcv, argMap );
    good = good && AMP::Utilities::approx_equal( tcv[0], param[0] );

    if ( good )
        ut.passes( "basic tests of parameterized FixedClad" );
    else
        ut.failure( "basic tests of parameterized FixedClad" );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
