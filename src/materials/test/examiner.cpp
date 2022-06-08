/*
 * examiner.cc
 *
 *  Created on: Sep 17, 2010
 *      Author: gad
 */

#include "AMP/utils/Utilities.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <valarray>
#include <vector>

#include "AMP/materials/Material.h"
#include "AMP/materials/Property.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Factory.h"

// Allow external materials to include additional headers in the test
// Note: this includes 1 additional include header that is passed from the command line:
//   Ex:  -D EXTRA_MATERIAL_HEADER='"materials/FuelMaterial.h"'
#ifdef EXTRA_MATERIAL_HEADER
    #include EXTRA_MATERIAL_HEADER
#endif


/**
 * Examine the values put out by a material model.
 */

size_t nhelp   = 26;
char helpmsg[] = { "usage: examiner -h\n"
                   "	print this message\n"
                   "\n"
                   "usage: examines filename\n"
                   "	read contents of filename and print material property evaluations\n"
                   "\n"
                   "output is to stdout\n"
                   "\n"
                   "Input file has the form: (comments in parentheses)\n"
                   "\n"
                   "Material name\n"
                   "Property name\n"
                   "Count_NAME = value\n"
                   "Low_NAME = value\n"
                   "High_NAME = value\n"
                   "(where NAME is the name of a material argument)\n"
                   "(an alternate specification of evaluation grid is)\n"
                   "Grid_NAME = [value0 value1 ...]\n"
                   "Format = TSV | CSV | Mathematica\n"
                   "(format of the output)"
                   "\n"
                   "You can mix these two forms of grid specs.\n"
                   "If Count_* is specified and Low * and High * are not, then the whole "
                   "material argument range is used.\n"
                   "If Count_* is one, then the midpoint of the range is used.\n"
                   "TSV=tab-separated values\n"
                   "CSV=comma-separated values\n" };

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );

    // help message
    if ( argc == 2 && std::string( argv[1] ) == "-h" ) {
        std::cerr << helpmsg << std::endl;
        return 1;
    }

    // input section
    std::string infile( "inputExaminer" );
    if ( argc > 1 ) {
        infile = std::string( argv[1] );
    }
    auto inDb = AMP::Database::parseInputFile( infile );

    auto format = inDb->getWithDefault<std::string>( "Format", "TSV" );
    AMP_INSIST( format == "TSV" || format == "CSV" || format == "Mathematica",
                "invalid format specified" );

    AMP_INSIST( inDb->keyExists( "Material" ), "must specify material" );
    auto matname = inDb->getString( "Material" );

    AMP_INSIST( inDb->keyExists( "Property" ), "must specify material property" );
    auto propname = inDb->getString( "Property" );

    // Get the material named by matname
    auto material = AMP::Materials::getMaterial( matname );

    // get argument names and ranges
    auto names   = material->property( propname )->get_arguments();
    auto ranges  = material->property( propname )->get_arg_ranges();
    size_t nargs = names.size();
    std::vector<size_t> narg( nargs );
    std::vector<double> lowarg( nargs ), hiarg( nargs );

    // Create a map that will hold the input variable name and a corresponding pointer to a vector
    // of input values
    std::map<std::string, std::shared_ptr<std::vector<double>>> argMap;
    for ( size_t i = 0; i < nargs; i++ ) {
        argMap.insert( std::make_pair( names[i], std::make_shared<std::vector<double>>( 1 ) ) );
    }

    // Fill in the argument value grid
    std::vector<std::vector<double>> args( nargs );
    for ( size_t iarg = 0; iarg < nargs; iarg++ ) {
        auto keyNumber = std::string( "Count_" ) + names[iarg];
        if ( inDb->keyExists( keyNumber ) ) {
            narg[iarg] = inDb->getScalar<int>( keyNumber );
            AMP_INSIST( narg[iarg] >= 1, "must have" + keyNumber + " >= 1" );

            bool haveLow  = inDb->keyExists( std::string( "Low_" ) + names[iarg] );
            bool haveHi   = inDb->keyExists( std::string( "High_" ) + names[iarg] );
            bool haveBoth = haveLow && haveHi;
            if ( haveLow || haveHi )
                AMP_INSIST( haveBoth,
                            std::string( "must specify Low and High " ) + names[iarg] +
                                std::string( " together" ) );
            if ( haveBoth ) {
                lowarg[iarg] = inDb->getScalar<double>( std::string( "Low_" ) + names[iarg] );
                hiarg[iarg]  = inDb->getScalar<double>( std::string( "High_" ) + names[iarg] );
            } else {
                lowarg[iarg] = ranges[iarg][0];
                hiarg[iarg]  = ranges[iarg][1];
            }
            args[iarg].resize( narg[iarg] );
            if ( narg[iarg] == 1 ) {
                args[iarg][0] = .5 * ( lowarg[iarg] + hiarg[iarg] );
            } else {
                for ( size_t i = 0; i < narg[iarg]; i++ )
                    args[iarg][i] =
                        lowarg[iarg] + i * ( hiarg[iarg] - lowarg[iarg] ) / ( narg[iarg] - 1 );
            }
        } else {
            AMP_INSIST( inDb->keyExists( "Grid_" + names[iarg] ),
                        std::string( "must specify a Grid for " ) + names[iarg] );
            args[iarg] = inDb->getVector<double>( "Grid_" + names[iarg] );
            narg[iarg] = args[iarg].size();
        }
    }
    size_t nargTotal = 1;
    std::vector<size_t> argJump( nargs );
    argJump[nargs - 1] = 1;
    for ( size_t i = nargs - 1; i >= 1; i-- ) {
        argJump[i - 1] = argJump[i] * narg[i];
        nargTotal *= narg[i];
    }
    nargTotal *= narg[0];

    // output section, arbitrary dimension loop
    std::string separator;
    if ( format == "TSV" )
        separator = " ";
    if ( format == "CSV" )
        separator = ",";
    if ( format == "Mathematica" )
        separator = ",";
    if ( format == "Mathematica" ) {
        std::cout << "(* material = " << matname << ", property = " << propname << " *)"
                  << std::endl;
        std::cout << "sizes={";
        for ( size_t i = 0; i < nargs; i++ ) {
            std::cout << narg[i];
            if ( i < nargs - 1 )
                std::cout << ",";
        }
        std::cout << "};" << std::endl << std::endl;
        std::cout << "values={" << std::endl;
    }
    for ( size_t m = 0; m < nargTotal; m++ ) {
        std::vector<size_t> indices( nargs );
        indices[0] = m / argJump[0];
        for ( size_t i = 1; i < nargs; i++ )
            indices[i] = ( m - indices[i - 1] * argJump[i - 1] ) / argJump[i];
        for ( size_t i = 0; i < nargs; i++ ) {
            ( *argMap[names[i]] )[0] = args[i][indices[i]];
        }
        std::vector<double> value( 1 );
        auto property = material->property( propname );
        property->evalv( value, {}, argMap );
        if ( format == "Mathematica" )
            std::cout << "{";
        for ( size_t i = 0; i < nargs; i++ )
            std::cout << args[i][indices[i]] << separator;
        std::cout << value[0];
        if ( format == "Mathematica" ) {
            std::cout << "}";
            if ( m < nargTotal - 1 )
                std::cout << ",";
            else
                std::cout << ";";
        }
        std::cout << std::endl;
    }
    if ( format == "Mathematica" )
        std::cout << "};" << std::endl;

    AMP::AMPManager::shutdown();
}
