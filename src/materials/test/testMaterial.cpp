#include "AMP/materials/Material.h"
#include "AMP/materials/Property.h"
#include "AMP/materials/TensorProperty.h"
#include "AMP/materials/VectorProperty.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Factory.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/VectorBuilder.h"

#include "testMaterialHelpers.h"

#include <cmath>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <valarray>
#include <vector>


// Allow external materials to include additional headers in the test
// Note: this includes 1 additional include header that is passed from the command line:
//   Ex:  -D EXTRA_MATERIAL_HEADER='"materials/FuelMaterial.h"'
#ifdef EXTRA_MATERIAL_HEADER
#include EXTRA_MATERIAL_HEADER
#endif


// Test getting a value out of range
static inline void set( std::vector<double> &vec, size_t index, double value )
{
    vec[index] = value;
}
static inline void set( AMP::LinearAlgebra::Vector &vec, size_t index, double value )
{
    vec.setValuesByLocalID( 1, &index, &value );
}
template<class VEC>
static inline VEC &find( std::map<std::string, VEC> &args, const std::string &var )
{
    return args.find( var )->second;
}
static inline AMP::LinearAlgebra::Vector::shared_ptr
find( std::shared_ptr<AMP::LinearAlgebra::MultiVector> &multivec,
      const std::shared_ptr<AMP::LinearAlgebra::Variable> &var )
{
    return multivec->subsetVectorForVariable( var );
}
template<class PROP, class ARGS, class VAR, class VAL>
void testArgsRange( const PROP &property,
                    ARGS &args,
                    const VAR &var,
                    VAL &value,
                    double outOfRange,
                    double justRight,
                    bool &success,
                    bool &unknown )
{
    auto vec = find( args, var );
    try {
        set( *vec, 5, outOfRange );
        property->evalv( value, args );
        set( *vec, 5, justRight );
        success = false;
    } catch ( const std::exception & ) {
        set( *vec, 5, justRight );
        success = true;
    } catch ( ... ) {
        success = false;
        unknown = true;
    }
}


// Test a given material
MatTestResult testMaterial( std::string &name )
{
    MatTestResult results;
    results.name = name;

    // create material object
    std::shared_ptr<AMP::Materials::Material> mat;
    try {
        mat = AMP::voodoo::Factory<AMP::Materials::Material>::instance().create( name );
        results.creationGood = true;
    } catch ( const std::exception & ) {
        results.creationGood = false;
    } catch ( ... ) {
        results.unknown = true;
    }

    // check for undefined property
    try {
        mat->property( "RiDiCuLoUs#!$^&*Name" );
    } catch ( const std::exception & ) {
        results.undefined = true;
    } catch ( ... ) {
        results.undefined = false;
        results.unknown   = true;
    }

    // test property evaluations
    std::vector<std::string> proplist( mat->list() );
    size_t nprop = proplist.size();
    results.propResults.resize( nprop );
    for ( size_t type = 0; type < proplist.size(); type++ ) {

        auto propname     = proplist[type];
        auto &propResults = results.propResults[type];
        propResults.name  = propname;
        auto property     = mat->property( propname );

        // test parameter get and set
        try {
            auto params = property->get_parameters();
            if ( params.size() > 0 ) {
                params *= 10.0;
                property->set_parameters( &params[0], params.size() );
                params /= 10.0;
                property->set_parameters( &params[0], params.size() );
                auto nparams = property->get_parameters();
                bool good    = abs( nparams - params ).max() <= 1.e-10 * abs( params ).max();
                if ( good )
                    propResults.params = true;
                else
                    propResults.params = false;
            } else {
                propResults.params = true;
            }
        } catch ( const std::exception & ) {
            propResults.params = false;
        } catch ( ... ) {
            propResults.params  = false;
            propResults.unknown = true;
        }

        // get argument info
        std::vector<std::string> argnames( property->get_arguments() );
        size_t nargs = property->get_number_arguments();

        // get min and max arg values
        auto ranges    = property->get_arg_ranges();
        size_t npoints = 10;
        std::vector<std::vector<double>> toosmall( nargs, std::vector<double>( npoints ) );
        std::vector<std::vector<double>> justright( nargs, std::vector<double>( npoints ) );
        std::vector<std::vector<double>> toobig( nargs, std::vector<double>( npoints ) );
        for ( size_t j = 0; j < npoints; j++ ) {
            for ( size_t i = 0; i < ranges.size(); i++ ) {
                if ( ranges[i][0] > 0. )
                    toosmall[i][j] = 0.;
                else if ( ranges[i][0] < 0. )
                    toosmall[i][j] = 2. * ranges[i][0];
                else
                    toosmall[i][j] = -1.;

                justright[i][j] = .5 * ( ranges[i][1] + ranges[i][0] );

                if ( ranges[i][1] > 0. )
                    toobig[i][j] = 2. * ranges[i][1];
                else if ( ranges[i][1] < 0. )
                    toobig[i][j] = 0.;
                else
                    toobig[i][j] = 1.;
            }
        }

        // set up AMP::SimpleVector versions of above
        std::vector<std::string> justrightname( nargs );
        std::vector<AMP::LinearAlgebra::Vector::shared_ptr> toosmallVec( nargs );
        std::vector<AMP::LinearAlgebra::Vector::shared_ptr> justrightVec( nargs );
        std::vector<AMP::LinearAlgebra::Vector::shared_ptr> toobigVec( nargs );
        for ( size_t i = 0; i < nargs; i++ ) {
            std::string istr  = std::to_string( i );
            justrightname[i]  = "justright" + istr;
            auto toosmallVar  = std::make_shared<AMP::LinearAlgebra::Variable>( "toosmall" + istr );
            auto justrightVar = std::make_shared<AMP::LinearAlgebra::Variable>( justrightname[i] );
            auto toobigVar = std::make_shared<AMP::LinearAlgebra::Variable>( "toobigVar" + istr );
            toosmallVec[i] = AMP::LinearAlgebra::createSimpleVector<double>( npoints, toosmallVar );
            justrightVec[i] =
                AMP::LinearAlgebra::createSimpleVector<double>( npoints, justrightVar );
            toobigVec[i] = AMP::LinearAlgebra::createSimpleVector<double>( npoints, toobigVar );
            for ( size_t j = 0; j < npoints; j++ ) {
                toosmallVec[i]->setValuesByLocalID( 1, &j, &toosmall[i][j] );
                justrightVec[i]->setValuesByLocalID( 1, &j, &justright[i][j] );
                toobigVec[i]->setValuesByLocalID( 1, &j, &toobig[i][j] );
            }
        }

        // set up std::vector arguments to evalv
        std::vector<double> value( npoints ), nominal;
        std::map<std::string, std::shared_ptr<std::vector<double>>> args;
        for ( size_t i = 0; i < nargs; i++ ) {
            args.insert( std::make_pair( argnames[i],
                                         std::make_shared<std::vector<double>>( justright[i] ) ) );
        }

        // set up AMP::Vector arguments to evalv
        auto valueVar   = std::make_shared<AMP::LinearAlgebra::Variable>( "value" );
        auto valueVec   = AMP::LinearAlgebra::createSimpleVector<double>( npoints, valueVar );
        auto nominalVar = std::make_shared<AMP::LinearAlgebra::Variable>( "nominal" );
        auto nominalVec = AMP::LinearAlgebra::createSimpleVector<double>( npoints, nominalVar );
        std::map<std::string, AMP::LinearAlgebra::Vector::shared_ptr> argsVec;
        for ( size_t i = 0; i < nargs; i++ ) {
            argsVec.insert( std::make_pair( argnames[i], justrightVec[i] ) );
        }

        // set up AMP::MultiVector arguments to evalv
        auto argsMultiVec =
            AMP::LinearAlgebra::MultiVector::create( "argsMultiVec", AMP_COMM_SELF );
        for ( size_t i = 0; i < nargs; i++ ) {
            argsMultiVec->addVector( toosmallVec[i] );  // extra junk, should be ignored
            argsMultiVec->addVector( justrightVec[i] ); // paydirt
            argsMultiVec->addVector( toobigVec[i] );    // extra junk, should be ignored
        }
        std::map<std::string, std::string> xlator;
        int count = 0;
        for ( size_t i = 0; i < nargs; i++ ) {
            std::string jname = justrightname[i];
            xlator.insert( std::make_pair( argnames[i], jname ) );
            count++;
        }
        auto nominalMultiVar = std::make_shared<AMP::LinearAlgebra::Variable>( "nominalMulti" );
        auto nominalMultiVec =
            AMP::LinearAlgebra::createSimpleVector<double>( npoints, nominalMultiVar );

        // test material range functions
        std::vector<double> range( 2 );
        bool pass = nargs == argnames.size();
        for ( size_t i = 0; i < argnames.size(); i++ ) {
            range = property->get_arg_range( argnames[i] );
            pass  = pass && range[0] <= range[1];
            pass  = pass && ( range[0] == ranges[i][0] ) && ( range[1] == ranges[i][1] );
            pass  = pass && property->in_range( argnames[i], justright[i] );
            pass  = pass && !property->in_range( argnames[i], toosmall[i] );
            pass  = pass && !property->in_range( argnames[i], toobig[i] );
            pass  = pass && property->in_range( argnames[i], justright[i][0] );
            pass  = pass && !property->in_range( argnames[i], toosmall[i][0] );
            pass  = pass && !property->in_range( argnames[i], toobig[i][0] );
            pass  = pass && property->in_range( argnames[i], *justrightVec[i] );
            pass  = pass && !property->in_range( argnames[i], *toosmallVec[i] );
            pass  = pass && !property->in_range( argnames[i], *toobigVec[i] );
        }
        if ( pass )
            propResults.range = true;

        // test defaults get and set
        try {
            auto prop = property;
            std::vector<double> defin( nargs );
            for ( size_t i = 0; i < nargs; i++ )
                defin[i] = justright[i][0];
            prop->set_defaults( defin );
            std::vector<double> defaults( prop->get_defaults() );
            if ( defaults == defin )
                propResults.nargeval[0] = true;
            else
                propResults.nargeval[0] = false;
        } catch ( const std::exception & ) {
            propResults.nargeval[0] = false;
        } catch ( ... ) {
            propResults.nargeval[0] = false;
            propResults.unknown     = true;
        }

        /////////////////////////////////////////////////////////////////////////////////////////
        // Scalar Property
        /////////////////////////////////////////////////////////////////////////////////////////
        if ( property->isScalar() ) {

            // all in range, std::vector
            try {
                property->evalv( value, args );
                nominal                = value;
                propResults.success[0] = true;
            } catch ( const std::exception & ) {
                propResults.success[0] = false;
            } catch ( ... ) {
                propResults.success[0] = false;
                propResults.unknown    = true;
            }

            // first out of range low/hi, std::vector
            if ( !args.empty() ) {
                testArgsRange( property,
                               args,
                               argnames[0],
                               value,
                               toosmall[0][5],
                               justright[0][5],
                               propResults.success[1],
                               propResults.unknown );
                testArgsRange( property,
                               args,
                               argnames[0],
                               value,
                               toobig[0][5],
                               justright[0][5],
                               propResults.success[2],
                               propResults.unknown );
            } else {
                propResults.success[1] = true;
                propResults.success[2] = true;
            }

            // all in range, AMP::Vector
            try {
                property->evalv( valueVec, argsVec );
                nominalVec->copyVector( valueVec );
                propResults.success[4] = true;
            } catch ( const std::exception & ) {
                propResults.success[4] = false;
            } catch ( ... ) {
                propResults.success[4] = false;
                propResults.unknown    = true;
            }

            // first out of range low/hi, AMP::Vector
            if ( nargs > 0 ) {
                testArgsRange( property,
                               argsVec,
                               argnames[0],
                               valueVec,
                               toosmall[0][5],
                               justright[0][5],
                               propResults.success[5],
                               propResults.unknown );
                testArgsRange( property,
                               argsVec,
                               argnames[0],
                               valueVec,
                               toobig[0][5],
                               justright[0][5],
                               propResults.success[6],
                               propResults.unknown );
            } else {
                propResults.success[5] = true;
                propResults.success[6] = true;
            }

            // test make_map, first without setting a translator or setting an empty translator
            std::map<std::string, std::string> currentXlator = property->get_translator();
            if ( !currentXlator.empty() ) {
                currentXlator.clear();
                property->set_translator( currentXlator );
            }
            bool xlateGood = false;
            if ( nargs > 0 ) {
                try {
                    auto testMap           = property->make_map( argsMultiVec );
                    propResults.success[7] = false;
                } catch ( const std::exception & ) {
                    xlateGood = true;
                } catch ( ... ) {
                    propResults.success[7] = false;
                    propResults.unknown    = true;
                }
            } else {
                xlateGood = true;
            }
            property->set_translator( xlator );
            auto testXlatorGet = property->get_translator();
            if ( testXlatorGet == xlator && xlateGood ) {
                propResults.success[7] = true;
            }

            // test make_map, now with a translator
            try {
                auto testMap = property->make_map( argsMultiVec );
                bool good    = true;
                for ( size_t i = 0; i < nargs; i++ ) {
                    auto vec1It = testMap.find( argnames[i] );
                    auto vec2It = argsVec.find( argnames[i] );
                    bool goodIt = true;
                    if ( vec1It == testMap.end() ) {
                        goodIt = false; // make_map missed an argument
                    }
                    if ( vec2It == argsVec.end() ) {
                        AMP_INSIST( false, "argsVec composed incorrectly" );
                    }
                    if ( goodIt ) {
                        good = good && vec1It->second == vec2It->second;
                    }
                }
                if ( good )
                    propResults.success[8] = true;
                else
                    propResults.success[8] = false;
            } catch ( const std::exception & ) {
                propResults.success[8] = false;
            } catch ( ... ) {
                propResults.success[8] = false;
                propResults.unknown    = true;
            }

            // all in range, AMP::MultiVector
            try {
                property->evalv( valueVec, argsMultiVec );
                nominalMultiVec->copyVector( valueVec );
                propResults.success[9] = true;
            } catch ( const std::exception & ) {
                propResults.success[9] = false;
            } catch ( ... ) {
                propResults.success[9] = false;
                propResults.unknown    = true;
            }

            // first out of range low/hi, AMP::MultiVector
            if ( nargs > 0 ) {
                testArgsRange( property,
                               argsMultiVec,
                               justrightVec[0]->getVariable(),
                               valueVec,
                               toosmall[0][5],
                               justright[0][5],
                               propResults.success[10],
                               propResults.unknown );
                testArgsRange( property,
                               argsMultiVec,
                               justrightVec[0]->getVariable(),
                               valueVec,
                               toobig[0][5],
                               justright[0][5],
                               propResults.success[11],
                               propResults.unknown );
            } else {
                propResults.success[10] = true;
                propResults.success[11] = true;
            }

            // check vector, Vector, MultiVector all agree
            pass = true;
            if ( propResults.success[0] && propResults.success[4] && propResults.success[9] ) {
                for ( size_t i = 0; i < npoints; i++ ) {
                    double vstd      = nominal[i];
                    double vVec      = nominalVec->getValueByLocalID( i );
                    double vMultiVec = nominalMultiVec->getValueByLocalID( i );
                    pass             = pass && ( vstd == vVec && vVec == vMultiVec );
                }
                if ( pass )
                    propResults.success[3] = true;
            } else {
                propResults.success[3] = false;
            }

            // set up reduced argument list
            std::map<std::string, std::shared_ptr<std::vector<double>>> argsm( args );
            if ( nargs > 0 ) {
                auto argend = argsm.end();
                --argend;
                argsm.erase( argend );
            }

            // check that evalv with fewer than normal number of arguments works
            if ( propResults.success[0] ) {
                try {
                    property->evalv( value, argsm );
                    propResults.nargeval[1] = true;
                } catch ( const std::exception & ) {
                    propResults.nargeval[1] = false;
                } catch ( ... ) {
                    propResults.nargeval[1] = false;
                    propResults.unknown     = true;
                }
            } else {
                propResults.nargeval[1] = false;
            }


            /////////////////////////////////////////////////////////////////////////////////////////////
            // Vector Property
            /////////////////////////////////////////////////////////////////////////////////////////////
        } else if ( property->isVector() ) {

            propResults.isVector = true;

            auto vectorProperty =
                std::dynamic_pointer_cast<AMP::Materials::VectorProperty<double>>( property );

            // check that scalar nature is not signaled
            if ( vectorProperty->isScalar() ) {
                propResults.vector[2] = false;
            } else {
                propResults.vector[2] = true;
            }

            // check scalar evaluator for std::vector disabled
            try {
                vectorProperty->evalv( value, args );
                propResults.vector[3] = false;
            } catch ( const std::exception & ) {
                propResults.vector[3] = true;
            } catch ( ... ) {
                propResults.vector[3] = false;
                propResults.unknown   = true;
            }

            // check scalar evaluator for AMP::Vector disabled
            try {
                vectorProperty->evalv( valueVec, argsVec );
                propResults.vector[4] = false;
            } catch ( const std::exception & ) {
                propResults.vector[4] = true;
            } catch ( ... ) {
                propResults.vector[4] = false;
                propResults.unknown   = true;
            }

            // test make_map, first without setting a translator or setting an empty translator
            auto currentXlator = vectorProperty->get_translator();
            if ( !currentXlator.empty() ) {
                currentXlator.clear();
                vectorProperty->set_translator( currentXlator );
            }
            bool xlateGood = false;
            if ( nargs > 0 ) {
                try {
                    auto testMap           = vectorProperty->make_map( argsMultiVec );
                    propResults.success[7] = false;
                } catch ( const std::exception & ) {
                    xlateGood = true;
                } catch ( ... ) {
                    propResults.success[7] = false;
                    propResults.unknown    = true;
                }
            } else {
                xlateGood = true;
            }
            vectorProperty->set_translator( xlator );
            auto testXlatorGet = vectorProperty->get_translator();
            if ( testXlatorGet == xlator && xlateGood ) {
                propResults.success[7] = true;
            }

            // test make_map, now with a translator
            try {
                auto testMap = vectorProperty->make_map( argsMultiVec );
                bool good    = true;
                for ( size_t i = 0; i < nargs; i++ ) {
                    auto vec1It = testMap.find( argnames[i] );
                    auto vec2It = argsVec.find( argnames[i] );
                    bool goodIt = true;
                    if ( vec1It == testMap.end() ) {
                        goodIt = false; // make_map missed an argument
                    }
                    if ( vec2It == argsVec.end() ) {
                        AMP_INSIST( false, "argsVec composed incorrectly" );
                    }
                    if ( goodIt ) {
                        good = good && vec1It->second == vec2It->second;
                    }
                }
                if ( good )
                    propResults.success[8] = true;
                else
                    propResults.success[8] = false;
            } catch ( const std::exception & ) {
                propResults.success[8] = false;
            } catch ( ... ) {
                propResults.success[8] = false;
                propResults.unknown    = true;
            }

            // check scalar evaluator for AMP::MultiVector disabled
            try {
                vectorProperty->evalv( valueVec, argsMultiVec );
                propResults.vector[5] = false;
            } catch ( const std::exception & ) {
                propResults.vector[5] = true;
            } catch ( ... ) {
                propResults.vector[5] = false;
                propResults.unknown   = true;
            }

            // prepare results vector, check for reasonable size info
            size_t nvec = 0;
            try {
                nvec                  = vectorProperty->get_dimension();
                propResults.vector[0] = true;
            } catch ( const std::exception & ) {
                propResults.vector[0] = false;
            } catch ( ... ) {
                propResults.vector[0] = false;
                propResults.unknown   = true;
            }
            std::vector<std::shared_ptr<std::vector<double>>> stdEval( nvec );
            std::vector<std::shared_ptr<std::vector<double>>> nominalEval( nvec );
            for ( size_t i = 0; i < nvec; i++ ) {
                stdEval[i]     = std::make_shared<std::vector<double>>( npoints );
                nominalEval[i] = std::make_shared<std::vector<double>>( npoints );
            }

            // check that number of components is positive
            if ( propResults.vector[0] && nvec > 0 ) {
                propResults.vector[1] = true;
            } else {
                propResults.vector[1] = false;
            }

            // all in range, std::vector
            try {
                vectorProperty->evalv( stdEval, args );
                nominalEval            = stdEval;
                propResults.success[0] = true;
            } catch ( const std::exception & ) {
                propResults.success[0] = false;
            } catch ( ... ) {
                propResults.success[0] = false;
                propResults.unknown    = true;
            }

            // first out of range low/hi, std::vector
            if ( !args.empty() ) {
                testArgsRange( vectorProperty,
                               args,
                               argnames[0],
                               stdEval,
                               toosmall[0][5],
                               justright[0][5],
                               propResults.success[1],
                               propResults.unknown );
                testArgsRange( vectorProperty,
                               args,
                               argnames[0],
                               stdEval,
                               toobig[0][5],
                               justright[0][5],
                               propResults.success[2],
                               propResults.unknown );
            } else {
                propResults.success[1] = true;
                propResults.success[2] = true;
            }

            // setup AMP::Vector evalv results
            std::vector<AMP::LinearAlgebra::Vector::shared_ptr> ampEval( nvec );
            std::vector<AMP::LinearAlgebra::Vector::shared_ptr> nominalAmpEval( nvec );
            std::vector<AMP::LinearAlgebra::Vector::shared_ptr> nominalMultiEval( nvec );
            for ( size_t i = 0; i < nvec; i++ ) {
                auto istr     = std::to_string( i );
                auto evalVar  = std::make_shared<AMP::LinearAlgebra::Variable>( "amp" + istr );
                auto nomVar   = std::make_shared<AMP::LinearAlgebra::Variable>( "nominal" + istr );
                auto multiVar = std::make_shared<AMP::LinearAlgebra::Variable>( "multivec" + istr );
                ampEval[i]    = AMP::LinearAlgebra::createSimpleVector<double>( npoints, evalVar );
                nominalAmpEval[i] =
                    AMP::LinearAlgebra::createSimpleVector<double>( npoints, nomVar );
                nominalMultiEval[i] =
                    AMP::LinearAlgebra::createSimpleVector<double>( npoints, multiVar );
            }

            // all in range, AMP::Vector
            try {
                vectorProperty->evalv( ampEval, argsVec );
                for ( size_t i = 0; i < nvec; i++ )
                    nominalAmpEval[i]->copyVector( ampEval[i] );
                propResults.success[4] = true;
            } catch ( const std::exception & ) {
                propResults.success[4] = false;
            } catch ( ... ) {
                propResults.success[4] = false;
                propResults.unknown    = true;
            }

            // first out of range low/hi, AMP::Vector
            if ( nargs > 0 ) {
                testArgsRange( vectorProperty,
                               argsVec,
                               argnames[0],
                               ampEval,
                               toosmall[0][5],
                               justright[0][5],
                               propResults.success[5],
                               propResults.unknown );
                testArgsRange( vectorProperty,
                               argsVec,
                               argnames[0],
                               ampEval,
                               toobig[0][5],
                               justright[0][5],
                               propResults.success[6],
                               propResults.unknown );
            } else {
                propResults.success[5] = true;
                propResults.success[6] = true;
            }

            // all in range, AMP::MultiVector
            try {
                vectorProperty->evalv( ampEval, argsMultiVec );
                for ( size_t i = 0; i < nvec; i++ )
                    nominalMultiEval[i]->copyVector( ampEval[i] );
                propResults.success[9] = true;
            } catch ( const std::exception & ) {
                propResults.success[9] = false;
            } catch ( ... ) {
                propResults.success[9] = false;
                propResults.unknown    = true;
            }

            // first out of range low/hi, AMP::MultiVector
            if ( nargs > 0 ) {
                testArgsRange( vectorProperty,
                               argsMultiVec,
                               justrightVec[0]->getVariable(),
                               ampEval,
                               toosmall[0][5],
                               justright[0][5],
                               propResults.success[10],
                               propResults.unknown );
                testArgsRange( vectorProperty,
                               argsMultiVec,
                               justrightVec[0]->getVariable(),
                               ampEval,
                               toobig[0][5],
                               justright[0][5],
                               propResults.success[11],
                               propResults.unknown );
            } else {
                propResults.success[10] = true;
                propResults.success[11] = true;
            }

            // check vector, Vector, MultiVector all agree
            pass = true;
            if ( propResults.success[0] && propResults.success[4] && propResults.success[9] ) {
                for ( size_t j = 0; j < nvec; j++ ) {
                    for ( size_t i = 0; i < npoints; i++ ) {
                        double vstd      = ( *nominalEval[j] )[i];
                        double vVec      = nominalAmpEval[j]->getValueByLocalID( i );
                        double vMultiVec = nominalMultiEval[j]->getValueByLocalID( i );
                        pass             = pass && ( vstd == vVec && vVec == vMultiVec );
                    }
                }
                if ( pass )
                    propResults.success[3] = true;
            } else {
                propResults.success[3] = false;
            }

            // set up reduced argument list
            std::map<std::string, std::shared_ptr<std::vector<double>>> argsm( args );
            if ( nargs > 0 ) {
                auto argend = argsm.end();
                --argend;
                argsm.erase( argend );
            }

            // check that evalv with fewer than normal number of arguments works
            if ( propResults.success[0] ) {
                try {
                    vectorProperty->evalv( stdEval, argsm );
                    propResults.nargeval[1] = true;
                } catch ( const std::exception & ) {
                    propResults.nargeval[1] = false;
                } catch ( ... ) {
                    propResults.nargeval[1] = false;
                    propResults.unknown     = true;
                }
            } else {
                propResults.nargeval[1] = false;
            }

            /////////////////////////////////////////////////////////////////////////////////////////////
            // Tensor Property
            /////////////////////////////////////////////////////////////////////////////////////////////
        } else if ( property->isTensor() ) {

            propResults.isTensor = true;

            auto tensorProperty =
                std::dynamic_pointer_cast<AMP::Materials::TensorProperty<double>>( property );

            // check that scalar nature is not signaled
            if ( tensorProperty->isScalar() ) {
                propResults.tensor[2] = false;
            } else {
                propResults.tensor[2] = true;
            }

            // check scalar evaluator for std::vector disabled
            try {
                tensorProperty->evalv( value, args );
                propResults.tensor[3] = false;
            } catch ( const std::exception & ) {
                propResults.tensor[3] = true;
            } catch ( ... ) {
                propResults.tensor[3] = false;
                propResults.unknown   = true;
            }

            // check scalar evaluator for AMP::Vector disabled
            try {
                tensorProperty->evalv( valueVec, argsVec );
                propResults.tensor[4] = false;
            } catch ( const std::exception & ) {
                propResults.tensor[4] = true;
            } catch ( ... ) {
                propResults.tensor[4] = false;
                propResults.unknown   = true;
            }

            // test make_map, first without setting a translator or setting an empty translator
            auto currentXlator = tensorProperty->get_translator();
            if ( !currentXlator.empty() ) {
                currentXlator.clear();
                tensorProperty->set_translator( currentXlator );
            }
            bool xlateGood = false;
            if ( nargs > 0 ) {
                try {
                    auto testMap           = tensorProperty->make_map( argsMultiVec );
                    propResults.success[7] = false;
                } catch ( const std::exception & ) {
                    xlateGood = true;
                } catch ( ... ) {
                    propResults.success[7] = false;
                    propResults.unknown    = true;
                }
            } else {
                xlateGood = true;
            }
            tensorProperty->set_translator( xlator );
            auto testXlatorGet = tensorProperty->get_translator();
            if ( testXlatorGet == xlator && xlateGood ) {
                propResults.success[7] = true;
            }

            // test make_map, now with a translator
            try {
                auto testMap = tensorProperty->make_map( argsMultiVec );
                bool good    = true;
                for ( size_t i = 0; i < nargs; i++ ) {
                    auto vec1It = testMap.find( argnames[i] );
                    auto vec2It = argsVec.find( argnames[i] );
                    bool goodIt = true;
                    if ( vec1It == testMap.end() ) {
                        goodIt = false; // make_map missed an argument
                    }
                    if ( vec2It == argsVec.end() ) {
                        AMP_INSIST( false, "argsVec composed incorrectly" );
                    }
                    if ( goodIt ) {
                        good = good && vec1It->second == vec2It->second;
                    }
                }
                if ( good )
                    propResults.success[8] = true;
                else
                    propResults.success[8] = false;
            } catch ( const std::exception & ) {
                propResults.success[8] = false;
            } catch ( ... ) {
                propResults.success[8] = false;
                propResults.unknown    = true;
            }

            // check scalar evaluator for AMP::MultiVector disabled
            try {
                tensorProperty->evalv( valueVec, argsMultiVec );
                propResults.tensor[5] = false;
            } catch ( const std::exception & ) {
                propResults.tensor[5] = true;
            } catch ( ... ) {
                propResults.tensor[5] = false;
                propResults.unknown   = true;
            }

            // prepare results vector, check for reasonable size info
            std::vector<size_t> nvecs( 2, 0U );
            try {
                nvecs                 = tensorProperty->get_dimensions();
                propResults.tensor[0] = true;
            } catch ( const std::exception & ) {
                propResults.tensor[0] = false;
            } catch ( ... ) {
                propResults.tensor[0] = false;
                propResults.unknown   = true;
            }
            std::vector<std::vector<std::shared_ptr<std::vector<double>>>> stdEval(
                nvecs[0], std::vector<std::shared_ptr<std::vector<double>>>( nvecs[1] ) );
            std::vector<std::vector<std::shared_ptr<std::vector<double>>>> nominalEval(
                nvecs[0], std::vector<std::shared_ptr<std::vector<double>>>( nvecs[1] ) );
            for ( size_t i = 0; i < nvecs[0]; i++ )
                for ( size_t j = 0; j < nvecs[1]; j++ ) {
                    stdEval[i][j]     = std::make_shared<std::vector<double>>( npoints );
                    nominalEval[i][j] = std::make_shared<std::vector<double>>( npoints );
                }

            // check that number of components is positive
            if ( propResults.tensor[0] && nvecs[0] > 0 && nvecs[1] > 0 && nvecs.size() == 2 ) {
                propResults.tensor[1] = true;
            } else {
                propResults.tensor[1] = false;
            }

            // all in range, std::vector
            try {
                tensorProperty->evalv( stdEval, args );
                nominalEval            = stdEval;
                propResults.success[0] = true;
            } catch ( const std::exception & ) {
                propResults.success[0] = false;
            } catch ( ... ) {
                propResults.success[0] = false;
                propResults.unknown    = true;
            }

            // first out of range low/hi, std::vector
            if ( !args.empty() ) {
                testArgsRange( tensorProperty,
                               args,
                               argnames[0],
                               stdEval,
                               toosmall[0][5],
                               justright[0][5],
                               propResults.success[1],
                               propResults.unknown );
                testArgsRange( tensorProperty,
                               args,
                               argnames[0],
                               stdEval,
                               toobig[0][5],
                               justright[0][5],
                               propResults.success[2],
                               propResults.unknown );
            } else {
                propResults.success[1] = true;
                propResults.success[2] = true;
            }

            // setup AMP::Vector evalv results
            std::vector<std::vector<AMP::LinearAlgebra::Vector::shared_ptr>> ampEval( nvecs[0] ),
                nominalAmpEval( nvecs[0] ), nominalMultiEval( nvecs[0] );
            for ( size_t i = 0; i < nvecs[0]; i++ ) {
                ampEval[i].resize( nvecs[1] );
                nominalAmpEval[i].resize( nvecs[1] );
                nominalMultiEval[i].resize( nvecs[1] );
                for ( size_t j = 0; j < nvecs[1]; j++ ) {
                    auto istr    = std::to_string( i );
                    auto evalVar = std::make_shared<AMP::LinearAlgebra::Variable>( "amp" + istr );
                    auto nomVar =
                        std::make_shared<AMP::LinearAlgebra::Variable>( "nominal" + istr );
                    auto multiVar =
                        std::make_shared<AMP::LinearAlgebra::Variable>( "multivec" + istr );
                    ampEval[i][j] =
                        AMP::LinearAlgebra::createSimpleVector<double>( npoints, evalVar );
                    nominalAmpEval[i][j] =
                        AMP::LinearAlgebra::createSimpleVector<double>( npoints, nomVar );
                    nominalMultiEval[i][j] =
                        AMP::LinearAlgebra::createSimpleVector<double>( npoints, multiVar );
                }
            }

            // all in range, AMP::Vector
            try {
                tensorProperty->evalv( ampEval, argsVec );
                for ( size_t i = 0; i < nvecs[0]; i++ )
                    for ( size_t j = 0; j < nvecs[1]; j++ )
                        nominalAmpEval[i][j]->copyVector( ampEval[i][j] );
                propResults.success[4] = true;
            } catch ( const std::exception & ) {
                propResults.success[4] = false;
            } catch ( ... ) {
                propResults.success[4] = false;
                propResults.unknown    = true;
            }

            // first out of range low/hi, AMP::Vector
            if ( !args.empty() ) {
                testArgsRange( tensorProperty,
                               argsVec,
                               argnames[0],
                               ampEval,
                               toosmall[0][5],
                               justright[0][5],
                               propResults.success[5],
                               propResults.unknown );
                testArgsRange( tensorProperty,
                               argsVec,
                               argnames[0],
                               ampEval,
                               toobig[0][5],
                               justright[0][5],
                               propResults.success[6],
                               propResults.unknown );
            } else {
                propResults.success[5] = true;
                propResults.success[6] = true;
            }

            // all in range, AMP::MultiVector
            try {
                tensorProperty->evalv( ampEval, argsMultiVec );
                for ( size_t i = 0; i < nvecs[0]; i++ )
                    for ( size_t j = 0; j < nvecs[1]; j++ )
                        nominalMultiEval[i][j]->copyVector( ampEval[i][j] );
                propResults.success[9] = true;
            } catch ( const std::exception & ) {
                propResults.success[9] = false;
            } catch ( ... ) {
                propResults.success[9] = false;
                propResults.unknown    = true;
            }

            // first out of range low/hi, AMP::MultiVector
            if ( nargs > 0 ) {
                testArgsRange( tensorProperty,
                               argsMultiVec,
                               justrightVec[0]->getVariable(),
                               ampEval,
                               toosmall[0][5],
                               justright[0][5],
                               propResults.success[10],
                               propResults.unknown );
                testArgsRange( tensorProperty,
                               argsMultiVec,
                               justrightVec[0]->getVariable(),
                               ampEval,
                               toobig[0][5],
                               justright[0][5],
                               propResults.success[11],
                               propResults.unknown );
            } else {
                propResults.success[10] = true;
                propResults.success[11] = true;
            }

            // check vector, Vector, MultiVector all agree
            pass = true;
            if ( propResults.success[0] && propResults.success[4] && propResults.success[9] ) {
                for ( size_t k = 0; k < nvecs[0]; k++ ) {
                    for ( size_t j = 0; j < nvecs[1]; j++ ) {
                        for ( size_t i = 0; i < npoints; i++ ) {
                            double vstd      = ( *nominalEval[k][j] )[i];
                            double vVec      = nominalAmpEval[k][j]->getValueByLocalID( i );
                            double vMultiVec = nominalMultiEval[k][j]->getValueByLocalID( i );
                            pass             = pass && ( vstd == vVec && vVec == vMultiVec );
                        }
                    }
                }
                if ( pass )
                    propResults.success[3] = true;
            } else {
                propResults.success[3] = false;
            }

            // set up reduced argument list
            std::map<std::string, std::shared_ptr<std::vector<double>>> argsm( args );
            if ( nargs > 0 ) {
                auto argend = argsm.end();
                --argend;
                argsm.erase( argend );
            }

            // check that evalv with fewer than normal number of arguments works
            if ( propResults.success[0] ) {
                try {
                    tensorProperty->evalv( stdEval, argsm );
                    propResults.nargeval[1] = true;
                } catch ( const std::exception & ) {
                    propResults.nargeval[1] = false;
                } catch ( ... ) {
                    propResults.nargeval[1] = false;
                    propResults.unknown     = true;
                }
            } else {
                propResults.nargeval[1] = false;
            }
        }
    }

    return results;
}


int main( int argc, char **argv )
{
    AMP::AMPManagerProperties amprops;
    amprops.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, amprops );
    AMP::UnitTest ut;

    { // Limit scope


        // test all materials and all properties and print report
        auto matlist = AMP::voodoo::Factory<AMP::Materials::Material>::instance().getKeys();
        std::vector<MatTestResult> scoreCard;
        std::cout
            << "In the following output the labels have the following meaning:\n"
            << "creation:  yes = material was created successfully\n"
            << "undefined: undefined material create was correctly detected\n"
            << "range:     yes=material range functions were successfully tested\n"
            << "params:    yes=get and set property parameters successful\n"
            << "nevalv:     # errors reported/max #, for calls to evalv\n"
            << "nargsize:  # errors reported/max #, for incorrect number of arguments to evalv\n"
            << "unknown:   yes=an unknown error occurred during property tests\n\n\n";
        std::cout << "number of materials = " << matlist.size() << std::endl;
        std::cout << "materials = ";
        for ( auto &elem : matlist )
            std::cout << elem << " ";
        std::cout << std::endl;
        for ( auto &elem : matlist ) {
            auto score = testMaterial( elem );
            scoreCard.push_back( score );
            score.record( ut );
            score.print();
            std::cout << std::endl << std::endl;
        }

        // check that undefined material name is caught
        try {
            std::string name( "flubber" );
            auto mat = AMP::voodoo::Factory<AMP::Materials::Material>::instance().create( name );
            ut.failure( "Failed to catch unknown material" );
        } catch ( const StackTrace::abort_error &err ) {
            if ( err.message == "Unregistered creator" )
                ut.passes( "detected undefined material" );
            else
                ut.failure( "did not detect undefined material" );
        } catch ( const std::exception &err ) {
            ut.failure( "Caught unknown exception type" );
        }
    }

    ut.report();
    int num_failed = ut.NumFailGlobal();
    ut.reset();

    AMP::AMPManager::shutdown();
    return num_failed;
}
