#include <cmath>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <valarray>
#include <vector>


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
#include "AMP/vectors/SimpleVector.h"

#include "testMaterialHelpers.h"


// Allow external materials to include additional headers in the test
// Note: this includes 1 additional include header that is passed from the command line:
//   Ex:  -D EXTRA_MATERIAL_HEADER='"materials/FuelMaterial.h"'
#ifdef EXTRA_MATERIAL_HEADER
#include EXTRA_MATERIAL_HEADER
#endif


MatTestResult testMaterial( std::string &name )
{
    using namespace AMP::Materials;
    MatTestResult results;
    results.name = name;

    // create material object
    Material::shared_ptr mat;
    try {
        mat = AMP::voodoo::Factory<AMP::Materials::Material>::instance().create( name );
        results.creationGood = true;
    } catch ( std::exception ) {
        results.creationGood = false;
    } catch ( ... ) {
        results.unknown = true;
    }

    // check for undefined property
    try {
        mat->property( "RiDiCuLoUs#!$^&*Name" );
    } catch ( std::exception & ) {
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

        auto propname                  = proplist[type];
        results.propResults[type].name = propname;
        auto property                  = mat->property( propname );

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
                    results.propResults[type].params = true;
                else
                    results.propResults[type].params = false;
            } else {
                results.propResults[type].params = true;
            }
        } catch ( std::exception & ) {
            results.propResults[type].params = false;
        } catch ( ... ) {
            results.propResults[type].params  = false;
            results.propResults[type].unknown = true;
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
        std::vector<AMP::LinearAlgebra::Variable::shared_ptr> toosmallVar( nargs );
        std::vector<AMP::LinearAlgebra::Variable::shared_ptr> justrightVar( nargs );
        std::vector<AMP::LinearAlgebra::Variable::shared_ptr> toobigVar( nargs );
        std::vector<AMP::LinearAlgebra::Vector::shared_ptr> toosmallVec( nargs );
        std::vector<AMP::LinearAlgebra::Vector::shared_ptr> justrightVec( nargs );
        std::vector<AMP::LinearAlgebra::Vector::shared_ptr> toobigVec( nargs );
        for ( size_t i = 0; i < nargs; i++ ) {
            std::stringstream istr;
            istr << i;
            toosmallVar[i].reset( new AMP::LinearAlgebra::Variable( "toosmall" + istr.str() ) );
            justrightVar[i].reset( new AMP::LinearAlgebra::Variable( "justright" + istr.str() ) );
            toobigVar[i].reset( new AMP::LinearAlgebra::Variable( "toobig" + istr.str() ) );
            toosmallVec[i] =
                AMP::LinearAlgebra::SimpleVector<double>::create( npoints, toosmallVar[i] );
            justrightVec[i] =
                AMP::LinearAlgebra::SimpleVector<double>::create( npoints, justrightVar[i] );
            toobigVec[i] =
                AMP::LinearAlgebra::SimpleVector<double>::create( npoints, toobigVar[i] );
            for ( size_t j = 0; j < npoints; j++ ) {
                toosmallVec[i]->setValueByLocalID( j, toosmall[i][j] );
                justrightVec[i]->setValueByLocalID( j, justright[i][j] );
                toobigVec[i]->setValueByLocalID( j, toobig[i][j] );
            }
        }

        // set up std::vector arguments to evalv
        std::vector<double> value( npoints ), nominal;
        std::map<std::string, AMP::shared_ptr<std::vector<double>>> args;
        for ( size_t i = 0; i < nargs; i++ ) {
            args.insert( std::make_pair( argnames[i],
                                         AMP::make_shared<std::vector<double>>( justright[i] ) ) );
        }

        // set up AMP::Vector arguments to evalv
        auto valueVar   = std::make_shared<AMP::LinearAlgebra::Variable>( "value" );
        auto valueVec   = AMP::LinearAlgebra::SimpleVector<double>::create( npoints, valueVar );
        auto nominalVar = std::make_shared<AMP::LinearAlgebra::Variable>( "nominal" );
        auto nominalVec = AMP::LinearAlgebra::SimpleVector<double>::create( npoints, nominalVar );
        std::map<std::string, AMP::LinearAlgebra::Vector::shared_ptr> argsVec;
        for ( size_t i = 0; i < nargs; i++ ) {
            argsVec.insert( std::make_pair( argnames[i], justrightVec[i] ) );
        }

        // set up AMP::MultiVector arguments to evalv
        auto argsMultiVecVec =
            AMP::LinearAlgebra::MultiVector::create( "argsMultiVec", AMP_COMM_SELF );
        auto argsMultiVec =
            AMP::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVector>( argsMultiVecVec );
        for ( size_t i = 0; i < nargs; i++ ) {
            argsMultiVec->addVector( toosmallVec[i] );  // extra junk, should be ignored
            argsMultiVec->addVector( justrightVec[i] ); // paydirt
            argsMultiVec->addVector( toobigVec[i] );    // extra junk, should be ignored
        }
        std::map<std::string, std::string> xlator;
        int count = 0;
        for ( size_t i = 0; i < nargs; i++ ) {
            std::string name = justrightVar[i]->getName();
            xlator.insert( std::make_pair( argnames[i], name ) );
            count++;
        }
        auto nominalMultiVar = std::make_shared<AMP::LinearAlgebra::Variable>( "nominalMulti" );
        auto nominalMultiVec =
            AMP::LinearAlgebra::SimpleVector<double>::create( npoints, nominalMultiVar );

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
            results.propResults[type].range = true;

        // test defaults get and set
        try {
            auto prop = property;
            std::vector<double> defin( nargs );
            for ( size_t i = 0; i < nargs; i++ )
                defin[i] = justright[i][0];
            prop->set_defaults( defin );
            std::vector<double> defaults( prop->get_defaults() );
            if ( defaults == defin )
                results.propResults[type].nargeval[0] = true;
            else
                results.propResults[type].nargeval[0] = false;
        } catch ( std::exception & ) {
            results.propResults[type].nargeval[0] = false;
        } catch ( ... ) {
            results.propResults[type].nargeval[0] = false;
            results.propResults[type].unknown     = true;
        }

        /////////////////////////////////////////////////////////////////////////////////////////
        // Scalar Property
        /////////////////////////////////////////////////////////////////////////////////////////
        if ( property->isScalar() ) {

            // all in range, std::vector
            try {
                property->evalv( value, args );
                nominal                              = value;
                results.propResults[type].success[0] = true;
            } catch ( std::exception & ) {
                results.propResults[type].success[0] = false;
            } catch ( ... ) {
                results.propResults[type].success[0] = false;
                results.propResults[type].unknown    = true;
            }

            // first out of range low, std::vector
            if ( !args.empty() ) {
                try {
                    args.find( argnames[0] )->second->operator[]( 5 ) = toosmall[0][5];
                    property->evalv( value, args );
                    args.find( argnames[0] )->second->operator[]( 5 ) = justright[0][5];
                    results.propResults[type].success[1]              = false;
                } catch ( std::exception & ) {
                    results.propResults[type].success[1]              = true;
                    args.find( argnames[0] )->second->operator[]( 5 ) = justright[0][5];
                } catch ( ... ) {
                    results.propResults[type].success[1] = false;
                    results.propResults[type].unknown    = true;
                }
            } else {
                results.propResults[type].success[1] = true;
            }

            // first out of range hi, std::vector
            if ( !args.empty() ) {
                try {
                    args.find( argnames[0] )->second->operator[]( 5 ) = toobig[0][5];
                    property->evalv( value, args );
                    args.find( argnames[0] )->second->operator[]( 5 ) = justright[0][5];
                    results.propResults[type].success[2]              = false;
                } catch ( std::exception & ) {
                    results.propResults[type].success[2]              = true;
                    args.find( argnames[0] )->second->operator[]( 5 ) = justright[0][5];
                } catch ( ... ) {
                    results.propResults[type].success[2] = false;
                    results.propResults[type].unknown    = true;
                }
            } else {
                results.propResults[type].success[2] = true;
            }

            // all in range, AMP::Vector
            try {
                property->evalv( valueVec, argsVec );
                nominalVec->copyVector( valueVec );
                results.propResults[type].success[4] = true;
            } catch ( std::exception & ) {
                results.propResults[type].success[4] = false;
            } catch ( ... ) {
                results.propResults[type].success[4] = false;
                results.propResults[type].unknown    = true;
            }

            // first out of range low, AMP::Vector
            if ( nargs > 0 ) {
                try {
                    argsVec.find( argnames[0] )->second->setValueByLocalID( 5, toosmall[0][5] );
                    property->evalv( valueVec, argsVec );
                    argsVec.find( argnames[0] )->second->setValueByLocalID( 5, justright[0][5] );
                    results.propResults[type].success[5] = false;
                } catch ( std::exception & ) {
                    results.propResults[type].success[5] = true;
                    argsVec.find( argnames[0] )->second->setValueByLocalID( 5, justright[0][5] );
                } catch ( ... ) {
                    results.propResults[type].success[5] = false;
                    results.propResults[type].unknown    = true;
                }
            } else {
                results.propResults[type].success[5] = true;
            }

            // first out of range hi, AMP::Vector
            if ( nargs > 0 ) {
                try {
                    argsVec.find( argnames[0] )->second->setValueByLocalID( 5, toobig[0][5] );
                    property->evalv( valueVec, argsVec );
                    argsVec.find( argnames[0] )->second->setValueByLocalID( 5, justright[0][5] );
                    results.propResults[type].success[6] = false;
                } catch ( std::exception & ) {
                    results.propResults[type].success[6] = true;
                    argsVec.find( argnames[0] )->second->setValueByLocalID( 5, justright[0][5] );
                } catch ( ... ) {
                    results.propResults[type].success[6] = false;
                    results.propResults[type].unknown    = true;
                }
            } else {
                results.propResults[type].success[6] = true;
            }

            // test make_map, first without setting a translator or setting an empty translator
            std::map<std::string, AMP::LinearAlgebra::Vector::shared_ptr> testMap;
            std::map<std::string, std::string> currentXlator = property->get_translator();
            if ( !currentXlator.empty() ) {
                currentXlator.clear();
                property->set_translator( currentXlator );
            }
            bool xlateGood = false;
            if ( nargs > 0 ) {
                try {
                    testMap                              = property->make_map( argsMultiVec );
                    results.propResults[type].success[7] = false;
                } catch ( std::exception & ) {
                    xlateGood = true;
                } catch ( ... ) {
                    results.propResults[type].success[7] = false;
                    results.propResults[type].unknown    = true;
                }
            } else {
                xlateGood = true;
            }
            property->set_translator( xlator );
            std::map<std::string, std::string> testXlatorGet = property->get_translator();
            if ( testXlatorGet == xlator && xlateGood ) {
                results.propResults[type].success[7] = true;
            }

            // test make_map, now with a translator
            try {
                testMap   = property->make_map( argsMultiVec );
                bool good = true;
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
                    results.propResults[type].success[8] = true;
                else
                    results.propResults[type].success[8] = false;
            } catch ( std::exception & ) {
                results.propResults[type].success[8] = false;
            } catch ( ... ) {
                results.propResults[type].success[8] = false;
                results.propResults[type].unknown    = true;
            }

            // all in range, AMP::MultiVector
            try {
                property->evalv( valueVec, argsMultiVec );
                nominalMultiVec->copyVector( valueVec );
                results.propResults[type].success[9] = true;
            } catch ( std::exception & ) {
                results.propResults[type].success[9] = false;
            } catch ( ... ) {
                results.propResults[type].success[9] = false;
                results.propResults[type].unknown    = true;
            }

            // first out of range low, AMP::MultiVector
            if ( nargs > 0 ) {
                try {
                    argsMultiVec->subsetVectorForVariable( justrightVec[0]->getVariable() )
                        ->setValueByLocalID( 5, toosmall[0][5] );
                    property->evalv( valueVec, argsMultiVec );
                    argsMultiVec->subsetVectorForVariable( justrightVec[0]->getVariable() )
                        ->setValueByLocalID( 5, justright[0][5] );
                    results.propResults[type].success[10] = false;
                } catch ( std::exception & ) {
                    results.propResults[type].success[10] = true;
                    argsMultiVec->subsetVectorForVariable( justrightVec[0]->getVariable() )
                        ->setValueByLocalID( 5, justright[0][5] );
                } catch ( ... ) {
                    results.propResults[type].success[10] = false;
                    results.propResults[type].unknown     = true;
                }
            } else {
                results.propResults[type].success[10] = true;
            }

            // first out of range hi, AMP::MultiVector
            if ( nargs > 0 ) {
                try {
                    argsMultiVec->subsetVectorForVariable( justrightVec[0]->getVariable() )
                        ->setValueByLocalID( 5, toobig[0][5] );
                    property->evalv( valueVec, argsMultiVec );
                    argsMultiVec->subsetVectorForVariable( justrightVec[0]->getVariable() )
                        ->setValueByLocalID( 5, justright[0][5] );
                    results.propResults[type].success[11] = false;
                } catch ( std::exception & ) {
                    results.propResults[type].success[11] = true;
                    argsMultiVec->subsetVectorForVariable( justrightVec[0]->getVariable() )
                        ->setValueByLocalID( 5, justright[0][5] );
                } catch ( ... ) {
                    results.propResults[type].success[11] = false;
                    results.propResults[type].unknown     = true;
                }
            } else {
                results.propResults[type].success[11] = true;
            }

            // check vector, Vector, MultiVector all agree
            pass = true;
            if ( results.propResults[type].success[0] && results.propResults[type].success[4] &&
                 results.propResults[type].success[9] ) {
                for ( size_t i = 0; i < npoints; i++ ) {
                    double vstd      = nominal[i];
                    double vVec      = nominalVec->getValueByLocalID( i );
                    double vMultiVec = nominalMultiVec->getValueByLocalID( i );
                    pass             = pass && ( vstd == vVec && vVec == vMultiVec );
                }
                if ( pass )
                    results.propResults[type].success[3] = true;
            } else {
                results.propResults[type].success[3] = false;
            }

            // set up reduced argument list
            std::map<std::string, AMP::shared_ptr<std::vector<double>>> argsm( args );
            if ( nargs > 0 ) {
                auto argend = argsm.end();
                --argend;
                argsm.erase( argend );
            }

            // check that evalv with fewer than normal number of arguments works
            if ( results.propResults[type].success[0] ) {
                try {
                    property->evalv( value, argsm );
                    results.propResults[type].nargeval[1] = true;
                } catch ( std::exception & ) {
                    results.propResults[type].nargeval[1] = false;
                } catch ( ... ) {
                    results.propResults[type].nargeval[1] = false;
                    results.propResults[type].unknown     = true;
                }
            } else {
                results.propResults[type].nargeval[1] = false;
            }


            /////////////////////////////////////////////////////////////////////////////////////////////
            // Vector Property
            /////////////////////////////////////////////////////////////////////////////////////////////
        } else if ( property->isVector() ) {

            results.propResults[type].isVector = true;

            auto vectorProperty =
                AMP::dynamic_pointer_cast<AMP::Materials::VectorProperty<double>>( property );

            // check that scalar nature is not signaled
            if ( vectorProperty->isScalar() ) {
                results.propResults[type].vector[2] = false;
            } else {
                results.propResults[type].vector[2] = true;
            }

            // check scalar evaluator for std::vector disabled
            try {
                vectorProperty->evalv( value, args );
                results.propResults[type].vector[3] = false;
            } catch ( std::exception & ) {
                results.propResults[type].vector[3] = true;
            } catch ( ... ) {
                results.propResults[type].vector[3] = false;
                results.propResults[type].unknown   = true;
            }

            // check scalar evaluator for AMP::Vector disabled
            try {
                vectorProperty->evalv( valueVec, argsVec );
                results.propResults[type].vector[4] = false;
            } catch ( std::exception & ) {
                results.propResults[type].vector[4] = true;
            } catch ( ... ) {
                results.propResults[type].vector[4] = false;
                results.propResults[type].unknown   = true;
            }

            // test make_map, first without setting a translator or setting an empty translator
            std::map<std::string, AMP::LinearAlgebra::Vector::shared_ptr> testMap;
            auto currentXlator = vectorProperty->get_translator();
            if ( !currentXlator.empty() ) {
                currentXlator.clear();
                vectorProperty->set_translator( currentXlator );
            }
            bool xlateGood = false;
            if ( nargs > 0 ) {
                try {
                    testMap                              = vectorProperty->make_map( argsMultiVec );
                    results.propResults[type].success[7] = false;
                } catch ( std::exception & ) {
                    xlateGood = true;
                } catch ( ... ) {
                    results.propResults[type].success[7] = false;
                    results.propResults[type].unknown    = true;
                }
            } else {
                xlateGood = true;
            }
            vectorProperty->set_translator( xlator );
            auto testXlatorGet = vectorProperty->get_translator();
            if ( testXlatorGet == xlator && xlateGood ) {
                results.propResults[type].success[7] = true;
            }

            // test make_map, now with a translator
            try {
                testMap   = vectorProperty->make_map( argsMultiVec );
                bool good = true;
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
                    results.propResults[type].success[8] = true;
                else
                    results.propResults[type].success[8] = false;
            } catch ( std::exception & ) {
                results.propResults[type].success[8] = false;
            } catch ( ... ) {
                results.propResults[type].success[8] = false;
                results.propResults[type].unknown    = true;
            }

            // check scalar evaluator for AMP::MultiVector disabled
            try {
                vectorProperty->evalv( valueVec, argsMultiVec );
                results.propResults[type].vector[5] = false;
            } catch ( std::exception & ) {
                results.propResults[type].vector[5] = true;
            } catch ( ... ) {
                results.propResults[type].vector[5] = false;
                results.propResults[type].unknown   = true;
            }

            // prepare results vector, check for reasonable size info
            size_t nvec = 0;
            try {
                nvec                                = vectorProperty->get_dimension();
                results.propResults[type].vector[0] = true;
            } catch ( std::exception & ) {
                results.propResults[type].vector[0] = false;
            } catch ( ... ) {
                results.propResults[type].vector[0] = false;
                results.propResults[type].unknown   = true;
            }
            std::vector<AMP::shared_ptr<std::vector<double>>> stdEval( nvec );
            std::vector<AMP::shared_ptr<std::vector<double>>> nominalEval( nvec );
            for ( size_t i = 0; i < nvec; i++ ) {
                stdEval[i]     = AMP::make_shared<std::vector<double>>( npoints );
                nominalEval[i] = AMP::make_shared<std::vector<double>>( npoints );
            }

            // check that number of components is positive
            if ( results.propResults[type].vector[0] && nvec > 0 ) {
                results.propResults[type].vector[1] = true;
            } else {
                results.propResults[type].vector[1] = false;
            }

            // all in range, std::vector
            try {
                vectorProperty->evalv( stdEval, args );
                nominalEval                          = stdEval;
                results.propResults[type].success[0] = true;
            } catch ( std::exception & ) {
                results.propResults[type].success[0] = false;
            } catch ( ... ) {
                results.propResults[type].success[0] = false;
                results.propResults[type].unknown    = true;
            }

            // first out of range low, std::vector
            if ( !args.empty() ) {
                try {
                    args.find( argnames[0] )->second->operator[]( 5 ) = toosmall[0][5];
                    vectorProperty->evalv( stdEval, args );
                    args.find( argnames[0] )->second->operator[]( 5 ) = justright[0][5];
                    results.propResults[type].success[1]              = false;
                } catch ( std::exception & ) {
                    results.propResults[type].success[1]              = true;
                    args.find( argnames[0] )->second->operator[]( 5 ) = justright[0][5];
                } catch ( ... ) {
                    results.propResults[type].success[1] = false;
                    results.propResults[type].unknown    = true;
                }
            } else {
                results.propResults[type].success[1] = true;
            }

            // first out of range hi, std::vector
            if ( !args.empty() ) {
                try {
                    args.find( argnames[0] )->second->operator[]( 5 ) = toobig[0][5];
                    vectorProperty->evalv( stdEval, args );
                    args.find( argnames[0] )->second->operator[]( 5 ) = justright[0][5];
                    results.propResults[type].success[2]              = false;
                } catch ( std::exception & ) {
                    results.propResults[type].success[2]              = true;
                    args.find( argnames[0] )->second->operator[]( 5 ) = justright[0][5];
                } catch ( ... ) {
                    results.propResults[type].success[2] = false;
                    results.propResults[type].unknown    = true;
                }
            } else {
                results.propResults[type].success[2] = true;
            }

            // setup AMP::Vector evalv results
            std::vector<AMP::LinearAlgebra::Variable::shared_ptr> ampEvalVar( nvec );
            std::vector<AMP::LinearAlgebra::Vector::shared_ptr> ampEval( nvec );
            std::vector<AMP::LinearAlgebra::Variable::shared_ptr> nominalAmpEvalVar( nvec );
            std::vector<AMP::LinearAlgebra::Vector::shared_ptr> nominalAmpEval( nvec );
            std::vector<AMP::LinearAlgebra::Variable::shared_ptr> nominalMultiEvalVar( nvec );
            std::vector<AMP::LinearAlgebra::Vector::shared_ptr> nominalMultiEval( nvec );
            for ( size_t i = 0; i < nvec; i++ ) {
                std::stringstream istr;
                istr << i;
                ampEvalVar[i].reset( new AMP::LinearAlgebra::Variable( "ampEval" + istr.str() ) );
                ampEval[i] =
                    AMP::LinearAlgebra::SimpleVector<double>::create( npoints, ampEvalVar[i] );
                nominalAmpEvalVar[i].reset(
                    new AMP::LinearAlgebra::Variable( "nominalAmpEval" + istr.str() ) );
                nominalAmpEval[i] = AMP::LinearAlgebra::SimpleVector<double>::create(
                    npoints, nominalAmpEvalVar[i] );
                nominalMultiEvalVar[i].reset(
                    new AMP::LinearAlgebra::Variable( "nominalMultiEval" + istr.str() ) );
                nominalMultiEval[i] = AMP::LinearAlgebra::SimpleVector<double>::create(
                    npoints, nominalMultiEvalVar[i] );
            }

            // all in range, AMP::Vector
            try {
                vectorProperty->evalv( ampEval, argsVec );
                for ( size_t i = 0; i < nvec; i++ )
                    nominalAmpEval[i]->copyVector( ampEval[i] );
                results.propResults[type].success[4] = true;
            } catch ( std::exception & ) {
                results.propResults[type].success[4] = false;
            } catch ( ... ) {
                results.propResults[type].success[4] = false;
                results.propResults[type].unknown    = true;
            }

            // first out of range low, AMP::Vector
            if ( nargs > 0 ) {
                try {
                    argsVec.find( argnames[0] )->second->setValueByLocalID( 5, toosmall[0][5] );
                    vectorProperty->evalv( ampEval, argsVec );
                    argsVec.find( argnames[0] )->second->setValueByLocalID( 5, justright[0][5] );
                    results.propResults[type].success[5] = false;
                } catch ( std::exception & ) {
                    results.propResults[type].success[5] = true;
                    argsVec.find( argnames[0] )->second->setValueByLocalID( 5, justright[0][5] );
                } catch ( ... ) {
                    results.propResults[type].success[5] = false;
                    results.propResults[type].unknown    = true;
                }
            } else {
                results.propResults[type].success[5] = true;
            }

            // first out of range hi, AMP::Vector
            if ( nargs > 0 ) {
                try {
                    argsVec.find( argnames[0] )->second->setValueByLocalID( 5, toobig[0][5] );
                    vectorProperty->evalv( ampEval, argsVec );
                    argsVec.find( argnames[0] )->second->setValueByLocalID( 5, justright[0][5] );
                    results.propResults[type].success[6] = false;
                } catch ( std::exception & ) {
                    results.propResults[type].success[6] = true;
                    argsVec.find( argnames[0] )->second->setValueByLocalID( 5, justright[0][5] );
                } catch ( ... ) {
                    results.propResults[type].success[6] = false;
                    results.propResults[type].unknown    = true;
                }
            } else {
                results.propResults[type].success[6] = true;
            }

            // all in range, AMP::MultiVector
            try {
                vectorProperty->evalv( ampEval, argsMultiVec );
                for ( size_t i = 0; i < nvec; i++ )
                    nominalMultiEval[i]->copyVector( ampEval[i] );
                results.propResults[type].success[9] = true;
            } catch ( std::exception & ) {
                results.propResults[type].success[9] = false;
            } catch ( ... ) {
                results.propResults[type].success[9] = false;
                results.propResults[type].unknown    = true;
            }

            // first out of range low, AMP::MultiVector
            if ( nargs > 0 ) {
                try {
                    argsMultiVec->subsetVectorForVariable( justrightVec[0]->getVariable() )
                        ->setValueByLocalID( 5, toosmall[0][5] );
                    vectorProperty->evalv( ampEval, argsMultiVec );
                    argsMultiVec->subsetVectorForVariable( justrightVec[0]->getVariable() )
                        ->setValueByLocalID( 5, justright[0][5] );
                    results.propResults[type].success[10] = false;
                } catch ( std::exception & ) {
                    results.propResults[type].success[10] = true;
                    argsMultiVec->subsetVectorForVariable( justrightVec[0]->getVariable() )
                        ->setValueByLocalID( 5, justright[0][5] );
                } catch ( ... ) {
                    results.propResults[type].success[10] = false;
                    results.propResults[type].unknown     = true;
                }
            } else {
                results.propResults[type].success[10] = true;
            }

            // first out of range hi, AMP::MultiVector
            if ( nargs > 0 ) {
                try {
                    argsMultiVec->subsetVectorForVariable( justrightVec[0]->getVariable() )
                        ->setValueByLocalID( 5, toobig[0][5] );
                    vectorProperty->evalv( ampEval, argsMultiVec );
                    argsMultiVec->subsetVectorForVariable( justrightVec[0]->getVariable() )
                        ->setValueByLocalID( 5, justright[0][5] );
                    results.propResults[type].success[11] = false;
                } catch ( std::exception & ) {
                    results.propResults[type].success[11] = true;
                    argsMultiVec->subsetVectorForVariable( justrightVec[0]->getVariable() )
                        ->setValueByLocalID( 5, justright[0][5] );
                } catch ( ... ) {
                    results.propResults[type].success[11] = false;
                    results.propResults[type].unknown     = true;
                }
            } else {
                results.propResults[type].success[11] = true;
            }

            // check vector, Vector, MultiVector all agree
            pass = true;
            if ( results.propResults[type].success[0] && results.propResults[type].success[4] &&
                 results.propResults[type].success[9] ) {
                for ( size_t j = 0; j < nvec; j++ ) {
                    for ( size_t i = 0; i < npoints; i++ ) {
                        double vstd      = ( *nominalEval[j] )[i];
                        double vVec      = nominalAmpEval[j]->getValueByLocalID( i );
                        double vMultiVec = nominalMultiEval[j]->getValueByLocalID( i );
                        pass             = pass && ( vstd == vVec && vVec == vMultiVec );
                    }
                }
                if ( pass )
                    results.propResults[type].success[3] = true;
            } else {
                results.propResults[type].success[3] = false;
            }

            // set up reduced argument list
            std::map<std::string, AMP::shared_ptr<std::vector<double>>> argsm( args );
            if ( nargs > 0 ) {
                auto argend = argsm.end();
                --argend;
                argsm.erase( argend );
            }

            // check that evalv with fewer than normal number of arguments works
            if ( results.propResults[type].success[0] ) {
                try {
                    vectorProperty->evalv( stdEval, argsm );
                    results.propResults[type].nargeval[1] = true;
                } catch ( std::exception & ) {
                    results.propResults[type].nargeval[1] = false;
                } catch ( ... ) {
                    results.propResults[type].nargeval[1] = false;
                    results.propResults[type].unknown     = true;
                }
            } else {
                results.propResults[type].nargeval[1] = false;
            }

            /////////////////////////////////////////////////////////////////////////////////////////////
            // Tensor Property
            /////////////////////////////////////////////////////////////////////////////////////////////
        } else if ( property->isTensor() ) {

            results.propResults[type].isTensor = true;

            auto tensorProperty =
                AMP::dynamic_pointer_cast<AMP::Materials::TensorProperty<double>>( property );

            // check that scalar nature is not signaled
            if ( tensorProperty->isScalar() ) {
                results.propResults[type].tensor[2] = false;
            } else {
                results.propResults[type].tensor[2] = true;
            }

            // check scalar evaluator for std::vector disabled
            try {
                tensorProperty->evalv( value, args );
                results.propResults[type].tensor[3] = false;
            } catch ( std::exception & ) {
                results.propResults[type].tensor[3] = true;
            } catch ( ... ) {
                results.propResults[type].tensor[3] = false;
                results.propResults[type].unknown   = true;
            }

            // check scalar evaluator for AMP::Vector disabled
            try {
                tensorProperty->evalv( valueVec, argsVec );
                results.propResults[type].tensor[4] = false;
            } catch ( std::exception & ) {
                results.propResults[type].tensor[4] = true;
            } catch ( ... ) {
                results.propResults[type].tensor[4] = false;
                results.propResults[type].unknown   = true;
            }

            // test make_map, first without setting a translator or setting an empty translator
            std::map<std::string, AMP::LinearAlgebra::Vector::shared_ptr> testMap;
            std::map<std::string, std::string> currentXlator = tensorProperty->get_translator();
            if ( !currentXlator.empty() ) {
                currentXlator.clear();
                tensorProperty->set_translator( currentXlator );
            }
            bool xlateGood = false;
            if ( nargs > 0 ) {
                try {
                    testMap                              = tensorProperty->make_map( argsMultiVec );
                    results.propResults[type].success[7] = false;
                } catch ( std::exception & ) {
                    xlateGood = true;
                } catch ( ... ) {
                    results.propResults[type].success[7] = false;
                    results.propResults[type].unknown    = true;
                }
            } else {
                xlateGood = true;
            }
            tensorProperty->set_translator( xlator );
            std::map<std::string, std::string> testXlatorGet = tensorProperty->get_translator();
            if ( testXlatorGet == xlator && xlateGood ) {
                results.propResults[type].success[7] = true;
            }

            // test make_map, now with a translator
            try {
                testMap   = tensorProperty->make_map( argsMultiVec );
                bool good = true;
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
                    results.propResults[type].success[8] = true;
                else
                    results.propResults[type].success[8] = false;
            } catch ( std::exception & ) {
                results.propResults[type].success[8] = false;
            } catch ( ... ) {
                results.propResults[type].success[8] = false;
                results.propResults[type].unknown    = true;
            }

            // check scalar evaluator for AMP::MultiVector disabled
            try {
                tensorProperty->evalv( valueVec, argsMultiVec );
                results.propResults[type].tensor[5] = false;
            } catch ( std::exception & ) {
                results.propResults[type].tensor[5] = true;
            } catch ( ... ) {
                results.propResults[type].tensor[5] = false;
                results.propResults[type].unknown   = true;
            }

            // prepare results vector, check for reasonable size info
            std::vector<size_t> nvecs( 2, 0U );
            try {
                nvecs                               = tensorProperty->get_dimensions();
                results.propResults[type].tensor[0] = true;
            } catch ( std::exception & ) {
                results.propResults[type].tensor[0] = false;
            } catch ( ... ) {
                results.propResults[type].tensor[0] = false;
                results.propResults[type].unknown   = true;
            }
            std::vector<std::vector<AMP::shared_ptr<std::vector<double>>>> stdEval(
                nvecs[0], std::vector<AMP::shared_ptr<std::vector<double>>>( nvecs[1] ) );
            std::vector<std::vector<AMP::shared_ptr<std::vector<double>>>> nominalEval(
                nvecs[0], std::vector<AMP::shared_ptr<std::vector<double>>>( nvecs[1] ) );
            for ( size_t i = 0; i < nvecs[0]; i++ )
                for ( size_t j = 0; j < nvecs[1]; j++ ) {
                    stdEval[i][j]     = AMP::make_shared<std::vector<double>>( npoints );
                    nominalEval[i][j] = AMP::make_shared<std::vector<double>>( npoints );
                }

            // check that number of components is positive
            if ( results.propResults[type].tensor[0] && nvecs[0] > 0 && nvecs[1] > 0 &&
                 nvecs.size() == 2 ) {
                results.propResults[type].tensor[1] = true;
            } else {
                results.propResults[type].tensor[1] = false;
            }

            // all in range, std::vector
            try {
                tensorProperty->evalv( stdEval, args );
                nominalEval                          = stdEval;
                results.propResults[type].success[0] = true;
            } catch ( std::exception & ) {
                results.propResults[type].success[0] = false;
            } catch ( ... ) {
                results.propResults[type].success[0] = false;
                results.propResults[type].unknown    = true;
            }

            // first out of range low, std::vector
            if ( !args.empty() ) {
                try {
                    args.find( argnames[0] )->second->operator[]( 5 ) = toosmall[0][5];
                    tensorProperty->evalv( stdEval, args );
                    args.find( argnames[0] )->second->operator[]( 5 ) = justright[0][5];
                    results.propResults[type].success[1]              = false;
                } catch ( std::exception & ) {
                    results.propResults[type].success[1]              = true;
                    args.find( argnames[0] )->second->operator[]( 5 ) = justright[0][5];
                } catch ( ... ) {
                    results.propResults[type].success[1] = false;
                    results.propResults[type].unknown    = true;
                }
            } else {
                results.propResults[type].success[1] = true;
            }

            // first out of range hi, std::vector
            if ( !args.empty() ) {
                try {
                    args.find( argnames[0] )->second->operator[]( 5 ) = toobig[0][5];
                    tensorProperty->evalv( stdEval, args );
                    args.find( argnames[0] )->second->operator[]( 5 ) = justright[0][5];
                    results.propResults[type].success[2]              = false;
                } catch ( std::exception & ) {
                    results.propResults[type].success[2]              = true;
                    args.find( argnames[0] )->second->operator[]( 5 ) = justright[0][5];
                } catch ( ... ) {
                    results.propResults[type].success[2] = false;
                    results.propResults[type].unknown    = true;
                }
            } else {
                results.propResults[type].success[2] = true;
            }

            // setup AMP::Vector evalv results
            std::vector<std::vector<AMP::LinearAlgebra::Variable::shared_ptr>> ampEvalVar(
                nvecs[0], std::vector<AMP::LinearAlgebra::Variable::shared_ptr>( nvecs[1] ) );
            std::vector<std::vector<AMP::LinearAlgebra::Vector::shared_ptr>> ampEval(
                nvecs[0], std::vector<AMP::LinearAlgebra::Vector::shared_ptr>( nvecs[1] ) );
            std::vector<std::vector<AMP::LinearAlgebra::Variable::shared_ptr>> nominalAmpEvalVar(
                nvecs[0], std::vector<AMP::LinearAlgebra::Variable::shared_ptr>( nvecs[1] ) );
            std::vector<std::vector<AMP::LinearAlgebra::Vector::shared_ptr>> nominalAmpEval(
                nvecs[0], std::vector<AMP::LinearAlgebra::Vector::shared_ptr>( nvecs[1] ) );
            std::vector<std::vector<AMP::LinearAlgebra::Variable::shared_ptr>> nominalMultiEvalVar(
                nvecs[0], std::vector<AMP::LinearAlgebra::Variable::shared_ptr>( nvecs[1] ) );
            std::vector<std::vector<AMP::LinearAlgebra::Vector::shared_ptr>> nominalMultiEval(
                nvecs[0], std::vector<AMP::LinearAlgebra::Vector::shared_ptr>( nvecs[1] ) );
            for ( size_t i = 0; i < nvecs[0]; i++ )
                for ( size_t j = 0; j < nvecs[1]; j++ ) {
                    std::stringstream istr;
                    istr << i;
                    ampEvalVar[i][j].reset(
                        new AMP::LinearAlgebra::Variable( "ampEval" + istr.str() ) );
                    ampEval[i][j] = AMP::LinearAlgebra::SimpleVector<double>::create(
                        npoints, ampEvalVar[i][j] );
                    nominalAmpEvalVar[i][j].reset(
                        new AMP::LinearAlgebra::Variable( "nominalAmpEval" + istr.str() ) );
                    nominalAmpEval[i][j] = AMP::LinearAlgebra::SimpleVector<double>::create(
                        npoints, nominalAmpEvalVar[i][j] );
                    nominalMultiEvalVar[i][j].reset(
                        new AMP::LinearAlgebra::Variable( "nominalMultiEval" + istr.str() ) );
                    nominalMultiEval[i][j] = AMP::LinearAlgebra::SimpleVector<double>::create(
                        npoints, nominalMultiEvalVar[i][j] );
                }

            // all in range, AMP::Vector
            try {
                tensorProperty->evalv( ampEval, argsVec );
                for ( size_t i = 0; i < nvecs[0]; i++ )
                    for ( size_t j = 0; j < nvecs[1]; j++ )
                        nominalAmpEval[i][j]->copyVector( ampEval[i][j] );
                results.propResults[type].success[4] = true;
            } catch ( std::exception & ) {
                results.propResults[type].success[4] = false;
            } catch ( ... ) {
                results.propResults[type].success[4] = false;
                results.propResults[type].unknown    = true;
            }

            // first out of range low, AMP::Vector
            if ( nargs > 0 ) {
                try {
                    argsVec.find( argnames[0] )->second->setValueByLocalID( 5, toosmall[0][5] );
                    tensorProperty->evalv( ampEval, argsVec );
                    argsVec.find( argnames[0] )->second->setValueByLocalID( 5, justright[0][5] );
                    results.propResults[type].success[5] = false;
                } catch ( std::exception & ) {
                    results.propResults[type].success[5] = true;
                    argsVec.find( argnames[0] )->second->setValueByLocalID( 5, justright[0][5] );
                } catch ( ... ) {
                    results.propResults[type].success[5] = false;
                    results.propResults[type].unknown    = true;
                }
            } else {
                results.propResults[type].success[5] = true;
            }

            // first out of range hi, AMP::Vector
            if ( nargs > 0 ) {
                try {
                    argsVec.find( argnames[0] )->second->setValueByLocalID( 5, toobig[0][5] );
                    tensorProperty->evalv( ampEval, argsVec );
                    argsVec.find( argnames[0] )->second->setValueByLocalID( 5, justright[0][5] );
                    results.propResults[type].success[6] = false;
                } catch ( std::exception & ) {
                    results.propResults[type].success[6] = true;
                    argsVec.find( argnames[0] )->second->setValueByLocalID( 5, justright[0][5] );
                } catch ( ... ) {
                    results.propResults[type].success[6] = false;
                    results.propResults[type].unknown    = true;
                }
            } else {
                results.propResults[type].success[6] = true;
            }

            // all in range, AMP::MultiVector
            try {
                tensorProperty->evalv( ampEval, argsMultiVec );
                for ( size_t i = 0; i < nvecs[0]; i++ )
                    for ( size_t j = 0; j < nvecs[1]; j++ )
                        nominalMultiEval[i][j]->copyVector( ampEval[i][j] );
                results.propResults[type].success[9] = true;
            } catch ( std::exception & ) {
                results.propResults[type].success[9] = false;
            } catch ( ... ) {
                results.propResults[type].success[9] = false;
                results.propResults[type].unknown    = true;
            }

            // first out of range low, AMP::MultiVector
            if ( nargs > 0 ) {
                try {
                    argsMultiVec->subsetVectorForVariable( justrightVec[0]->getVariable() )
                        ->setValueByLocalID( 5, toosmall[0][5] );
                    tensorProperty->evalv( ampEval, argsMultiVec );
                    argsMultiVec->subsetVectorForVariable( justrightVec[0]->getVariable() )
                        ->setValueByLocalID( 5, justright[0][5] );
                    results.propResults[type].success[10] = false;
                } catch ( std::exception & ) {
                    results.propResults[type].success[10] = true;
                    argsMultiVec->subsetVectorForVariable( justrightVec[0]->getVariable() )
                        ->setValueByLocalID( 5, justright[0][5] );
                } catch ( ... ) {
                    results.propResults[type].success[10] = false;
                    results.propResults[type].unknown     = true;
                }
            } else {
                results.propResults[type].success[10] = true;
            }

            // first out of range hi, AMP::MultiVector
            if ( nargs > 0 ) {
                try {
                    argsMultiVec->subsetVectorForVariable( justrightVec[0]->getVariable() )
                        ->setValueByLocalID( 5, toobig[0][5] );
                    tensorProperty->evalv( ampEval, argsMultiVec );
                    argsMultiVec->subsetVectorForVariable( justrightVec[0]->getVariable() )
                        ->setValueByLocalID( 5, justright[0][5] );
                    results.propResults[type].success[11] = false;
                } catch ( std::exception & ) {
                    results.propResults[type].success[11] = true;
                    argsMultiVec->subsetVectorForVariable( justrightVec[0]->getVariable() )
                        ->setValueByLocalID( 5, justright[0][5] );
                } catch ( ... ) {
                    results.propResults[type].success[11] = false;
                    results.propResults[type].unknown     = true;
                }
            } else {
                results.propResults[type].success[11] = true;
            }

            // check vector, Vector, MultiVector all agree
            pass = true;
            if ( results.propResults[type].success[0] && results.propResults[type].success[4] &&
                 results.propResults[type].success[9] ) {
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
                    results.propResults[type].success[3] = true;
            } else {
                results.propResults[type].success[3] = false;
            }

            // set up reduced argument list
            std::map<std::string, AMP::shared_ptr<std::vector<double>>> argsm( args );
            if ( nargs > 0 ) {
                auto argend = argsm.end();
                --argend;
                argsm.erase( argend );
            }

            // check that evalv with fewer than normal number of arguments works
            if ( results.propResults[type].success[0] ) {
                try {
                    tensorProperty->evalv( stdEval, argsm );
                    results.propResults[type].nargeval[1] = true;
                } catch ( std::exception & ) {
                    results.propResults[type].nargeval[1] = false;
                } catch ( ... ) {
                    results.propResults[type].nargeval[1] = false;
                    results.propResults[type].unknown     = true;
                }
            } else {
                results.propResults[type].nargeval[1] = false;
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

        using namespace AMP::Materials;

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
        } catch ( std::exception &err ) {
            ut.failure( "Caught unknown exception type" );
        }
    }

    ut.report();
    int num_failed = ut.NumFailGlobal();
    ut.reset();

    AMP::AMPManager::shutdown();
    return num_failed;
}
