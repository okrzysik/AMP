/*
 * testMaterial.cc
 *
 *  Created on: Feb 11, 2010
 *      Author: gad
 */

#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <valarray>
#include <map>
#include <sstream>
using std::cout;
using std::endl;
using std::exception;
using std::vector;
using std::string;
using std::valarray;
using std::map;

#include "materials/Material.h"
#include "materials/Property.h"
#include "utils/Factory.h"
#include "utils/Utilities.h"
#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/UnitTest.h"
#include "vectors/SimpleVector.h"
#include "vectors/MultiVector.h"
#include "vectors/CommunicationList.h"

// Allow external materials to include additional headers in the test
// Note: this includes 1 additional include header that is passed from the command line:
//   Ex:  -D EXTRA_MATERIAL_HEADER='"materials/FuelMaterial.h"'
#ifdef EXTRA_MATERIAL_HEADER
    #include EXTRA_MATERIAL_HEADER
#endif

static const size_t NSUCCESS = 12;
static const size_t NARGEVAL = 2;

class PropTestResult {
public:
	PropTestResult() :
		range(false), params(false), unknown(false), name("none") {
		for (size_t i = 0; i < NSUCCESS; i++)
			success[i] = false;
		for (size_t i = 0; i < NARGEVAL; i++)
			nargeval[i] = false;
	}
	bool range;
	bool success[NSUCCESS];
	bool params;
	bool nargeval[NARGEVAL];
	bool unknown;
	string name;
};

class MatTestResult {
public:
	MatTestResult() :
		undefined(false), unknown(false), creationGood(false), name("none") {
	}
	bool undefined;
	bool unknown;
	bool creationGood;
	string name;
	vector<PropTestResult> propResults;
};

MatTestResult testMaterial(string &name) {
	using namespace AMP::Materials;
	MatTestResult results;

	// create material object
	Material::shared_ptr mat;
	try {
		mat	= AMP::voodoo::Factory<AMP::Materials::Material>::instance().create(name);
		results.creationGood = true;
	} catch (std::exception) {
		results.creationGood = false;
	} catch(...) {
		results.unknown = true;
	}

	// check for undefined property
	try {
		mat->property("RiDiCuLoUs#!$^&*Name");
	} catch (std::exception &) {
		results.undefined = true;
	} catch (...) {
		results.undefined = false;
		results.unknown = true;
	}

	// set up communication list for AMP::Vectors
	AMP::AMP_MPI comm(AMP_COMM_SELF);
	AMP::LinearAlgebra::CommunicationListParameters::shared_ptr commlistparam(new AMP::LinearAlgebra::CommunicationListParameters);
	commlistparam->d_comm = comm;
	commlistparam->d_localsize = 1;
	AMP::LinearAlgebra::CommunicationList::shared_ptr commlist(new AMP::LinearAlgebra::CommunicationList(commlistparam));

	// test property evaluations
	vector<string> proplist(mat->list());
	size_t nprop = proplist.size();
	results.propResults.resize(nprop);
	for (size_t type = 0; type < proplist.size(); type++) {

		string propname = proplist[type];
		results.propResults[type].name = propname;
		vector<string> argnames(mat->property(propname)->get_arguments());
		size_t nargs = mat->property(propname)->get_number_arguments();

		// get min and max arg values
		vector<vector<double> > ranges =
				mat->property(propname)->get_arg_ranges();
		const size_t npoints = 10;
		vector<vector<double> > toosmall(nargs, vector<double> (npoints));
		vector<vector<double> > justright(nargs, vector<double> (npoints));
		vector<vector<double> > toobig(nargs, vector<double> (npoints));
		for (size_t j = 0; j < npoints; j++) {
			for (size_t i = 0; i < ranges.size(); i++) {
				if (ranges[i][0] > 0.)
					toosmall[i][j] = 0.;
				else if (ranges[i][0] < 0.)
					toosmall[i][j] = 2. * ranges[i][0];
				else
					toosmall[i][j] = -1.;

				justright[i][j] = .5 * (ranges[i][1] + ranges[i][0]);

				if (ranges[i][1] > 0.)
					toobig[i][j] = 2. * ranges[i][1];
				else if (ranges[i][1] < 0.)
					toobig[i][j] = 0.;
				else
					toobig[i][j] = 1.;
			}
		}

		// set up AMP::SimpleVector versions of above
		std::vector<AMP::LinearAlgebra::Variable::shared_ptr> toosmallVar(nargs);
		std::vector<AMP::LinearAlgebra::Variable::shared_ptr> justrightVar(nargs);
		std::vector<AMP::LinearAlgebra::Variable::shared_ptr> toobigVar(nargs);
		std::vector<AMP::LinearAlgebra::Vector::shared_ptr> toosmallVec(nargs);
		std::vector<AMP::LinearAlgebra::Vector::shared_ptr> justrightVec(nargs);
		std::vector<AMP::LinearAlgebra::Vector::shared_ptr> toobigVec(nargs);
		for (size_t i=0; i<nargs; i++) {
			std::stringstream istr;
			istr << i;
			toosmallVar[i] .reset(new AMP::LinearAlgebra::Variable("toosmall"+istr.str()));
			justrightVar[i].reset(new AMP::LinearAlgebra::Variable("justright"+istr.str()));
			toobigVar[i]   .reset(new AMP::LinearAlgebra::Variable("toobig"+istr.str()));
			toosmallVec[i]  = AMP::LinearAlgebra::SimpleVector::create(npoints, toosmallVar[i]);
			justrightVec[i] = AMP::LinearAlgebra::SimpleVector::create(npoints, justrightVar[i]);
			toobigVec[i]    = AMP::LinearAlgebra::SimpleVector::create(npoints, toobigVar[i]);
			for (size_t j=0; j<npoints; j++) {
				toosmallVec[i] ->setValueByLocalID(j,toosmall[i][j]);
				justrightVec[i]->setValueByLocalID(j,justright[i][j]);
				toobigVec[i]   ->setValueByLocalID(j,toobig[i][j]);
			}
			toosmallVec[i]->setCommunicationList(commlist);
			justrightVec[i]->setCommunicationList(commlist);
			toobigVec[i]->setCommunicationList(commlist);
		}

		// set up std::vector arguments to evalv
		std::vector<double> value(npoints), nominal;
		map<string, boost::shared_ptr<vector<double> > > args;
		for (size_t i = 0; i < nargs; i++) {
			args.insert(
					std::make_pair(argnames[i],
							new vector<double> (justright[i])));
		}

		// set up AMP::Vector arguments to evalv
		AMP::LinearAlgebra::Variable::shared_ptr valueVar(new AMP::LinearAlgebra::Variable("value"));
		AMP::LinearAlgebra::Vector::shared_ptr valueVec = AMP::LinearAlgebra::SimpleVector::create(npoints, valueVar);
		AMP::LinearAlgebra::Variable::shared_ptr nominalVar(new AMP::LinearAlgebra::Variable("nominal"));
		AMP::LinearAlgebra::Vector::shared_ptr nominalVec = AMP::LinearAlgebra::SimpleVector::create(npoints, nominalVar);
		map<string, AMP::LinearAlgebra::Vector::shared_ptr > argsVec;
		for (size_t i = 0; i < nargs; i++) {
			argsVec.insert(std::make_pair(argnames[i], justrightVec[i]));
		}

		// set up AMP::MultiVector arguments to evalv
		AMP::LinearAlgebra::Vector::shared_ptr argsMultiVecVec =
				AMP::LinearAlgebra::MultiVector::create("argsMultiVec", AMP_COMM_SELF);
		boost::shared_ptr<AMP::LinearAlgebra::MultiVector> argsMultiVec =
				boost::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVector>(argsMultiVecVec);
		argsMultiVec->setCommunicationList(commlist);
		for (size_t i=0; i<nargs; i++) {
			argsMultiVec->addVector(toosmallVec[i]); // extra junk, should be ignored
			argsMultiVec->addVector(justrightVec[i]); // paydirt
			argsMultiVec->addVector(toobigVec[i]); // extra junk, should be ignored
		}
		std::map<std::string, std::string> xlator;
		int count = 0;
		for (size_t i=0; i<nargs; i++) {
			std::string name =  justrightVar[i]->getName();
			xlator.insert(std::make_pair(argnames[i], name));
			count ++;
		}
		AMP::LinearAlgebra::Variable::shared_ptr nominalMultiVar(new AMP::LinearAlgebra::Variable("nominalMulti"));
		AMP::LinearAlgebra::Vector::shared_ptr nominalMultiVec = AMP::LinearAlgebra::SimpleVector::create(npoints, nominalMultiVar);

		// test material range functions
		vector<double> range(2);
		bool good = true;
		good = good and nargs == argnames.size();
		for (size_t i = 0; i < argnames.size(); i++) {
			range = mat->property(propname)->get_arg_range(argnames[i]);
			good = good and range[0] <= range[1];
			good = good and range[0] == ranges[i][0] and range[1] == ranges[i][1];
			good = good and     mat->property(propname)->in_range(argnames[i], justright[i]);
			good = good and not mat->property(propname)->in_range(argnames[i], toosmall[i]);
			good = good and not mat->property(propname)->in_range(argnames[i], toobig[i]);
			good = good and     mat->property(propname)->in_range(argnames[i], justright[i][0]);
			good = good and not mat->property(propname)->in_range(argnames[i], toosmall[i][0]);
			good = good and not mat->property(propname)->in_range(argnames[i], toobig[i][0]);

			good = good and     mat->property(propname)->in_range(argnames[i], *justrightVec[i]);
			good = good and not mat->property(propname)->in_range(argnames[i], *toosmallVec[i]);
			good = good and not mat->property(propname)->in_range(argnames[i], *toobigVec[i]);
		}
		if (good)
			results.propResults[type].range = true;

		// all in range
		try {
			mat->property(propname)->evalv(value, args);
			nominal = value;
			results.propResults[type].success[0] = true;
		} catch (std::exception &) {
			results.propResults[type].success[0] = false;
		} catch (...) {
			results.propResults[type].success[0] = false;
			results.propResults[type].unknown = true;
		}

		// first out of range low
		if(args.size()>0) {
			try {
				args.find(argnames[0])->second->operator[](5) = toosmall[0][5];
				mat->property(propname)->evalv(value, args);
				args.find(argnames[0])->second->operator[](5) = justright[0][5];
				results.propResults[type].success[1] = false;
			} catch (std::exception &) {
				results.propResults[type].success[1] = true;
				args.find(argnames[0])->second->operator[](5) = justright[0][5];
			} catch (...) {
				results.propResults[type].success[1] = false;
				results.propResults[type].unknown = true;
			}
		} else {
			results.propResults[type].success[1] = true;
		}

		// first out of range hi
		if(args.size() >0) {
			try {
				args.find(argnames[0])->second->operator[](5) = toobig[0][5];
				mat->property(propname)->evalv(value, args);
				args.find(argnames[0])->second->operator[](5) = justright[0][5];
				results.propResults[type].success[2] = false;
			} catch (std::exception &) {
				results.propResults[type].success[2] = true;
				args.find(argnames[0])->second->operator[](5) = justright[0][5];
			} catch (...) {
				results.propResults[type].success[2] = false;
			results.propResults[type].unknown = true;
			}
		} else {
			results.propResults[type].success[2] = true;
		}

		// all in range, AMP::Vector
		try {
			mat->property(propname)->evalv(valueVec, argsVec);
			nominalVec->copyVector(valueVec);
			results.propResults[type].success[4] = true;
		} catch (std::exception &) {
			results.propResults[type].success[4] = false;
		} catch (...) {
			results.propResults[type].success[4] = false;
			results.propResults[type].unknown = true;
		}

		// first out of range low, AMP::Vector
		if( nargs>0 ) {
			try {
				argsVec.find(argnames[0])->second->setValueByLocalID(5,toosmall[0][5]);
				mat->property(propname)->evalv(valueVec, argsVec);
				argsVec.find(argnames[0])->second->setValueByLocalID(5,justright[0][5]);
				results.propResults[type].success[5] = false;
			} catch (std::exception &) {
				results.propResults[type].success[5] = true;
				argsVec.find(argnames[0])->second->setValueByLocalID(5,justright[0][5]);
			} catch (...) {
				results.propResults[type].success[5] = false;
				results.propResults[type].unknown = true;
			}
		} else {
			results.propResults[type].success[5] = true;
		}

		// first out of range hi, AMP::Vector
		if( nargs>0 ) {
			try {
				argsVec.find(argnames[0])->second->setValueByLocalID(5,toobig[0][5]);
				mat->property(propname)->evalv(valueVec, argsVec);
				argsVec.find(argnames[0])->second->setValueByLocalID(5,justright[0][5]);
				results.propResults[type].success[6] = false;
			} catch (std::exception &) {
				results.propResults[type].success[6] = true;
				argsVec.find(argnames[0])->second->setValueByLocalID(5,justright[0][5]);
			} catch (...) {
				results.propResults[type].success[6] = false;
				results.propResults[type].unknown = true;
			}
		} else {
			results.propResults[type].success[6] = true;
		}

		// test make_map, first without setting a translator or setting an empty translator
		std::map<std::string, AMP::LinearAlgebra::Vector::shared_ptr> testMap;
		std::map<std::string, std::string> currentXlator = mat->property(propname)->get_translator();
		if (not currentXlator.empty()) {
			currentXlator.clear();
			mat->property(propname)->set_translator(currentXlator);
		}
		bool xlateGood = false;
		if ( nargs>0 ) {
			try {
				testMap = mat->property(propname)->make_map(argsMultiVec);
				results.propResults[type].success[7] = false;
			} catch (std::exception &) {
				xlateGood = true;
			} catch (...) {
				results.propResults[type].success[7] = false;
				results.propResults[type].unknown = true;
			}
		} else {
			xlateGood = true;
		}
		mat->property(propname)->set_translator(xlator);
		std::map<std::string, std::string> testXlatorGet = mat->property(propname)->get_translator();
		if (testXlatorGet == xlator and xlateGood) {
			results.propResults[type].success[7] = true;
		}

		// test make_map, now with a translator
		try {
			testMap = mat->property(propname)->make_map(argsMultiVec);
			bool good = true;
			for (size_t i=0; i<nargs; i++) {
				map<std::string, AMP::LinearAlgebra::Vector::shared_ptr>::iterator vec1It = testMap.find(argnames[i]);
				map<std::string, AMP::LinearAlgebra::Vector::shared_ptr>::iterator vec2It = argsVec.find(argnames[i]);
				bool goodIt=true;
				if (vec1It == testMap.end()) {
					goodIt = false; // make_map missed an argument
				}
				if (vec2It == argsVec.end()) {
					AMP_INSIST(false, "argsVec composed incorrectly");
				}
				if (goodIt) {
					good = good and vec1It->second == vec2It->second;
				}
			}
			if (good) results.propResults[type].success[8] = true;
			else results.propResults[type].success[8] = false;
		} catch (std::exception &) {
			results.propResults[type].success[8] = false;
		} catch (...) {
			results.propResults[type].success[8] = false;
			results.propResults[type].unknown = true;
		}

		// all in range, AMP::MultiVector
		try {
			mat->property(propname)->evalv(valueVec, argsMultiVec);
			nominalMultiVec->copyVector(valueVec);
			results.propResults[type].success[9] = true;
		} catch (std::exception &) {
			results.propResults[type].success[9] = false;
		} catch (...) {
			results.propResults[type].success[9] = false;
			results.propResults[type].unknown = true;
		}

		// first out of range low, AMP::MultiVector
		if( nargs>0 ) {
			try {
				argsMultiVec->subsetVectorForVariable(justrightVec[0]->getVariable())->setValueByLocalID(5,toosmall[0][5]);
				mat->property(propname)->evalv(valueVec, argsMultiVec);
				argsMultiVec->subsetVectorForVariable(justrightVec[0]->getVariable())->setValueByLocalID(5,justright[0][5]);
				results.propResults[type].success[10] = false;
			} catch (std::exception &) {
				results.propResults[type].success[10] = true;
				argsMultiVec->subsetVectorForVariable(justrightVec[0]->getVariable())->setValueByLocalID(5,justright[0][5]);
			} catch (...) {
				results.propResults[type].success[10] = false;
				results.propResults[type].unknown = true;
			}
		} else {
			results.propResults[type].success[10] = true;
		}

		// first out of range hi, AMP::MultiVector
		if( nargs>0 ) {
			try {
				argsMultiVec->subsetVectorForVariable(justrightVec[0]->getVariable())->setValueByLocalID(5,toobig[0][5]);
				mat->property(propname)->evalv(valueVec, argsMultiVec);
				argsMultiVec->subsetVectorForVariable(justrightVec[0]->getVariable())->setValueByLocalID(5,justright[0][5]);
				results.propResults[type].success[11] = false;
			} catch (std::exception &) {
				results.propResults[type].success[11] = true;
				argsMultiVec->subsetVectorForVariable(justrightVec[0]->getVariable())->setValueByLocalID(5,justright[0][5]);
			} catch (...) {
				results.propResults[type].success[11] = false;
				results.propResults[type].unknown = true;
			}
		} else {
			results.propResults[type].success[11] = true;
		}

		// check vector, Vector, MultiVector all agree
		good=true;
		for (size_t i=0; i<npoints; i++) {
			double vstd = nominal[i];
			double vVec = nominalVec->getValueByLocalID(i);
			double vMultiVec = nominalMultiVec->getValueByLocalID(i);
			good = good and (vstd == vVec and vVec == vMultiVec);
		}
		if (good) results.propResults[type].success[3] = true;

		// test parameter get and set
		try {
			PropertyPtr prop = mat->property(propname);
			valarray<double> params = prop->get_parameters();
			if (params.size() > 0) {
				params *= 10.0;
				prop->set_parameters(&params[0], params.size());
				params /= 10.0;
				prop->set_parameters(&params[0], params.size());
				valarray<double> nparams = prop->get_parameters();
				bool good = abs(nparams - params).max() <= 1.e-10
						* abs(params).max();
				if (good)
					results.propResults[type].params = true;
				else
					results.propResults[type].params = false;
			} else {
				results.propResults[type].params = true;
			}
		} catch (std::exception &) {
			results.propResults[type].params = false;
		} catch (...) {
			results.propResults[type].params = false;
			results.propResults[type].unknown = true;
		}

		// test defaults get and set
		try {
			PropertyPtr prop = mat->property(propname);
			vector<double> defin(nargs);
			for (size_t i = 0; i < nargs; i++)
				defin[i] = justright[i][0];
			prop->set_defaults(defin);
			vector<double> defaults(prop->get_defaults());
			if (defaults == defin)
				results.propResults[type].nargeval[0] = true;
			else
				results.propResults[type].nargeval[0] = false;
		} catch (std::exception &) {
			results.propResults[type].nargeval[0] = false;
		} catch (...) {
			results.propResults[type].nargeval[0] = false;
			results.propResults[type].unknown = true;
		}

		map<string, boost::shared_ptr<vector<double> > > argsm(args);
		if( nargs>0 ) {
			map<string, boost::shared_ptr<vector<double> > >::iterator argend = argsm.end();
			argend--;
			argsm.erase(argend);
		}

		// check that evalv with fewer than normal number of arguments works
		try {
			mat->property(propname)->evalv(value, argsm);
			results.propResults[type].nargeval[1] = true;
		} catch (std::exception &) {
			results.propResults[type].nargeval[1] = false;
		} catch (...) {
			results.propResults[type].nargeval[1] = false;
			results.propResults[type].unknown = true;
		}
	}

	return results;
}

string xlate(bool val) {
	string y("yes"), n("no "), r;
	r = val ? y : n;
	return r;
}

int main(int argc, char **argv) {
	AMP::AMPManagerProperties amprops;
	amprops.use_MPI_Abort = false;
	AMP::AMPManager::startup(argc, argv, amprops);
	AMP::UnitTest ut;

	bool good = true;

	try {
		using namespace AMP::Materials;

		// test all materials and all properties and print report
		vector<string>
				matlist =
						AMP::voodoo::Factory<AMP::Materials::Material>::instance().getKeys();
		vector<MatTestResult> scoreCard;
		cout
				<< "In the following output the labels have the following meaning:\n"
				<< "creation:  yes = material was created successfully\n"
				<< "undefined: undefined material create was correctly detected\n"
				<< "range:     yes=material range functions were successfully tested\n"
				<< "params:    yes=get and set property parameters successful\n"
				<< "nevalv:     # errors reported/max #, for calls to evalv\n"
				<< "nargsize:  # errors reported/max #, for incorrect number of arguments to evalv\n"
				<< "unknown:   yes=an unknown error occurred during property tests\n\n\n";
		cout << "number of materials = " << matlist.size() << endl;
		cout << "materials = ";
		for (size_t i = 0; i < matlist.size(); ++i)
			cout << matlist[i] << " ";
		cout << endl;
		for (size_t i = 0; i < matlist.size(); i++) {
			MatTestResult score = testMaterial(matlist[i]);
			score.name = matlist[i];
			scoreCard.push_back(score);
			cout << "for material " << matlist[i] << ": ";
			cout << "creation=" << xlate(score.creationGood) << " ";
			cout << "undefined=" << xlate(score.undefined) << " ";
			cout << "unknown=" << xlate(score.unknown) << " ";
			cout << endl;
			cout
					<< "    property name                           range params nevalv nargsize unknown"
					<< endl;
			for (vector<PropTestResult>::iterator j = score.propResults.begin(); j
					!= score.propResults.end(); ++j) {
				cout << "    ";
				unsigned int osize = cout.width();
				cout.width(29);
				cout << j->name << " ";
				cout.width(osize);
				cout << "          ";
				cout << xlate(j->range) << "   ";
				cout << xlate(j->params) << "    ";
				unsigned int nsuccess = 0, nargeval = 0;
				for (size_t k = 0; k < NSUCCESS; k++)
					if (j->success[k])
						nsuccess++;
				cout << nsuccess << "/" << NSUCCESS << "    ";
				for (size_t k = 0; k < NARGEVAL; k++)
					if (j->nargeval[k])
						nargeval++;
				cout << nargeval << "/" << NARGEVAL << "      ";
				cout << xlate(j->unknown) << "     ";
				cout << endl;
			}
			cout << endl << endl << endl;
		}

		size_t maxpassed = 0;
		for (size_t i = 0; i < scoreCard.size(); i++) {
			MatTestResult score = scoreCard[i];
			string msg = "material " + score.name + " ";
			if (score.creationGood)
				ut.passes(msg + "created");
			else
				ut.failure(msg + "created");
			maxpassed += 1;
			for (vector<PropTestResult>::iterator j = score.propResults.begin(); j
					!= score.propResults.end(); ++j) {
				msg = "material " + score.name + " property" + " " + j->name
						+ " ";
				if (j->params)
					ut.passes(msg + "get/set parameters");
				else
					ut.failure(msg + "get/set parameters");

				if (not j->unknown)
					ut.passes(msg + "unknown error");
				else
					ut.failure(msg + "unknown error");

				if (j->range)
					ut.passes(msg + "in_range std::vector");
				else
					ut.failure(msg + "in_range std::vector");

				if (j->success[0])
					ut.passes(msg + "evalv std::vector");
				else
					ut.failure(msg + "evalv std::vector");

				if (j->success[1])
					ut.passes(msg + "evalv std::vector out of range lo 1");
				else
					ut.failure(msg + "evalv std::vector out of range lo 1");

				if (j->success[2])
					ut.passes(msg + "evalv std::vector out of range hi 1");
				else
					ut.failure(msg + "evalv std::vector out of range hi 1");

				if (j->success[4])
					ut.passes(msg + "evalv AMP::Vector");
				else
					ut.failure(msg + "evalv AMP::Vector");

				if (j->success[5])
					ut.passes(msg + "evalv AMP::Vector out of range lo 1");
				else
					ut.failure(msg + "evalv AMP::Vector out of range lo 1");

				if (j->success[6])
					ut.passes(msg + "evalv AMP::Vector out of range hi 1");
				else
					ut.failure(msg + "evalv AMP::Vector out of range hi 1");

				if (j->success[7])
					ut.passes(msg + "AMP::Multivector translator");
				else
					ut.failure(msg + "AMP::Multivector translator");

				if (j->success[8])
					ut.passes(msg + "make_map");
				else
					ut.failure(msg + "make_map");

				if (j->success[9])
					ut.passes(msg + "evalv AMP::MultiVector");
				else
					ut.failure(msg + "evalv AMP::MultiVector");

				if (j->success[10])
					ut.passes(msg + "evalv AMP::MultiVector out of range lo 1");
				else
					ut.failure(msg + "evalv AMP::MultiVector out of range lo 1");

				if (j->success[11])
					ut.passes(msg + "evalv AMP::MultiVector out of range hi 1");
				else
					ut.failure(msg + "evalv AMP::MultiVector out of range hi 1");

				if (j->success[3])
					ut.passes(msg + "evalv agrees std::vector, AMP::Vector, AMP::MultiVector");
				else
					ut.failure(msg + "evalv agrees std::vector, AMP::Vector, AMP::MultiVector");

				if (j->nargeval[0])
					ut.passes(msg + "set/get defaults");
				else
					ut.failure(msg + "set/get defaults");

				if (j->nargeval[1])
					ut.passes(msg + "evalv with missing arguments");
				else
					ut.failure(msg + "evalv with missing arguments");
				maxpassed += 17;
			}
		}
		cout << endl << endl << endl;

		// check that undefined material name is caught
		try {
			string name("flubber");
			Material::shared_ptr mat = AMP::voodoo::Factory<
					AMP::Materials::Material>::instance().create(name);
			if(mat != NULL) maxpassed += 1;
		} catch (std::exception &err) {
			string msg = err.what();
			bool check = msg == "Unregistered creator";
			good = good and check;
			if (good)
				ut.passes("detected undefined material");
			else
				ut.failure("did not detect undefined material");
			maxpassed += 1;
			good = true;
		}

		cout << endl << endl << endl;
		cout << "number of tests passed = " << ut.NumPassLocal() << "/"
				<< maxpassed << " possible" << endl;
		cout << endl << endl << endl;
	} catch (std::exception &err) {
		cout << "ERROR: While testing " << argv[0] << err.what() << std::endl;
		ut.failure("ERROR: While testing");
	} catch (...) {
		cout << "ERROR: While testing " << argv[0]
				<< "An unknown exception was thrown" << endl;
		ut.failure("ERROR: While testing");
	}

	ut.report(2);

	int num_failed = ut.NumFailGlobal();
	AMP::AMPManager::shutdown();
	return num_failed;
}
