/*
 * examiner.cc
 *
 *  Created on: Sep 17, 2010
 *      Author: gad
 */

#include "utils/Utilities.h"

#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <valarray>
#include <algorithm>
#include <map>
using std::cout;
using std::endl;
using std::exception;
using std::vector;
using std::string;
using std::valarray;
using std::map;

#include "materials/Material.h"
#include "materials/Property.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/Factory.h"

// Allow external materials to include additional headers in the test
// Note: this includes 1 additional include header that is passed from the command line:
//   Ex:  -D EXTRA_MATERIAL_HEADER='"materials/FuelMaterial.h"'
#ifdef EXTRA_MATERIAL_HEADER
    #include EXTRA_MATERIAL_HEADER
#endif

using namespace AMP;

/**
 * Examine the values put out by a material model.
 */

size_t nhelp=26;
string helpmsg[]={
  "usage: examine -h",
  "	print this message",
  "",
  "usage: examine filename",
  "	read contents of filename and print material property evaluations",
  "",
  "Input file has the form: ",
  "",
  "Material name",
  "Property name",
  "NumberOfTemperatures value",
  "LowTemperature value",
  "HighTemperature value",
  "(likewise for concentration and burnup)",
  "(an alternate specification of evaluation grid is)",
  "Temperatures value0 value1 ...",
  "(likewise for concentration and burnup)",
  "Format TSV | CSV | Mathematica",
  "",
  "You can mix these two forms of grid specs.",
  "If NumberOf* is specified and Low* and High* are not, then the whole material argument range is used.",
  "If NumberOf* is one, then the midpoint of the range is used.",
  "TSV=tab-separated values",
  "CSV=comma-separated values",
  "Currently only valid for properties taking 3 arguments.",
  ""
};

int main(int argc, char* argv[]) {
  AMP::AMPManager::startup(argc, argv);

  // help message
  if (argc==2 and string(argv[1])=="-h") {
    for (size_t i=0; i<nhelp; i++) cout << helpmsg[i]<<endl;
    exit(0);
  }

  // input section
  string infile("inputExaminer");
  if (argc>1) {
    infile = string(argv[1]);
  }
  boost::shared_ptr<AMP::InputDatabase> inDb(new AMP::InputDatabase("inDb"));
  AMP::InputManager::getManager()->parseInputFile(infile, inDb);

  string format;
  format = inDb->getStringWithDefault("Format", "TSV");
  AMP_INSIST(format=="TSV" or format=="CSV" or format=="Mathematica", "invalid format specified");

  AMP_INSIST(inDb->keyExists("Material"), "must specify material");
  string matname = inDb->getString("Material");

  AMP_INSIST(inDb->keyExists("Property"), "must specify material property");
  string propname = inDb->getString("Property");

  // Use the material factory to grab an instance of the material named by matname
  AMP::Materials::Material::shared_ptr material = AMP::voodoo::Factory<AMP::Materials::Material>::instance().create(matname);

  // get argument ranges and indices
  vector<vector<double> > ranges = material->property(propname)->get_arg_ranges();
  vector<string>  names = material->property(propname)->get_arguments();
  vector<size_t> index(3);
  for (size_t i=0; i<names.size(); i++) {
    if (names[i] == "temperature")   index[0] = i;
    if (names[i] == "concentration") index[1] = i;
    if (names[i] == "burnup") 	     index[2] = i;
  }

  string ucnames[3]={"Temperature", "Concentration", "Burnup"};
  vector<size_t> narg(3);
  vector<double> lowarg(3), hiarg(3);
	
  // Setup boost shared pointers of vectors for the input values
  std::vector<double> &temperature(*new std::vector<double>(1));
  std::vector<double> &concentration(*new std::vector<double>(1));
  std::vector<double> &burnup(*new std::vector<double>(1));
	   	
  // Create a map that will hold the input variable name and a corresponding boost pointer to a vector of input values
  std::map<std::string, boost::shared_ptr<std::vector<double> > > argMap;
  vector<vector<double> > args(3);

  // Stuff the input vectors that we previously created into the map
  argMap.insert( std::make_pair( "temperature", 	&temperature ) );
  argMap.insert( std::make_pair( "concentration", &concentration ) );
  argMap.insert( std::make_pair( "burnup", 		&burnup ) );
	
  for (size_t iarg=0; iarg<3; iarg++) {
    if (inDb->keyExists(string("NumberOf")+ucnames[iarg]+string("s"))) {
      narg[iarg] = inDb->getInteger(string("NumberOf")+ucnames[iarg]+string("s"));
      AMP_INSIST(narg[iarg]>=1, string("must have NumberOf")+ucnames[iarg]+string("s >= 1"));

      bool haveLow = inDb->keyExists(string("Low")+ucnames[iarg]);
      bool haveHi = inDb->keyExists(string("High")+ucnames[iarg]);
      bool haveBoth = haveLow and haveHi;
      if (haveLow or haveHi) AMP_INSIST(haveBoth, string("must specify Low and High ")+ucnames[iarg]+string(" together"));
      if (haveBoth) {
	lowarg[iarg] = inDb->getDouble(string("Low")+ucnames[iarg]);
	hiarg[iarg]  = inDb->getDouble(string("High")+ucnames[iarg]);
      } else {
	lowarg[iarg]=ranges[index[iarg]][0];
	hiarg[iarg] =ranges[index[iarg]][1];
      }
      args[iarg].resize(narg[iarg]);
      if (narg[iarg]==1) {
	args[iarg][0] = .5*(lowarg[iarg]+ hiarg[iarg]);
      } else {
	for (size_t i=0; i<narg[iarg]; i++) args[iarg][i] =
	  lowarg[iarg] + i*(hiarg[iarg]-lowarg[iarg])/(narg[iarg]-1);
      }
    } else {
      AMP_INSIST(inDb->keyExists(ucnames[iarg]+string("s")), string("must specify a ")+ucnames[iarg]+string(" array"));
      args[iarg] = inDb->getDoubleArray(ucnames[iarg]+string("s"));
      narg[iarg] = args[iarg].size();
    }
  }

  // output section
  string separator;
  if (format=="TSV") separator = "\t";
  if (format=="CSV") separator = ",";
  if (format=="Mathematica") separator = ",";
  if (format=="Mathematica") {
    cout << "(* material = " << matname <<", property = " << propname <<" *)"<<endl;
    cout << "values={";
  }
  for (size_t k=0; k<narg[2]; k++) {
    if (format=="Mathematica") {
      cout << "{" << endl;
    }
    for (size_t j=0; j<narg[1]; j++) {
      if (format=="Mathematica") {
	cout << "{" << endl;
      }
      for (size_t i=0; i<narg[0]; i++) {
	// Initialize a vector of return values for evalv
	std::vector<double> value(1);
	if (format=="Mathematica") {
	  cout << "{";
	}
	temperature[0] 		= args[0][i];
	concentration[0] 	= args[1][j];
	burnup[0]		= args[2][k];
	material->property(propname)->evalv(value, argMap);
				
	cout << temperature[0] << separator << concentration[0] << separator << burnup[0] << separator << value[0];
	if (format=="Mathematica") {
	  cout << "}";
	  if (i<narg[0]-1) cout << ",";
	}
	cout << endl;
      }
      if (format=="Mathematica") {
	cout << "}";
	if (j<narg[1]-1) cout << ",";
      }
      cout <<endl;
    }
    if (format=="Mathematica") {
      cout << "}";
      if (k<narg[2]-1) cout << ",";
      cout <<endl;
    }
  }
  if (format=="Mathematica") {
    cout << "};"<<endl;
  }

  AMP::AMPManager::shutdown();
}
