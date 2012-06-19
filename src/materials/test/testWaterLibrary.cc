/*
 * \file materials/test/testWaterLibrary.cc
 * \brief test of WaterLibrary.cc class
 */


#include <string>
#include <iostream>
#include <valarray>
#include <vector>

#include "materials/Material.h"
#include "utils/Utilities.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"

using namespace std;

int main ( int argc , char **argv )
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

   using namespace AMP::Materials;

   bool good=true;

   try
   {
	   // test constructors
	   boost::shared_ptr<AMP::Materials::Material> mat =
			   AMP::voodoo::Factory<AMP::Materials::Material>::instance().create("WaterLibrary");
	   PropertyPtr prop=mat->property("Temperature");

	   // test property accessors
	   string tcname = prop->get_name();
	   string tcsorc = prop->get_source();
	   good = good and tcname == string("WaterLibrary_Temperature");
	   AMP::pout << "Temperature name is " << tcname << "\n";
	   AMP::pout << "Temperature source is " << tcsorc << "\n";
	   vector<string> args = prop->get_arguments();
	   good = good and args[0] == "enthalpy";
	   good = good and args[1] == "pressure";
	   AMP::pout << "arguments are " << args[0] << " " << args[1] <<"\n";
	   unsigned int nargs = prop->get_number_arguments();
	   good = good and nargs == 2;

	   // test material accessors, all arguments present
	   size_t n=10;
	   boost::shared_ptr<std::vector<double> > enthalpyValue(new std::vector<double>(n));
	   boost::shared_ptr<std::vector<double> > pressureValue(new std::vector<double>(n));
	   vector<double> temperatureValue(n);
	   for (size_t i=0; i<n; i++) 
	   {
			//(*enthalpyValue)[i]=1 + (double) i; 
      //(*pressureValue)[i]=1e5 (1 + (double) i);
			(*enthalpyValue)[i]=1;    // enthalpy 
      (*pressureValue)[i]=1e2;  // pressure
	   }
	   
	   // Block for temporary variables
	   {
		   std::map<std::string, boost::shared_ptr<std::vector<double> > > argMap;
		   argMap.insert( std::make_pair( "enthalpy", enthalpyValue ) );
		   argMap.insert( std::make_pair( "pressure", pressureValue ) );

		   std::vector<double> temperatureValue_mat(temperatureValue);
		   mat->property("Temperature")->evalv(temperatureValue_mat, argMap);
		   prop->evalv(temperatureValue, argMap);

		   // known solution is 1 + 2*H + 3*P
		   double knownSolution = 1 + 2*1 + 3*1e2;
		   if( !AMP::Utilities::approx_equal(temperatureValue[0], knownSolution) ) ut.failure("The answer is wrong.");
 
		   for (size_t i=0; i<n; i++) {good = good and AMP::Utilities::approx_equal(temperatureValue[i],temperatureValue_mat[i]);}
		   for (size_t i=0; i<n; i++) {good = good and AMP::Utilities::approx_equal(temperatureValue[0], temperatureValue[i]);}
	   }

	   // test material accessors, one argument present
	   // This sets the default values for enthalpy and temperature so that if it is not on the argMap, it uses the default.
	   std::vector<double> defaults(2);
	   defaults[0] = 563.4;  // enthalpy
     defaults[1] = 0.05;   // pressure
	   prop->set_defaults(defaults);
	   
	   // Block for temporary variables
	   {
		   // known solution is 1 + 2*H + 3*P
		   double knownSolution = 1 + 2* (*enthalpyValue)[0] + 3*defaults[1];
		   std::map<std::string, boost::shared_ptr<std::vector<double> > > argMap;
		   argMap.insert( std::make_pair( "enthalpy", enthalpyValue ) );
		   std::vector<double> temperatureValue_def(temperatureValue);
		   prop->evalv(temperatureValue_def, argMap);
		   for (size_t i=0; i<n; i++) {good = good and AMP::Utilities::approx_equal(temperatureValue_def[i], knownSolution);}
	   }

	   // test material accessors, no arguments present
	   {
		   // known solution is 1 + 2*H + 3*P
		   double knownSolution = 1 + 2* defaults[0] + 3*defaults[1];
		   std::map<std::string, boost::shared_ptr<std::vector<double> > > argMap;
		   std::vector<double> temperatureValue_def(temperatureValue);
		   prop->evalv(temperatureValue_def, argMap);
		   for (size_t i=0; i<n; i++) {good = good and AMP::Utilities::approx_equal(temperatureValue_def[i], knownSolution);}
	   }

	   if (good) ut.passes("basic tests of Material");
	   else ut.failure("basic tests of Material");
   }
   catch( std::exception &err )
   {
     AMP::pout << "ERROR: While testing " << argv[0] << err.what() << std::endl;
     ut.failure("ERROR: While testing");
   }
   catch( ... )
   {
     AMP::pout << "ERROR: While testing " << argv[0] <<  "An unknown exception was thrown" << endl;
     ut.failure("ERROR: While testing");
   }

   ut.report();

   int num_failed = ut.NumFailGlobal();
   AMP::AMPManager::shutdown();
   return num_failed;
}
