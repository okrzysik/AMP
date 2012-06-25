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
	   // test constructors for temperature
	   boost::shared_ptr<AMP::Materials::Material> mat =
		   AMP::voodoo::Factory<AMP::Materials::Material>::instance().create("WaterLibrary"); // get water library
	PropertyPtr liquidEnthalpyProperty	= mat->property("SaturatedLiquidEnthalpy");	// saturated liquid enthalpy property
	PropertyPtr temperatureProperty		= mat->property("Temperature");			// temperature property
	PropertyPtr volumeProperty		= mat->property("SpecificVolume");		// specific volume property
	PropertyPtr conductivityProperty	= mat->property("ThermalConductivity");		// thermal conductivity
	PropertyPtr viscosityProperty		= mat->property("DynamicViscosity");		// dynamic viscosity

	   // test property accessors for temperature
	   string tcname = temperatureProperty->get_name();
	   string tcsorc = temperatureProperty->get_source();
	   AMP::pout << "\n";
	   good = good and tcname == string("WaterLibrary_Temperature");
	   AMP::pout << "Temperature name is " << tcname << "\n";
	   AMP::pout << "Temperature source is " << tcsorc << "\n";
	   vector<string> args = temperatureProperty->get_arguments();
	   good = good and args[0] == "enthalpy";
	   good = good and args[1] == "pressure";
	   AMP::pout << "Temperature property arguments are " << args[0] << " and " << args[1] <<"\n\n";
	   unsigned int nargs = temperatureProperty->get_number_arguments();
	   good = good and nargs == 2;

	   // test material accessors, all arguments present
	   const size_t nIdentical=10; // size of input and output arrays for comparison with identical values
	   const size_t nTest=3; // size of input and output arrays for comparison with known thermodynamic values
	   boost::shared_ptr<std::vector<double> > enthalpyInput(new std::vector<double>(nTest));	// enthalpy input
	   boost::shared_ptr<std::vector<double> > pressureInput(new std::vector<double>(nTest));	// pressure input
	   boost::shared_ptr<std::vector<double> > temperatureInput(new std::vector<double>(nTest));	// temperature input
	   boost::shared_ptr<std::vector<double> > densityInput(new std::vector<double>(nTest));	// density input
	   boost::shared_ptr<std::vector<double> > enthalpyIdenticalInput(new std::vector<double>(nIdentical));// enthalpy input array with identical values
	   boost::shared_ptr<std::vector<double> > pressureIdenticalInput(new std::vector<double>(nIdentical));	// pressure input array with identical values
	   vector<double> liquidEnthalpyOutput(nTest);	// saturated liquid enthalpy output
	   vector<double> temperatureOutput(nTest);		// temperature output
	   vector<double> volumeOutput(nTest);		// specific volume output
	   vector<double> conductivityOutput(nTest);	// thermal conductivity output
	   vector<double> viscosityOutput(nTest);		// dynamic viscosity output
	   vector<double> temperatureIdenticalOutput(nIdentical);		// temperature output array with identical values
	   for (size_t i=0; i<nIdentical; i++) 
	   {
		(*enthalpyIdenticalInput)[i]=500e3;	// enthalpy: 500 kJ/kg 
		(*pressureIdenticalInput)[i]=1e6;	// pressure: 1 MPa
	   }
	(*enthalpyInput)[0] = 500e3;
	(*enthalpyInput)[1] = 1e6;
	(*enthalpyInput)[2] = 100e3;

	(*pressureInput)[0] = 1e6;
	(*pressureInput)[1] = 15e6;
	(*pressureInput)[2] = 30e3;

	(*temperatureInput)[0] = 400;
	(*temperatureInput)[1] = 600;
	(*temperatureInput)[2] = 300;

	(*densityInput)[0] = 937.871;
	(*densityInput)[1] = 659.388;
	(*densityInput)[2] = 996.526;
	   
	   // Block for temporary variables
	   {
		// saturated liquid enthalpy argument map
		   std::map<std::string, boost::shared_ptr<std::vector<double> > > liquidEnthalpyArgMap;
		   liquidEnthalpyArgMap.insert( std::make_pair( "pressure", pressureInput ) );
		// temperature argument map
		   std::map<std::string, boost::shared_ptr<std::vector<double> > > temperatureArgMap;
		   temperatureArgMap.insert( std::make_pair( "enthalpy", enthalpyInput ) );
		   temperatureArgMap.insert( std::make_pair( "pressure", pressureInput ) );
		// specific volume argument map
		   std::map<std::string, boost::shared_ptr<std::vector<double> > > volumeArgMap;
		   volumeArgMap.insert( std::make_pair( "enthalpy", enthalpyInput ) );
		   volumeArgMap.insert( std::make_pair( "pressure", pressureInput ) );
		// thermal conductivity argument map
		   std::map<std::string, boost::shared_ptr<std::vector<double> > > conductivityArgMap;
		   conductivityArgMap.insert( std::make_pair( "temperature", temperatureInput ) );
		   conductivityArgMap.insert( std::make_pair( "density", densityInput ) );
		// dynamic viscosity argument map
		   std::map<std::string, boost::shared_ptr<std::vector<double> > > viscosityArgMap;
		   viscosityArgMap.insert( std::make_pair( "temperature", temperatureInput ) );
		   viscosityArgMap.insert( std::make_pair( "density", densityInput ) );
		// temperature argument map for identical values case
		   std::map<std::string, boost::shared_ptr<std::vector<double> > > temperatureIdenticalArgMap;
		   temperatureIdenticalArgMap.insert( std::make_pair( "enthalpy", enthalpyIdenticalInput ) );
		   temperatureIdenticalArgMap.insert( std::make_pair( "pressure", pressureIdenticalInput ) );
	
		// evaluate saturated liquid enthalpy function
		   liquidEnthalpyProperty->evalv(liquidEnthalpyOutput, liquidEnthalpyArgMap);
		// evaluate temperature function
		   std::vector<double> temperatureOutput_mat(temperatureIdenticalOutput);
		   mat->property("Temperature")->evalv(temperatureOutput_mat, temperatureIdenticalArgMap);
		   temperatureProperty->evalv(temperatureOutput, temperatureArgMap);
		// evaluate specific volume function
		   volumeProperty->evalv(volumeOutput, volumeArgMap);
		// evaluate thermal conductivity function
		   conductivityProperty->evalv(conductivityOutput, conductivityArgMap);
		// evaluate dynamic viscosity function
		   viscosityProperty->evalv(viscosityOutput, viscosityArgMap);
		// evaluate temperature function for identical values case
		   temperatureProperty->evalv(temperatureIdenticalOutput, temperatureIdenticalArgMap);

		// known values for testing
		double liquidEnthalpyKnown[nTest] = {762.683e3,1610.15e3,289.229e3};
		double temperatureKnown[nTest] = {392.140,504.658,297.004};
		double volumeKnown[nTest] = {0.00105925,0.00119519,0.00100259};
		double conductivityKnown[nTest] = {0.684097,0.503998,0.610291};
		double viscosityKnown[nTest] = {0.000218794,7.72970e-5,0.000853838};

		// test property functions against known values
		for (size_t i=0; i<nTest; i++)
		{
			AMP::pout << "\nValue test " << i << ":\n=====================\n";
			// saturated liquid enthalpy
			if( !AMP::Utilities::approx_equal(liquidEnthalpyOutput[i], liquidEnthalpyKnown[i], 0.01) )
			{
				ut.failure("The answer is wrong.");
				AMP::pout << "The calculated saturated liquid enthalpy was " << liquidEnthalpyOutput[i] << " J/kg and actual is ";
				AMP::pout << liquidEnthalpyKnown[i] << " J/kg\n";
			}
			else AMP::pout << "saturated liquid enthalpy value is approximately equal to known value.\n";
			// temperature
			if( !AMP::Utilities::approx_equal(temperatureOutput[i], temperatureKnown[i], 0.01) )
			{
				ut.failure("The answer is wrong.");
				AMP::pout << "The calculated temperature was " << temperatureOutput[i] << " K and actual is ";
				AMP::pout << temperatureKnown[i] << " K\n";
			}
			else AMP::pout << "temperature value is approximately equal to known value.\n";
			// specific volume
			if( !AMP::Utilities::approx_equal(volumeOutput[i], volumeKnown[i], 0.01) )
			{
				ut.failure("The answer is wrong.");
				AMP::pout << "The calculated specific volume was " << volumeOutput[i] << " m3/kg and actual is ";
				AMP::pout << volumeKnown[i] << " m3/kg\n";
			}
			else AMP::pout << "specific volume value is approximately equal to known value.\n";
			// thermal conductivity
			if( !AMP::Utilities::approx_equal(conductivityOutput[i], conductivityKnown[i], 0.01) )
			{
				ut.failure("The answer is wrong.");
				AMP::pout << "The calculated thermal conductivity was " << conductivityOutput[i] << " W/m-K and actual is ";
				AMP::pout << conductivityKnown[i] << " W/m-K\n";
			}
			else AMP::pout << "thermal conductivity value is approximately equal to known value.\n";
			// dynamic viscosity
			if( !AMP::Utilities::approx_equal(viscosityOutput[i], viscosityKnown[i], 0.01) )
			{
				ut.failure("The answer is wrong.");
				AMP::pout << "The calculated dynamic viscosity was " << viscosityOutput[i] << " Pa-s and actual is ";
				AMP::pout << viscosityKnown[i] << " Pa-s\n";
			}
			else AMP::pout << "dynamic viscosity value is approximately equal to known value.\n";
		}

		// test temperature values against each other
		for (size_t i=0; i<nIdentical; i++)
		{
			if (!AMP::Utilities::approx_equal(temperatureIdenticalOutput[i], temperatureOutput_mat[i]))
			{
				AMP::pout << "Identical values temperature test 1 failed: 1st value: " << temperatureIdenticalOutput[i];
				AMP::pout << " and 2nd value: " << temperatureOutput_mat[i] << "\n";
				good = false;
			}
			if (!AMP::Utilities::approx_equal(temperatureIdenticalOutput[0], temperatureIdenticalOutput[i]))
			{
				AMP::pout << "Identical values temperature test 2 failed: 1st value: " << temperatureIdenticalOutput[0];
				AMP::pout << " and 2nd value: " << temperatureIdenticalOutput[i] << "\n";
				good = false;
			}
		}
	}

	// set defaults for temperature
	   std::vector<double> temperatureDefaults(2);
	   temperatureDefaults[0] = 100e3;	// enthalpy: 100 kJ/kg
	   temperatureDefaults[1] = 0.5e6;	// pressure: 0.5 MPa
	   temperatureProperty->set_defaults(temperatureDefaults);

	   // temperature test with one argument given: enthalpy
	   {
		   double knownSolution = 392.224;	// T [K]	@ {0.5 MPa, 500 kJ/kg}
		   std::map<std::string, boost::shared_ptr<std::vector<double> > > argMap;
		   argMap.insert( std::make_pair( "enthalpy", enthalpyIdenticalInput ) );
		   std::vector<double> temperatureOutput_def(temperatureIdenticalOutput);
		   temperatureProperty->evalv(temperatureOutput_def, argMap);
		for (size_t i=0; i<nIdentical; i++)
		{
			if (!AMP::Utilities::approx_equal(temperatureOutput_def[i], knownSolution, 0.01))
			{
				AMP::pout << "Temperature computed with one argument given is incorrect.\n";
				good = false;
			}
		}
	   }
	   // temperature test with no argument given
	   {
		   double knownSolution = 296.899;	// T [K]	@ {0.5 MPa, 100 kJ/kg}
		   std::map<std::string, boost::shared_ptr<std::vector<double> > > argMap;
		   std::vector<double> temperatureOutput_def(temperatureIdenticalOutput);
		   temperatureProperty->evalv(temperatureOutput_def, argMap);
		for (size_t i=0; i<nIdentical; i++)
		{
			if (!AMP::Utilities::approx_equal(temperatureOutput_def[i], knownSolution, 0.01))
			{
				AMP::pout << "Temperature computed with no argument given is incorrect.\n";
				good = false;
			}
		}
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
