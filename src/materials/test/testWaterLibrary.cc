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
	   size_t n=10; // size of input and output arrays
	   boost::shared_ptr<std::vector<double> > enthalpyInput(new std::vector<double>(n));		// enthalpy input
	   boost::shared_ptr<std::vector<double> > pressureInput(new std::vector<double>(n));		// pressure input
	   boost::shared_ptr<std::vector<double> > temperatureInput(new std::vector<double>(n));	// temperature input
	   boost::shared_ptr<std::vector<double> > densityInput(new std::vector<double>(n));		// density input
	   vector<double> liquidEnthalpyOutput(n);	// saturated liquid enthalpy output
	   vector<double> temperatureOutput(n);		// temperature output
	   vector<double> volumeOutput(n);		// specific volume output
	   vector<double> conductivityOutput(n);	// thermal conductivity output
	   vector<double> viscosityOutput(n);		// dynamic viscosity output
	   for (size_t i=0; i<n; i++) 
	   {
		(*enthalpyInput)[i]=500e3;	// enthalpy: 500 kJ/kg 
		(*pressureInput)[i]=1e6;	// pressure: 1 MPa
		(*temperatureInput)[i]=400;	// temperature: 400 K
		(*densityInput)[i]=937.871;	// density: 937.871 kg/m3
	   }
	   
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
	
		// evaluate saturated liquid enthalpy function
		   std::vector<double> liquidEnthalpyOutput_mat(liquidEnthalpyOutput);
		   mat->property("SaturatedLiquidEnthalpy")->evalv(liquidEnthalpyOutput_mat, liquidEnthalpyArgMap);
		   liquidEnthalpyProperty->evalv(liquidEnthalpyOutput, liquidEnthalpyArgMap);
		// evaluate temperature function
		   std::vector<double> temperatureOutput_mat(temperatureOutput);
		   mat->property("Temperature")->evalv(temperatureOutput_mat, temperatureArgMap);
		   temperatureProperty->evalv(temperatureOutput, temperatureArgMap);
		// evaluate specific volume function
		   std::vector<double> volumeOutput_mat(volumeOutput);
		   mat->property("SpecificVolume")->evalv(volumeOutput_mat, volumeArgMap);
		   volumeProperty->evalv(volumeOutput, volumeArgMap);
		// evaluate thermal conductivity function
		   std::vector<double> conductivityOutput_mat(conductivityOutput);
		   mat->property("ThermalConductivity")->evalv(conductivityOutput_mat, conductivityArgMap);
		   conductivityProperty->evalv(conductivityOutput, conductivityArgMap);
		// evaluate dynamic viscosity function
		   std::vector<double> viscosityOutput_mat(viscosityOutput);
		   mat->property("DynamicViscosity")->evalv(viscosityOutput_mat, viscosityArgMap);
		   viscosityProperty->evalv(viscosityOutput, viscosityArgMap);

		double liquidEnthalpyKnown = 762.683e3;	// Hf [J/kg]	@ {1 MPa}
		double temperatureKnown = 392.140;	// T [K]	@ {1 MPa, 500 kJ/kg}
		double volumeKnown = 0.00105925;	// V [m3/kg]	@ {1 MPa, 500 kJ/kg}
		double conductivityKnown = 0.684097;	// k [W/m-K]	@ {400 K, 937.871 m3/kg}; @ {400 K, 1 MPa} 
		double viscosityKnown = 0.000218794;	// u [Pa-s]	@ {400 K, 937.871 m3/kg}; @ {400 K, 1 MPa} 

		// test saturated liquid enthalpy function against known value
		   if( !AMP::Utilities::approx_equal(liquidEnthalpyOutput[0], liquidEnthalpyKnown, 0.01) )
		{
			ut.failure("The answer is wrong.");
			AMP::pout << "The calculated saturated liquid enthalpy was " << liquidEnthalpyOutput[0] << " J/kg ";
			AMP::pout << "and actual is " << liquidEnthalpyKnown << " J/kg\n";
			
		}
		else AMP::pout << "saturated liquid enthalpy value is approximately equal to known value.\n";
		// test temperature function against known value
		   if( !AMP::Utilities::approx_equal(temperatureOutput[0], temperatureKnown, 0.01) )
		{
			ut.failure("The answer is wrong.");
			AMP::pout << "The calculated temperature was " << temperatureOutput[0] << " K and actual is " << temperatureKnown << " K\n";
		}
		else AMP::pout << "temperature value is approximately equal to known value.\n";
		// test specific volume function against known value
		   if( !AMP::Utilities::approx_equal(volumeOutput[0], volumeKnown, 0.01) )
		{
			ut.failure("The answer is wrong.");
			AMP::pout << "The calculated specific volume was " << volumeOutput[0] << " m3/kg and actual is " << volumeKnown << " m3/kg\n";
		}
		else AMP::pout << "specific volume value is approximately equal to known value.\n";
		// test thermal conductivity function against known value
		   if( !AMP::Utilities::approx_equal(conductivityOutput[0], conductivityKnown, 0.01) )
		{
			ut.failure("The answer is wrong.");
			AMP::pout << "The calculated thermal conductivity was " << conductivityOutput[0] << " W/m-K and actual is " << conductivityKnown << " W/m-K\n";
		}
		else AMP::pout << "thermal conductivity value is approximately equal to known value.\n";
		   if( !AMP::Utilities::approx_equal(viscosityOutput[0], viscosityKnown, 0.01) )
		{
			ut.failure("The answer is wrong.");
			AMP::pout << "The calculated dynamic viscosity was " << viscosityOutput[0] << " Pa-s and actual is " << viscosityKnown << " Pa-s\n";
		}
		else AMP::pout << "dynamic viscosity value is approximately equal to known value.\n";

		// test saturated liquid enthalpy values against each other
		   for (size_t i=0; i<n; i++) {good = good and AMP::Utilities::approx_equal(liquidEnthalpyOutput[i], liquidEnthalpyOutput_mat[i]);}
		   for (size_t i=0; i<n; i++) {good = good and AMP::Utilities::approx_equal(liquidEnthalpyOutput[0], liquidEnthalpyOutput[i]);}
		// test temperature values against each other
		   for (size_t i=0; i<n; i++) {good = good and AMP::Utilities::approx_equal(temperatureOutput[i], temperatureOutput_mat[i]);}
		   for (size_t i=0; i<n; i++) {good = good and AMP::Utilities::approx_equal(temperatureOutput[0], temperatureOutput[i]);}
	   }

	// set defaults for saturated liquid enthalpy
	   std::vector<double> liquidEnthalpyDefaults(1);
	   liquidEnthalpyDefaults[0] = 0.5e6;	// pressure: 0.5 MPa
	   liquidEnthalpyProperty->set_defaults(liquidEnthalpyDefaults);

	// set defaults for temperature
	   std::vector<double> temperatureDefaults(2);
	   temperatureDefaults[0] = 100e3;	// enthalpy: 100 kJ/kg
	   temperatureDefaults[1] = 0.5e6;	// pressure: 0.5 MPa
	   temperatureProperty->set_defaults(temperatureDefaults);

	   // saturated liquid enthalpy test with no argument given
	   {
		   double knownSolution = 640.185e3;	// Hf [J/kg]	@ {0.5 MPa}
		   std::map<std::string, boost::shared_ptr<std::vector<double> > > argMap;
		   argMap.insert( std::make_pair( "enthalpy", enthalpyInput ) );
		   std::vector<double> liquidEnthalpyOutput_def(liquidEnthalpyOutput);
		   liquidEnthalpyProperty->evalv(liquidEnthalpyOutput_def, argMap);
		for (size_t i=0; i<n; i++)
		{
			if (!AMP::Utilities::approx_equal(liquidEnthalpyOutput_def[i], knownSolution, 0.01))
			{
				AMP::pout << "Saturated liquid enthalpy computed with no argument given is incorrect.\n";
				good = false;
			}
		}
	   }
	   
	   // temperature test with one argument given: enthalpy
	   {
		   double knownSolution = 392.224;	// T [K]	@ {0.5 MPa, 500 kJ/kg}
		   std::map<std::string, boost::shared_ptr<std::vector<double> > > argMap;
		   argMap.insert( std::make_pair( "enthalpy", enthalpyInput ) );
		   std::vector<double> temperatureOutput_def(temperatureOutput);
		   temperatureProperty->evalv(temperatureOutput_def, argMap);
		for (size_t i=0; i<n; i++)
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
		   std::vector<double> temperatureOutput_def(temperatureOutput);
		   temperatureProperty->evalv(temperatureOutput_def, argMap);
		for (size_t i=0; i<n; i++)
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
