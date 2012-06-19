/*
 * \file materials/WaterLibrary.h
 * \brief Material library that contains all the properties of water as a coolant for subchannel flow.
 */

#include "WaterLibrary.h"

#include "Property.h"
#include "Material.h"
#include "Utilities.h"

#include <string>
#include <valarray>

namespace AMP { 
namespace Materials {
namespace WaterLibrary_NS {

//=================== Constants =====================================================

	static const std::string name_base("WaterLibrary");
	static const std::string source("\
J.H. McFadden, et. al, RETRAN-02, Electric Power Research Institute, Technical Document NP-1850-CCMA, November 1984");

  // Temperature as a function of enthalpy and pressure from TLIQ in WATTPRO.h in COBRA-EN.
	static const unsigned int TempNumParams = 3;
	static const double       TempParams[3] = {1, 2, 3};
	static const unsigned int TempNumArgs   = 2;
	static const std::string  TempArgs[TempNumArgs] = {"enthalpy", "pressure"};
	static const double       TempHminVal = 0;  // I don't know.
	static const double       TempHmaxVal = 1e100; // I don't know.
	static const double       TempPminVal = 0; 
	static const double       TempPmaxVal = 906; // psi, need to convert to Pa.
	static const double       TempRanges[2][2]={{TempHminVal, TempHmaxVal}, {TempPminVal, TempPmaxVal}};

//=================== Classes =======================================================

	class TemperatureProp : public Property<double>{
	public:
		TemperatureProp() :
			Property<double> (	name_base + "_" + "Temperature",	// Name string
								source,										// Reference source
								TempParams,									// Property parameters
								TempNumParams, // Not sure why this isn't TempNumArgs; Number of parameters
								TempArgs,  								// Names of arguments to the eval function
								TempNumArgs,								// Number of arguments
								TempRanges ){}									// Range of variables

		virtual double eval( std::vector<double>& args );
	};

//=================== Functions =====================================================

	inline double TemperatureProp::eval( std::vector<double>& args ){
    double H            = args[0];  // local enthalpy in J/kg
    double P            = args[1];  // local pressure in Pa
    double T;                       // temperature in Kelvin
    
    AMP_ASSERT(H >= TempHminVal and H <= TempHmaxVal);
    AMP_ASSERT(P >= TempPminVal and P <= TempPmaxVal);
	   
	  std::valarray<double> Param = get_parameters();
		T = Param[0] + Param[1]*H + Param[2]*P; 
		return T;
	}

}

//=================== Materials =====================================================

WaterLibrary::WaterLibrary()
{
		d_propertyMap = new std::map<std::string, boost::shared_ptr<Property<double> > >();
		INSERT_PROPERTY_IN_MAP(Temperature,	WaterLibrary_NS);
}


} 
}


