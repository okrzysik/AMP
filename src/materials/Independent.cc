/*
 * Independent.h
 *
 *  Created on: Mar 11, 2010
 *	  Author: bm, gad
 */

#include "Independent.h"

#include "Property.h"
#include "Material.h"

#include <string>

namespace AMP { 
namespace Materials {

namespace Independent_NS {

//  =================== Constants =====================================================

	static const std::string name_base("Independent");
	static const std::string source("\
T=273, \
Bogdan Mihaila, Marius Stan, Juan Ramirez, \
Alek Zubelewicz, Petrica Cristea, \
Journal of Nuclear Materials 394 (2009) 182--189");

	/** \todo {the values of 1.0 below need to be filled in with defaults
	 *		 for T=273, U=0, B=0}
	 */

	static const double thermalval=1.;
	static const double fickval=1.;
	static const double soretval=1.;

	static const double densval=1.;
	static const double alphaval=1.;
	static const double heatcpval=1.;

	static const double youngsval=1.;
	static const double pratioval=0.290;

	static const double fickVectorVal[3]={1.,1.,1.};
	static const double fickTensorVal[3*3]={1.,1.,1.,1.,1.,1.,1.,1.,1.};

	static const std::string arguments[0];

//  =================== Classes =======================================================

	class ThermalConductivityProp : public Property<double>{
	public:
		ThermalConductivityProp() :
			Property<double> (	name_base + "_" + "ThermalConductivity",	// Name string
								source,										// Reference source
								&thermalval,								// Property parameters
								1U,											// Number of parameters
								arguments,									// Names of arguments to the eval function
								0 ){}										// Number of arguments

		virtual double eval( std::vector<double>& args );
	};

	class FickCoefficientProp : public Property<double>{
	public:
		FickCoefficientProp() :
			Property<double> (	name_base + "_" + "FickCoefficient",		// Name string
								source,										// Reference source
								&fickval,									// Property parameters
								1U,											// Number of parameters
								arguments,									// Names of arguments to the eval function
								0 ){}										// Number of arguments

		virtual double eval( std::vector<double>& args );
	};

	class SoretCoefficientProp : public Property<double>{
	public:
		SoretCoefficientProp() :
			Property<double> (	name_base + "_" + "FickCoefficient",		// Name string
								source,										// Reference source
								&fickval,									// Property parameters
								1U,											// Number of parameters
								arguments,									// Names of arguments to the eval function
								0 ){}										// Number of arguments

		virtual double eval( std::vector<double>& args );
	};

	class DensityProp : public Property<double>{
	public:
		DensityProp() :
			Property<double> (	name_base + "_" + "Density",				// Name string
								source,										// Reference source
								&densval,									// Property parameters
								1U,											// Number of parameters
								arguments,									// Names of arguments to the eval function
								0 ){}										// Number of arguments

		virtual double eval( std::vector<double>& args );
	};

	class ThermalExpansionProp : public Property<double>{
	public:
		ThermalExpansionProp() :
			Property<double> (	
								name_base + "_" + "ThermalExpansion",		// Name string
								source,										// Reference source
								&alphaval,									// Property parameters
								1U,											// Number of parameters
								arguments,									// Names of arguments to the eval function
								0 ){}										// Number of arguments

		virtual double eval( std::vector<double>& args );
	};

	class HeatCapacityPressureProp : public Property<double>{
	public:
		HeatCapacityPressureProp() :
			Property<double> (	
								name_base + "_" + "HeatCapacityPressure",	// Name string
								source,										// Reference source
								&heatcpval,									// Property parameters
								1U,											// Number of parameters
								arguments,									// Names of arguments to the eval function
								0 ){}										// Number of arguments
								
		virtual double eval( std::vector<double>& args );
	};

	class YoungsModulusProp : public Property<double>{
	public:
		YoungsModulusProp() :
			Property<double> (	
								name_base + "_" + "YoungsModulus",			// Name string
								source,										// Reference source
								&youngsval,									// Property parameters
								1U,											// Number of parameters
								arguments,									// Names of arguments to the eval function
								0 ){}										// Number of arguments

		virtual double eval( std::vector<double>& args );
	};

	class PoissonRatioProp : public Property<double>{
	public:
		PoissonRatioProp() :
			Property<double> (	name_base + "_" + "PoissonRatio",			// Name string
								source,										// Reference source
								&pratioval,									// Property parameters
								1U,											// Number of parameters
								arguments,									// Names of arguments to the eval function
								0 ){}										// Number of arguments

		virtual double eval( std::vector<double>& args );
	};

	class DTThermalConductivityProp : public Property<double>{
	public:
		DTThermalConductivityProp() :
			Property<double> (	name_base + "_" + "DTThermalConductivity",	// Name string
								source,										// Reference source
								&thermalval,								// Property parameters
								1U,											// Number of parameters
								arguments,									// Names of arguments to the eval function
								0 ){}										// Number of arguments

		virtual double eval( std::vector<double>& args );
	};

	class DTFickCoefficientProp : public Property<double>{
	public:
		DTFickCoefficientProp() :
			Property<double> (	name_base + "_" + "DTFickCoefficient",		// Name string
								source,										// Reference source
								&fickval,									// Property parameters
								1U,											// Number of parameters
								arguments,									// Names of arguments to the eval function
								0 ){}										// Number of arguments

		virtual double eval( std::vector<double>& args );
	};

	class DTSoretCoefficientProp : public Property<double>{
	public:
		DTSoretCoefficientProp() :
			Property<double> (	name_base + "_" + "DTSoretCoefficient",		// Name string
								source,										// Reference source
								&soretval,									// Property parameters
								1U,											// Number of parameters
								arguments,									// Names of arguments to the eval function
								0 ){}										// Number of arguments

		virtual double eval( std::vector<double>& args );
	};

	class DxThermalConductivityProp : public Property<double>{
	public:
		DxThermalConductivityProp() :
			Property<double> (	name_base + "_" + "DxThermalConductivity",	// Name string
								source,										// Reference source
								&thermalval,								// Property parameters
								1U,											// Number of parameters
								arguments,									// Names of arguments to the eval function
								0 ){}										// Number of arguments

		virtual double eval( std::vector<double>& args );
	};

	class DxFickCoefficientProp : public Property<double>{
	public:
		DxFickCoefficientProp() :
			Property<double> (	name_base + "_" + "DxFickCoefficient",		// Name string
								source,										// Reference source
								&fickval,									// Property parameters
								1U,											// Number of parameters
								arguments,									// Names of arguments to the eval function
								0 ){}										// Number of arguments

		virtual double eval( std::vector<double>& args );
	};

	class DxSoretCoefficientProp : public Property<double>{
	public:
		DxSoretCoefficientProp() :
			Property<double> (	name_base + "_" + "DxSoretCoefficient",		// Name string
								source,										// Reference source
								&soretval,									// Property parameters
								1U,											// Number of parameters
								arguments,									// Names of arguments to the eval function
								0 ){}										// Number of arguments

		virtual double eval( std::vector<double>& args );
	};

	class VectorFickCoefficientProp : public Property<double>{
	public:
		VectorFickCoefficientProp() :
			Property<double> (	name_base + "_" + "VectorFickCoefficient",	// Name string
								source,										// Reference source
								fickVectorVal,								// Property parameters
								3U,											// Number of parameters
								arguments,									// Names of arguments to the eval function
								0,											// Number of arguments
								(double(*)[2])(NULL), false, true, false){}

		virtual std::vector<double> evalVector( std::vector<double>& args );
	};

	class TensorFickCoefficientProp : public Property<double>{
	public:
		TensorFickCoefficientProp() :
			Property<double> (	name_base + "_" + "FickCoefficient",		// Name string
								source,										// Reference source
								fickTensorVal,								// Property parameters
								9U,											// Number of parameters
								arguments,									// Names of arguments to the eval function
								0,											// Number of arguments
								(double(*)[2])(NULL), false, false, true){}

		virtual std::vector<std::vector<double> > evalTensor( std::vector<double>& args );
	};

    static const unsigned int numberThDiffParams = 2;
    static double thermalDiffusionParams[numberThDiffParams] = {1.,1.};

	static const std::string thermDiffArgs[0];
	static const unsigned int numberThermDiffArgs = 0;
	static const  double thermDiffRanges[0][2]={};

#define THERMAL_DIFFUSION_DERIVATIVE
#include "ThermalDiffusionCoefficientProp.h"
#undef THERMAL_DIFFUSION_DERIVATIVE

//  =================== Functions =====================================================

	inline double ThermalConductivityProp::eval( std::vector<double>& args ){
		  return get_parameters()[0];
	}

	inline double FickCoefficientProp::eval( std::vector<double>& args ){
		  return get_parameters()[0];
	}

	inline double SoretCoefficientProp::eval( std::vector<double>& args ){
		  return get_parameters()[0];
	}

	inline double DTThermalConductivityProp::eval( std::vector<double>& args ){
		  return 0.;
	}

	inline double DxThermalConductivityProp::eval( std::vector<double>& args ){
		  return 0.;
	}

	inline double DTFickCoefficientProp::eval( std::vector<double>& args ){
		  return 0.;
	}

	inline double DxFickCoefficientProp::eval( std::vector<double>& args ){
		  return 0.;
	}

	inline double DTSoretCoefficientProp::eval( std::vector<double>& args ){
		  return 0.;
	}

	inline double DxSoretCoefficientProp::eval( std::vector<double>& args ){
		  return 0.;
	}

	inline double DensityProp::eval( std::vector<double>& args ){
		  return get_parameters()[0];
	}

	inline double ThermalExpansionProp::eval( std::vector<double>& args ){
		  return get_parameters()[0];
	}

	inline double HeatCapacityPressureProp::eval( std::vector<double>& args ){
		  return get_parameters()[0];
	}

	inline double YoungsModulusProp::eval( std::vector<double>& args ){
		  return get_parameters()[0];
	}

	inline double PoissonRatioProp::eval( std::vector<double>& args ){
		  return get_parameters()[0];
	}

	inline std::vector<double> VectorFickCoefficientProp::evalVector( std::vector<double>& args ){
		  std::vector<double> result(3);
		  for (size_t i=0; i<3; i++) result[i] = d_params[i];
		  return result;
	}

	inline std::vector<std::vector<double> > TensorFickCoefficientProp::evalTensor( std::vector<double>& args ){
		  std::vector<std::vector<double> > result(3, std::vector<double>(3));
		  for (size_t i=0; i<3; i++) for (size_t j=0; j<3; j++) result[i][j] = d_params[i*3+j];
		  return result;
	}
}

//  =================== Materials =====================================================

Independent::Independent()
{
		d_propertyMap = new std::map<std::string, boost::shared_ptr<Property<double> > >();
		INSERT_PROPERTY_IN_MAP(ThermalConductivity, 	Independent_NS);
		INSERT_PROPERTY_IN_MAP(FickCoefficient,	      	Independent_NS);
		INSERT_PROPERTY_IN_MAP(SoretCoefficient,	      	Independent_NS);
		INSERT_PROPERTY_IN_MAP(DTThermalConductivity, 	Independent_NS);
		INSERT_PROPERTY_IN_MAP(DTFickCoefficient,	      	Independent_NS);
		INSERT_PROPERTY_IN_MAP(DTSoretCoefficient,	      	Independent_NS);
		INSERT_PROPERTY_IN_MAP(DxThermalConductivity, 	Independent_NS);
		INSERT_PROPERTY_IN_MAP(DxFickCoefficient,	      	Independent_NS);
		INSERT_PROPERTY_IN_MAP(DxSoretCoefficient,	      	Independent_NS);
		INSERT_PROPERTY_IN_MAP(Density, 				Independent_NS);
		INSERT_PROPERTY_IN_MAP(HeatCapacityPressure, 	Independent_NS);
		INSERT_PROPERTY_IN_MAP(ThermalExpansion, 		Independent_NS);
		INSERT_PROPERTY_IN_MAP(YoungsModulus, 			Independent_NS);
		INSERT_PROPERTY_IN_MAP(PoissonRatio, 			Independent_NS);
		INSERT_PROPERTY_IN_MAP(ThermalDiffusionCoefficient,	Independent_NS);
}


} 
}


