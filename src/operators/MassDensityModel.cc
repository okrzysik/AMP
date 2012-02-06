/*
 * MassDensityModel.cc
 *
 *  Created on: Aug 4, 2010
 *      Author: gad
 */

#include "MassDensityModel.h"
#include "materials/TensorProperty.h"

namespace AMP {
namespace Operator {

MassDensityModel::MassDensityModel
(const boost::shared_ptr<MassDensityModelParameters>& params) :
    ElementPhysicsModel(params)
{
    AMP_INSIST((params->d_db->keyExists("Material")),
            "Mass Key ''Material'' is missing!");
    std::string matname = params->d_db->getStringWithDefault("Material",
            "Independent");
    d_material = AMP::voodoo::Factory<AMP::Materials::Material>::instance().create(matname);

    if (params->d_db->keyExists("Property")) {
    	d_PropertyName = params->d_db->getString("Property");
    	if (params->d_db->keyExists("Parameters")) {
    		d_Parameters = params->d_db->getDoubleArray("Parameters");
    	}
    } else {
    	d_PropertyName = "unspecified";
    }

    AMP_INSIST((params->d_db->keyExists("Equation")),
            "Mass Key ''Equation'' is missing!");
    std::string eqnname = params->d_db->getString("Equation");
    if (eqnname == "Mechanics")
        d_equation = Mechanics;
        // The mechanics mass matrix is multiplied by the density of the material.
    else if (eqnname == "ThermalSource")
        d_equation = Mechanics;
        // Because the specific power (Watts/gram) are defined from the NeutronicsSource,
        // The mass matrix for the right-hand-side of the thermal equation must include
        // the density (grams/cubic-centimeter) to the get units correct (Watts/cc)
    else if (eqnname == "Thermal")
        d_equation = Thermal;
    else if (eqnname == "Chemical")
        d_equation = Chemical;
    else if (eqnname == "ManufacturedSource")
        d_equation = Manufactured;
        // used for manufactured solution testing rhs
    else
        AMP_INSIST(false, "Mass Equation name is invalid");

    d_UseBilogScaling = params->d_db->getBoolWithDefault("UseBilogScaling", false);
    if (d_UseBilogScaling)
    {
        AMP_INSIST(params->d_db->keyExists("BilogVariable"),
                "must specify BilogVariable");
        d_BilogVariable = params->d_db->getStringWithDefault("BilogVariable",
                "NONE");

        if (d_equation == Thermal) {
            d_BilogRange = d_material->property("ThermalConductivity")->get_arg_range("temperature");
        } else if (d_equation == Chemical) {
            d_BilogRange = d_material->property("FickCoefficient")->get_arg_range("concentration");
        }
        AMP_INSIST(d_BilogRange[1] > d_BilogRange[0],
                "material argument upper bound == lower bound");

        std::vector<std::string> names;
        if (d_equation == Thermal) {
            names = d_material->property("ThermalConductivity")->get_arguments();
        } else if (d_equation == Chemical) {
            names = d_material->property("FickCoefficient")->get_arguments();
        }
        d_BilogIndex = 999999;
        for (size_t i = 0; i < names.size(); i++)
        {
            if (names[i] == d_BilogVariable)
            {
                d_BilogIndex = i;
                break;
            }
        }
        AMP_INSIST(d_BilogIndex < 999999, "Did not find " + d_BilogVariable
                + " in list of material argument names");

        if (eqnname == "Thermal")
        {
            AMP_INSIST(d_BilogVariable == "temperature",
                    "thermal equation requires bilog scaling of temperature");
        }
        if (eqnname == "Chemical")
        {
            AMP_INSIST(d_BilogVariable == "concentration",
                    "chemical equation requires bilog scaling of concentration");
        }
    }

    if (d_equation == Manufactured)
    {
        AMP_INSIST(params->d_db->keyExists("ManufacturedSourceEquation"),
                "ManufacturedSourceEquation is missing");
        std::string mfgeqn = params->d_db->getString(
                "ManufacturedSourceEquation");

        if (mfgeqn == "Thermal")
            d_ManufacturedEquation = ThermalSrc;
        else if (mfgeqn == "Fick")
            d_ManufacturedEquation = FickSrc;
        else if (mfgeqn == "Soret")
            d_ManufacturedEquation = SoretSrc;
        else if (mfgeqn == "FickSoret")
            d_ManufacturedEquation = FickSoretSrc;
        else
            AMP_INSIST(false, "invalid value for ManufacturedSourceEquation");

        AMP_INSIST(params->d_db->keyExists("ManufacturedVariable"),"must specify ManufacturedVariable");
        d_ManufacturedVariable = params->d_db->getString("ManufacturedVariable");
        AMP_INSIST(d_ManufacturedVariable == "Temperature" or d_ManufacturedVariable == "Concentration",
                "ManufacturedVariable must have the values Temperature or Concentration");
        if (d_ManufacturedVariable == "Temperature") d_ManufacturedUseTemp = true; else d_ManufacturedUseTemp = false;
        if (d_ManufacturedVariable == "Concentration") d_ManufacturedUseConc = true; else d_ManufacturedUseConc = false;

        boost::shared_ptr<Database> mfg_db = params->d_db->getDatabase("ManufacturedSolution");
        d_ManufacturedSolution.reset(new ManufacturedSolution(mfg_db));
    }

    if (d_equation == Mechanics)
    {
      AMP::Materials::PropertyPtr property = d_material->property("Density");

      // load and check defaults
      // initially set them to the minimum of the range plus a bit
      std::vector<double> defaults(property->get_number_arguments());
      std::vector<std::vector<double> > ranges = property->get_arg_ranges();
      for (size_t i=0; i<defaults.size(); ++i) { defaults[i] = ranges[i][0]*(1.0000001);}
      if (params->d_db->keyExists("Defaults")) {
        // check for correct names
        boost::shared_ptr<Database> defaults_db = params->d_db->getDatabase("Defaults");
        std::vector<std::string> defaultkeys = defaults_db->getAllKeys();
        AMP_INSIST(defaultkeys.size() == property->get_number_arguments(),
            "Incorrect number of defaults supplied.");
        std::vector<std::string> argnames = property-> get_arguments();
        for (std::vector<std::string>::iterator key= defaultkeys.begin(); key!=defaultkeys.end(); ++key)
        {
          std::vector<std::string>::iterator hit = std::find(argnames.begin(), argnames.end(), *key);
          AMP_INSIST(hit!=argnames.end(), std::string("Argument name ")+*key+std::string(" is invalid"));
        }

        // load defaults into the material property, checking range validity
        for (size_t i=0; i<argnames.size(); ++i) {
          defaults[i] = defaults_db->getDouble(argnames[i]);
          AMP_INSIST(property->in_range(argnames[i], defaults[i]),
          std::string("Default for argument ")+argnames[i]+std::string(" is out of range"));
        }
      property->set_defaults(defaults);
      }
    }
}

void MassDensityModel::getDensityMechanics(std::vector<double> & result,
        const std::vector<double>& T, const std::vector<double>& U,
        const std::vector<double>& B)
{
    AMP_ASSERT((T.size() == U.size()) && (U.size() == result.size()) && (B.size() == U.size()));
    std::map<std::string, boost::shared_ptr<std::vector<double> > > inputMaterialParameters;

    std::string temperatureString = "temperature"; // in the future get from input file
    std::string burnupString = "burnup"; // in the future get from input file
    std::string oxygenString = "concentration"; // in the future get from input file

    boost::shared_ptr<std::vector<double> > tempVec(new std::vector<double>(T));
    boost::shared_ptr<std::vector<double> > burnupVec(new std::vector<double>(B));
    boost::shared_ptr<std::vector<double> > oxygenVec(new std::vector<double>(U));

    inputMaterialParameters.insert( std::make_pair( temperatureString, tempVec) );
    inputMaterialParameters.insert( std::make_pair( burnupString, burnupVec) );
    inputMaterialParameters.insert( std::make_pair( oxygenString, oxygenVec) );

    std::string denString = "Density";

    d_material->property(denString)->evalv(result, inputMaterialParameters);
}

void MassDensityModel::getDensityThermal(std::vector<double> & result,
        const std::vector<double>& T, const std::vector<double>& U,
        const std::vector<double>& B)
{
    AMP_ASSERT((T.size() == U.size()) && (U.size() == result.size()) && (B.size() == U.size()));
    unsigned int n = result.size();
    std::vector<double> density(n), specificheat(n);
    std::map<std::string, boost::shared_ptr<std::vector<double> > > args;
    args.insert(std::make_pair("temperature", boost::shared_ptr<std::vector<double> >(new std::vector<double>(T))));
    args.insert(std::make_pair("concentration", boost::shared_ptr<std::vector<double> >(new std::vector<double>(U))));
    args.insert(std::make_pair("burnup", boost::shared_ptr<std::vector<double> >(new std::vector<double>(B))));
    d_material->property("Density")->evalv(density, args);
    d_material->property("HeatCapacityPressure")->evalv(specificheat, args);
    for (unsigned int i = 0; i < n; i++)
        result[i] = density[i] * specificheat[i];

    if (d_UseBilogScaling)
    {
        DiffusionTransportModel::bilogScale(result, d_BilogRange[0],
                d_BilogRange[1]);
    }
}

void MassDensityModel::getDensityChemical(std::vector<double> & result,
        const std::vector<double>& T, const std::vector<double>& U,
        const std::vector<double>& B)
{
    AMP_ASSERT((T.size() == U.size()) && (U.size() == result.size()) && (B.size() == U.size()));

    for (size_t i = 0; i < result.size(); i++)
        result[i] = 1.;

    if (d_UseBilogScaling)
    {
        DiffusionTransportModel::bilogScale(result, d_BilogRange[0],
                d_BilogRange[1]);
    }
}

void MassDensityModel::getDensityManufactured(std::vector<double> & result,
        const std::vector<double>& T, const std::vector<double>& U,
        const std::vector<double>& B, const std::vector<Point>& xyz)
{

    AMP_ASSERT((T.size() == U.size()) && (U.size() == result.size()) && (B.size() == U.size()));
    AMP_ASSERT(xyz.size() == result.size());

    std::valarray<double> soln(10);
    size_t neval = result.size();

    AMP::Materials::PropertyPtr sourceProp;
    AMP::Materials::PropertyPtr dSourceProp;
    bool needD = false;

    if (d_PropertyName == "unspecified") {
		if (d_ManufacturedEquation == ThermalSrc) {
			sourceProp = d_material->property("ThermalConductivity");
			if (d_ManufacturedUseConc) {dSourceProp = d_material->property("DxThermalConductivity"); needD=true;}
			if (d_ManufacturedUseTemp) {dSourceProp = d_material->property("DTThermalConductivity"); needD=true;}
		} else if (d_ManufacturedEquation == FickSrc) {
			sourceProp = d_material->property("FickCoefficient");
			if (d_ManufacturedUseConc) {dSourceProp = d_material->property("DxFickCoefficient"); needD=true;}
			if (d_ManufacturedUseTemp) {dSourceProp = d_material->property("DTFickCoefficient"); needD=true;}
		} else if (d_ManufacturedEquation == SoretSrc) {
			sourceProp = d_material->property("ThermalDiffusionCoefficient");
			if (d_ManufacturedUseConc) {dSourceProp = d_material->property("DxThermalDiffusionCoefficient"); needD=true;}
			if (d_ManufacturedUseTemp) {dSourceProp = d_material->property("DTThermalDiffusionCoefficient"); needD=true;}
		} else if (d_ManufacturedEquation == FickSoretSrc) {
			AMP_INSIST(false, "cannot do Fick-Soret yet");
		}
    } else {
    	sourceProp = d_material->property(d_PropertyName);
    	if (d_Parameters.size()>0 and sourceProp->variable_number_parameters()) {
    		sourceProp->set_parameters_and_number(&d_Parameters[0], d_Parameters.size());
    	} else if (d_Parameters.size()>0) {
    		sourceProp->set_parameters(&d_Parameters[0], d_Parameters.size());
    	}
    }

    std::map<std::string, boost::shared_ptr<std::vector<double> > > args;
    args.insert(std::make_pair("temperature", boost::shared_ptr<std::vector<double> >(new std::vector<double>(T))));
    args.insert(std::make_pair("concentration", boost::shared_ptr<std::vector<double> >(new std::vector<double>(U))));
    args.insert(std::make_pair("burnup", boost::shared_ptr<std::vector<double> >(new std::vector<double>(B))));

    if (sourceProp->isScalar()) {
        std::vector<double> coeff(neval), dCoeff(neval,0.);
		sourceProp->evalv(coeff, args);
		if (needD) {
			dSourceProp->evalv(dCoeff, args);
		}

		for (size_t i = 0; i < neval; i++)
		{
			d_ManufacturedSolution->evaluate(soln, xyz[i](0), xyz[i](1), xyz[i](2));

			result[i] = coeff[i]*(soln[4]+soln[7]+soln[9]) +
						dCoeff[i]*(soln[1]*soln[1]+soln[2]*soln[2]+soln[3]*soln[3]);
		}
    }

    std::string propname = sourceProp->get_name();
    std::string solnname = d_ManufacturedSolution->get_name();
    bool isCylindrical = propname.find("CylindricallySymmetric") < propname.size() and
    		solnname.find("Cylindrical") < solnname.size();

    if (sourceProp->isTensor() and not isCylindrical) {
    	boost::shared_ptr<Materials::TensorProperty<double> > sourceTensorProp =
    			boost::dynamic_pointer_cast<Materials::TensorProperty<double> >(sourceProp);
    	std::vector<size_t> dimensions = sourceTensorProp->get_dimensions();
    	std::vector<std::vector<boost::shared_ptr<std::vector<double> > > >
    		coeff(dimensions[0], std::vector<boost::shared_ptr<std::vector<double> > >(dimensions[1],
    				boost::shared_ptr<std::vector<double> >(new std::vector<double>(neval))));
			sourceTensorProp->evalv(coeff, args);

		// 4 + xx xy xz yy yz zz =
		//      4  5  6  7  8  9 =
		// xx xy xz
		// yx yy yz
		// zx zy zz
		size_t xlate[3][3]={{4, 5, 6}, {5, 7, 8}, {6, 8, 9}};

		for (size_t k = 0; k < neval; k++)
		{
			d_ManufacturedSolution->evaluate(soln, xyz[k](0), xyz[k](1), xyz[k](2));

			result[k] = 0.;
			for (size_t i=0; i<dimensions[0]; i++) for (size_t j=0; j<dimensions[1]; j++) {
				result[k] +=  (*coeff[i][j])[k] * soln[xlate[i][j]];
			}
		}

    } else if (sourceProp->isTensor() and isCylindrical) {
    	boost::shared_ptr<Materials::TensorProperty<double> > sourceTensorProp =
    			boost::dynamic_pointer_cast<Materials::TensorProperty<double> >(sourceProp);
    	std::vector<size_t> dimensions = sourceTensorProp->get_dimensions();
    	AMP_ASSERT(dimensions[0] == 3 and dimensions[1] == 3);
    	std::vector<std::vector<boost::shared_ptr<std::vector<double> > > >
    		coeff(dimensions[0], std::vector<boost::shared_ptr<std::vector<double> > >(dimensions[1],
    				boost::shared_ptr<std::vector<double> >(new std::vector<double>(neval))));
    	std::vector<std::vector<boost::shared_ptr<std::vector<double> > > >
    		coeffr(dimensions[0], std::vector<boost::shared_ptr<std::vector<double> > >(dimensions[1],
    				boost::shared_ptr<std::vector<double> >(new std::vector<double>(neval))));
    	std::vector<std::vector<boost::shared_ptr<std::vector<double> > > >
    		coeffz(dimensions[0], std::vector<boost::shared_ptr<std::vector<double> > >(dimensions[1],
    				boost::shared_ptr<std::vector<double> >(new std::vector<double>(neval))));

        args.insert(std::make_pair("radius", boost::shared_ptr<std::vector<double> >(new std::vector<double>(neval))));
        args.insert(std::make_pair("theta", boost::shared_ptr<std::vector<double> >(new std::vector<double>(neval))));
        args.insert(std::make_pair("zee", boost::shared_ptr<std::vector<double> >(new std::vector<double>(neval))));
        std::vector<double> &radius = (*args.find("radius")->second);
        std::vector<double> &theta = (*args.find("theta")->second);
        std::vector<double> &zee = (*args.find("zee")->second);

		double Pi=3.1415926535898;
		for (size_t k = 0; k < neval; k++)
		{
			double x=xyz[k](0), y=xyz[k](1), z=xyz[k](2);
			double r=sqrt(x*x+y*y), th=acos(x/r); if (y<0) th = 2*Pi-th;
			radius[k] = r;
			theta[k]  = th;
			zee[k]    = z;
		}

		sourceTensorProp->setAuxiliaryData("derivative", 0);
		sourceTensorProp->evalv(coeff, args);

		sourceTensorProp->setAuxiliaryData("derivative", 1);
		sourceTensorProp->evalv(coeffr, args);

		sourceTensorProp->setAuxiliaryData("derivative", 2);
		sourceTensorProp->evalv(coeffz, args);

		// compute div (K . grad u) = div K . grad u + K : grad grad u
		for (size_t k = 0; k < neval; k++)
		{
			double x=xyz[k](0), y=xyz[k](1), z=xyz[k](2);
			double r=sqrt(x*x+y*y), th=acos(x/r); if (y<0) th = 2*Pi-th;
			// soln is the set of all derivatives wrto r, th, and z of order <= 2
			d_ManufacturedSolution->evaluate(soln, r, th, z);

			// convert hessian wrto r, th, z to hessian wrto x, y, z
			double r3 = r*r*r;
			std::vector<std::vector<double> > hessian(3, std::vector<double>(3, 0.));
			if (r>0) {
				hessian[0][0] = (y*y * soln[1] + r * x*x * soln[4])/r3;
				hessian[0][1] = x*y*(-soln[1]+r*soln[4])/r3;
				hessian[0][2] = x*soln[5]/r;
				hessian[1][0] = x*y*(-soln[1]+r*soln[4])/r3;
				hessian[1][1] = (x*x*soln[1] + r*y*y*soln[4])/r3;
				hessian[1][2] = y*soln[5]/r;
				hessian[2][0] = x*soln[5]/r;
				hessian[2][1] = y*soln[5]/r;
			}
			hessian[2][2] = soln[9];

			// compute k : grad grad term
			result[k] = 0.;
			for (size_t i=0; i<3; i++) for (size_t j=0; j<3; j++) {
				result[k] +=  (*coeff[i][j])[k] * hessian[i][j];
			}

			// compute grad k . grad term
			std::vector<double> factor(3,0.);
			std::vector<double> divK(3,0.);
			if (r>0) {
				factor[0] = x/r;
				factor[1] = y/r;
			}
			factor[2] = 1.;
			divK[0] = factor[0]*(*coeffr[0][0])[k] + factor[1]*(*coeffr[1][0])[k] + (*coeffz[2][0])[k];
			divK[1] = factor[0]*(*coeffr[0][1])[k] + factor[1]*(*coeffr[1][1])[k] + (*coeffz[2][1])[k];
			divK[2] = factor[0]*(*coeffr[0][2])[k] + factor[1]*(*coeffr[1][2])[k] + (*coeffz[2][2])[k];
			result[k] += divK[0]*soln[1] + divK[1]*soln[2] + divK[2]*soln[3];
		}
    }

}

}
}

