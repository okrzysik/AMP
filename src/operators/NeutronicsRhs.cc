#include "NeutronicsRhs.h"

/*AMP Files */
#include "operators/Operator.h"
#include "NeutronicsRhs.h"
#include "NeutronicsRhsParameters.h"
#include "vectors/Vector.h"
//#include "ampmesh/MeshUtils.h"


#include "utils/InputDatabase.h"

/*Boost Files */
#include "boost/shared_ptr.hpp"

#include <vector>
#include <cmath>

namespace AMP {
namespace Operator {

  /*
  *************************************************************************
   * Constructor for NeutronicsRhs.  The constructor initializes the values   *
   * from the parameters.                                                  *
  *************************************************************************
  */
  NeutronicsRhs :: NeutronicsRhs(SP_Parameters parameters)
    : Operator(parameters),
      d_timeStep(0) {
      AMP_ASSERT(parameters);
      d_Mesh = parameters->d_Mesh;
      d_timeStepInSeconds = 0.;
      d_secondsPerDay = 86400.;      
      getFromInput(parameters->d_db);
      if( !d_useFixedValue ) { 
        int numValues; 
        numValues = parameters->d_db->getInteger("numValues");
        d_values.resize(d_numTimeSteps, Vec_Dbl(numValues, 0.0));
        char key[100];
        for( int step = 0; step < d_numTimeSteps; step++) {
          sprintf(key, "value_%d", step);
          AMP_INSIST( parameters->d_db->keyExists(key), "Key is missing!" );
          parameters->d_db->getDoubleArray(key, &d_values[step][0], numValues);
        }
      }

    }

  /* 
  *************************************************************************
   * Destructor.                                                           *
  *************************************************************************
  */
  NeutronicsRhs::~NeutronicsRhs()
  {
  }

  /*
  *************************************************************************
   * If simulation is not from restart, read data from input database.     *
   * Otherwise, override restart values for a subset of the data members   *
   * with those found in input.                                            *
  *************************************************************************
  */
  void 
    NeutronicsRhs::getFromInput(SP_Database db)
    {
AMP_ERROR("NeutronicsRhs is not converted yet");
/*
      AMP_ASSERT(db);

      // define the source type and create the output variable.
      std::string str = db->getStringWithDefault("type", "Power");
      d_type = str2id(str);

      std::string outVarName = db->getStringWithDefault("OutputVariable", str);
      d_outputVariable.reset ( new  HexGaussPointVariable (outVarName,d_Mesh) );

      // number of time steps
      d_numTimeSteps = db->getIntegerWithDefault("numTimeSteps", 1);
      AMP_ASSERT(d_numTimeSteps > 0);
      d_timeStepsInDays.resize(d_numTimeSteps);
      d_fixedValues.resize(d_numTimeSteps);

      if (db->keyExists("numTimeSteps")) {
        // time-step sizes
        if (db->keyExists("timeSteps")) {
          std::string tmp= "timeSteps";
          db->getDoubleArray(tmp, &d_timeStepsInDays[0], d_numTimeSteps);
        } else {
          // default value is only valid if the default number of time steps is used.
          AMP_ASSERT(d_numTimeSteps == 1);
          d_timeStepsInDays[0] = 100000000000.;
        }
        // default value is 1
      }

      // Power in Watts per gram
      d_useFixedValue = db->getBoolWithDefault("useFixedValue", 1);
      if( d_useFixedValue ) { 
        if (db->keyExists("fixedValues")) {
          std::string tmp= "fixedValues";
          db->getDoubleArray(tmp, &d_fixedValues[0], d_numTimeSteps);
        } else {
          // default value is only valid if the default number of time steps is used.
          AMP_ASSERT(d_numTimeSteps == 1);
          d_fixedValues[0] = 1.;
        }
      }
*/
    }

  /*
  *************************************************************************
   * Write out class version number and data members to database.          *
  *************************************************************************
  */
  void 
    NeutronicsRhs::putToDatabase(SP_Database db)
    {
      AMP_ASSERT(!db.use_count());
      db->putInteger("numTimeSteps", d_numTimeSteps);
      db->putDoubleArray("timeSteps", d_timeStepsInDays);
      db->putDoubleArray("fixedValues", d_fixedValues);
    }

  /*
  *************************************************************************
   * Print class data members to given output stream.                      *
  *************************************************************************
  */
  void 
    NeutronicsRhs::printClassData(std::ostream& os) const
    {
      os << "\nNeutronicsRhs::printClassData..." << std::endl;
      os << "d_numTimeSteps = " << d_numTimeSteps << std::endl;
    }

  /*
  *************************************************************************
   * Reset the class.                                                      *
  *************************************************************************
  */
  void 
    NeutronicsRhs :: reset(const SP_OperatorParameters & parameters) {

      AMP_ASSERT(parameters.get() != NULL);
      d_db = parameters->d_db;
      SP_Parameters params = boost::dynamic_pointer_cast<NeutronicsRhsParameters,
                                                         OperatorParameters         >(parameters);
      AMP_ASSERT(params.get() != NULL);
      AMP_ASSERT(((params->d_db).get())!=NULL);
      getFromInput(params->d_db);
      
      if( !d_useFixedValue ) { 
        int numValues; 
        numValues = params->d_db->getInteger("numValues");
        d_values.resize(d_numTimeSteps, Vec_Dbl(numValues, 0.0));
        char key[100];
        for( int step = 0; step < d_numTimeSteps; step++) {
          sprintf(key, "value_%d", step);
          AMP_INSIST( params->d_db->keyExists(key), "Key is missing!" );
          params->d_db->getDoubleArray(key, &d_values[step][0], numValues);
        }
      }

    }

  /*
  *************************************************************************
   * Provide the Apply function
  *************************************************************************
  */
  void 
    NeutronicsRhs :: apply(const  SP_Vector & f, 
                           const  SP_Vector & u, 
                                  SP_Vector & r,
                           const  double      a,
                           const  double      b) {
    (void) f; (void) u;
AMP_ERROR("NeutronicsRhs is not converted yet");
/*
      // NeutronicsRhs is made to provide a power, so a and b are not optional.
      AMP_ASSERT(AMP::Utilities::approx_equal(a,1.));
      AMP_ASSERT(AMP::Utilities::approx_equal(b,0.));

      SP_Vector rInternal = r->subsetVectorForVariable(d_outputVariable);

      AMP_ASSERT(rInternal!=NULL);
      
      // determine the present time
      int this_step = d_timeStep;

      // compute power distribution
      if(d_useFixedValue) {
        double value = d_fixedValues[this_step];
        rInternal->setToScalar(value);
      } else {
        AMP::Mesh::MeshManager::Adapter::ElementIterator  elem      = d_Mesh->beginElement();
        AMP::Mesh::MeshManager::Adapter::ElementIterator  end_elems = d_Mesh->endElement();

        int gp = 0;
        for( ; elem != end_elems; ++elem) {
          for( unsigned int i = 0; i < 8; gp++ , i++ ) {
            AMP::Mesh::DOFMap::shared_ptr  dof_map = d_Mesh->getDOFMap ( d_outputVariable );
            std::vector<unsigned int> ndx;
            std::vector<unsigned int> empty;
            dof_map->getDOFs ( *elem , ndx , empty );
            int  offset = ndx[i];
            rInternal->setValueByGlobalID ( offset, d_values[this_step][gp] );
          }//end for gauss-points
        }//end for elements
      }
*/
    }



  NeutronicsRhs::SourceType NeutronicsRhs::str2id(std::string str)
  {
    if      (str == "Power"        ) { return Power;         }
    else if (str == "Oxygen"       ) { return Oxygen;        }
    else if (str == "Metal"        ) { return Metal;         }
    else if (str == "FissionGas"   ) { return FissionGas;    }
    else if (str == "Isotopes"     ) { return Isotopes;    }
    else {
      std::string msg = "str2id could not find the right enumerated ID with string !" + str + "!.  Options are: Power, Oxygen, Metal, and FissionGas";
      AMP_INSIST(false,msg);
    } 
    return NUM_SOURCE_TYPES;
  }


      void NeutronicsRhs::setOutputVariableName(const std::string & name, int varId) {
        (void) varId;      
AMP_ERROR("NeutronicsRhs is not converted yet");
/*
        d_outputVariable->setName(name);
*/
      }

      AMP::LinearAlgebra::Variable::shared_ptr NeutronicsRhs::getOutputVariable() {
AMP_ERROR("NeutronicsRhs is not converted yet");
/*
        return d_outputVariable;
*/
      }


  /*SP_HexGaussPointVariable NeutronicsRhs::createOutputVariable (const std::string & name, int varId = -1) 
      {
        (void) varId;    
        SP_HexGaussPointVariable var( new HexGaussPointVariable (name) );
        return var;
      }*/


}
 }
