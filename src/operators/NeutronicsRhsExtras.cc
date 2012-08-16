
#include "NeutronicsRhsExtras.h"

/*AMP Files */
#include "operators/Operator.h"
#include "NeutronicsRhsExtras.h"
#include "NeutronicsRhsExtrasParameters.h"
#include "vectors/Vector.h"
#include "utils/Utilities.h"
#include "discretization/simpleDOF_Manager.h"

#include "utils/InputDatabase.h"

/*Boost Files */
#include "boost/shared_ptr.hpp"

#include <vector>
#include <cmath>

namespace AMP {
  namespace Operator {

    /*
     ***********************************************************************************
     * Constructor for NeutronicsRhsExtras.  The constructor initializes the values   *
     * from the parameters.                                                           *
     ***********************************************************************************
     */
    NeutronicsRhsExtras :: NeutronicsRhsExtras(SP_Parameters parameters)
      : Operator(parameters),
      d_timeStep(0),
      d_extrasId(0) 
    {
      AMP_ASSERT(parameters);
      d_Mesh = parameters->d_Mesh;
      d_numExtras = parameters->d_numExtras;
      d_extrasName = parameters->d_extrasName;
      d_timeStepInSeconds = 0.;
      d_secondsPerDay = 86400.;      
      getFromInput(parameters->d_db);
      if( !d_useFixedValue ) { 
        int numValues; 
        numValues = parameters->d_db->getInteger("numValues");
        Vec_Dbl1 tmp1(numValues,0);
        Vec_Dbl2 tmp2(d_numTimeSteps, tmp1);
        d_values.resize( d_numExtras, tmp2);

        for(int extras = 0 ; extras < d_numExtras ; extras++) {
          for(int t = 0 ; t < d_numTimeSteps ; t++) {
            char key[100];
            sprintf(key, "%s_%d", d_extrasName[extras].c_str(), t);
            parameters->d_db->getDoubleArray(key, &d_values[extras][t][0], numValues);
          }
        }

      }
    }


    /* 
     *************************************************************************
     * Destructor.                                                           *
     *************************************************************************
     */
    NeutronicsRhsExtras::~NeutronicsRhsExtras()
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
      NeutronicsRhsExtras::getFromInput(SP_Database db)
      {
        AMP_ASSERT(db);

        // define the source type and create the output variable.
        std::string str = db->getStringWithDefault("type", "Power");
        d_type = str2id(str);

        std::string outVarName = db->getStringWithDefault("OutputVariable", str);
        d_outputVariable.reset(new AMP::LinearAlgebra::Variable(outVarName)); 

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
      }

    /*
     *************************************************************************
     * Write out class version number and data members to database.          *
     *************************************************************************
     */
    void 
      NeutronicsRhsExtras::putToDatabase(SP_Database db)
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
      NeutronicsRhsExtras::printClassData(std::ostream& os) const
      {
        os << "\nNeutronicsRhsExtras::printClassData..." << std::endl;
        os << "d_numTimeSteps = " << d_numTimeSteps << std::endl;
      }

    /*
     *************************************************************************
     * Reset the class.                                                      *
     *************************************************************************
     */
    void 
      NeutronicsRhsExtras :: reset(const SP_OperatorParameters & parameters) {

        AMP_ASSERT(parameters.get() != NULL);
        d_db = parameters->d_db;
        SP_Parameters params = boost::dynamic_pointer_cast<NeutronicsRhsExtrasParameters,
                      OperatorParameters         >(parameters);
        AMP_ASSERT(params.get() != NULL);
        AMP_ASSERT(((params->d_db).get())!=NULL);
        getFromInput(params->d_db);

        if( !d_useFixedValue ) { 
          int numValues; 
          numValues = params->d_db->getInteger("numValues");
          d_values.resize( d_numExtras, Vec_Dbl2( d_numTimeSteps, Vec_Dbl1(numValues,0) ) );

          for(int extras = 0 ; extras < d_numExtras ; extras++) {
            for(int t = 0 ; t < d_numTimeSteps ; t++) {
              char key[100];
              sprintf(key, "%s_%d", d_extrasName[extras].c_str(), t);
              params->d_db->getDoubleArray(key, &d_values[extras][t][0], numValues);
            }
          }

        }

      }

    /*
     *************************************************************************
     * Provide the Apply function
     *************************************************************************
     */
    void 
      NeutronicsRhsExtras :: apply( AMP::LinearAlgebra::Vector::const_shared_ptr f, 
          AMP::LinearAlgebra::Vector::const_shared_ptr u, 
          AMP::LinearAlgebra::Vector::shared_ptr r,
          const  double      a,
          const  double      b) {
        (void) f; (void) u;
        // NeutronicsRhsExtras is made to provide a power, so a and b are not optional.
        AMP_ASSERT(AMP::Utilities::approx_equal(a,1.));
        AMP_ASSERT(AMP::Utilities::approx_equal(b,0.));

        AMP::LinearAlgebra::Vector::shared_ptr rInternal = r->subsetVectorForVariable(d_outputVariable);

        // determine the present time and extras
        int this_step = d_timeStep;
        int this_extrasId = d_extrasId;
        using namespace std;

        // compute power distribution
        if(d_useFixedValue) {
          double value = d_fixedValues[this_step];
          rInternal->setToScalar(value);
        } else {
          AMP::Mesh::MeshIterator  elem      = d_Mesh->getIterator(AMP::Mesh::Volume, 1);
          AMP::Mesh::MeshIterator  end_elems = elem.begin();

          int DOFsPerElement = 8;
          int ghostWidth = 1;
          bool split = true;
          AMP::Discretization::DOFManager::shared_ptr dof_map = AMP::Discretization::simpleDOFManager::create(d_Mesh, AMP::Mesh::Volume, ghostWidth, DOFsPerElement, split);

          int gp = 0;
          for( ; elem != end_elems; ++elem) {
            std::vector<size_t> gid;
            dof_map->getDOFs ( elem->globalID() , gid);
            for( unsigned int i = 0; i < gid.size(); gp++ , i++ ) {
              rInternal->setValueByGlobalID ( gid[i], d_values[this_extrasId][this_step][gp] );
            }//end for gauss-points
          }//end for elements

        }
      }


    NeutronicsRhsExtras::SourceType NeutronicsRhsExtras::str2id(std::string str)
    {
      if (str == "Isotopes"     ) { return Isotopes;    }
      else if (str == "Elements"     ) { return Elements;    }
      else {
        std::string msg = "str2id could not find the right enumerated ID with string !" + str + "!.  Options are: Isotopes and Elements";
        AMP_INSIST(false,msg);
      } 
      return NUM_SOURCE_TYPES;
    }

  }
}


