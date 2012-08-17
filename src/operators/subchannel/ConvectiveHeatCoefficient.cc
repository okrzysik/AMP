#include "ConvctiveHeatCoefficient.h"
#include "utils/Utilities.h"
#include <cmath>


namespace AMP {
  namespace Operator {


    ConvctiveHeatCoefficient :: ConvctiveHeatCoefficient(const
        boost::shared_ptr<RobinPhysicsModelParameters>& params) : RobinPhysicsModel (params) {

      AMP_INSIST( (params->d_db->keyExists("Material")), "Convective Heat  Coefficient Key ''Material'' is missing!" );
      std::string matname = params->d_db->getString("Material");

      d_material = AMP::voodoo::Factory<AMP::Materials::Material>::instance().create(matname);

      AMP_INSIST( (params->d_db->keyExists("Property")), "Convective Heat Coefficient Key ''Property'' is missing!" );
      std::string propname = params->d_db->getString("Property");
      d_property = d_material->property(propname);

      d_defaults.resize(d_property->get_number_arguments());
      std::vector<std::vector<double> > ranges = d_property->get_arg_ranges();
      for (size_t i=0; i<d_defaults.size(); ++i) { d_defaults[i] = ranges[i][0]*(1.0000001);}
      if (params->d_db->keyExists("Defaults")) {
        // check for correct names
        boost::shared_ptr<Database> defaults_db = params->d_db->getDatabase("Defaults");
        std::vector<std::string> defaultkeys = defaults_db->getAllKeys();
        AMP_INSIST(defaultkeys.size() == d_property->get_number_arguments(),
            "Incorrect number of defaults supplied.");
        d_argNames = d_property-> get_arguments();
        for (std::vector<std::string>::iterator key= defaultkeys.begin(); key!=defaultkeys.end(); ++key)
        {
          std::vector<std::string>::iterator hit = std::find(argnames.begin(), argnames.end(), *key);
          AMP_INSIST(hit!=argnames.end(), std::string("Argument name ")+*key+std::string(" is invalid"));
        }

        // load defaults into the material property, checking range validity
        for (size_t i=0; i<argnames.size(); ++i) {
          d_defaults[i] = defaults_db->getDouble(argnames[i]);
          AMP_INSIST(d_property->in_range(argnames[i], d_defaults[i]),
              std::string("Default for argument ")+argnames[i]+std::string(" is out of range"));
        }
      }
      d_property->set_defaults(d_defaults);

    }

    void ConvctiveHeatCoefficient :: getConductance(std::vector<double> & beta,
        std::vector<double> & gamma, const std::vector<std::vector <double> > & inputVectors)
    {
      size_t numArgs = inputVectors.size();
      size_t argSize = inputVectors[0].size();
      
      std::map<std::string, boost::shared_ptr<std::vector<double> > > argMap;
      for (size_t i=0; i<argSize; i++) {
        argMap.insert(std::make_pair(d_argNames[i],new std::vector<double>(inputVectors[i].begin(), inputVectors[i].end())));
      }

      d_property->evalv(beta, args);
      d_property->evalv(gamma, args);

    }

  }
}
