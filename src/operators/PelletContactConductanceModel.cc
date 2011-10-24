#include "PelletContactConductanceModel.h"
#include "utils/Utilities.h"
#include <cmath>
#include "ampmesh/MeshVariable.h"

namespace AMP {
namespace Operator {

  void PelletContactConductanceModel :: getConductance(std::vector<double> & beta,
          std::vector<double> & gamma, const std::vector<std::vector <double> > & inputVectors)
  {

    std::vector<std::string> variableNames;
    variableNames = getVariableName();

    std::vector<std::vector<double> > d_conductivity(2);

    std::vector<std::vector<double> > arguments;
    
    for (unsigned int i = 0; i < 2; i++)
    {
      arguments.resize(3);
      arguments[0].resize(inputVectors[i].size());
      for (size_t j=0; j<arguments[0].size(); j++){ 
        arguments[0][j] = inputVectors[i][j];
      }

      for (size_t l=1; l<inputVectors.size(); l++){
        for (size_t k=0; k<variableNames.size(); k++){
          if(variableNames[k]=="Concentration"){
            arguments[1].resize(inputVectors[l].size());
            for (size_t j=0; j<arguments[1].size(); j++){ 
              arguments[1][j] = inputVectors[l][j];
            }
          } else if(variableNames[k]=="Burnup"){
            arguments[2].resize(inputVectors[l].size());
            for (size_t j=0; j<arguments[2].size(); j++){ 
              arguments[2][j] = inputVectors[l][j];
            }
          }
        }
      }

      d_conductivity[i].resize(arguments[0].size());
      AMP_ERROR("FIX THIS: d_transportModels[0]->setDefaults");
      //d_transportModels[0]->setDefaults(arguments);
      AMP_ERROR("FIX THIS: d_transportModels[0]->getTransport");
      //d_transportModels[0]->getTransport(d_conductivity[i], arguments);
      arguments.resize(0);
    }

    double fkm, d_r;
    for (size_t l=1; l<inputVectors[0].size(); l++)
    {
      fkm = 2*d_conductivity[0][l]*d_conductivity[1][l]/(d_conductivity[0][l]+d_conductivity[1][l]);
      d_r = 1.414 * 1.4e-6;
      beta[l]  = 0.00125 * fkm / (d_r * exp(5.3778 - 0.528 * log(d_r*3.937*exp(7))) ) ;
      gamma[l] = 0.00125 * fkm / (d_r * exp(5.3778 - 0.528 * log(d_r*3.937*exp(7))) ) ;
    }
  }

}
}

