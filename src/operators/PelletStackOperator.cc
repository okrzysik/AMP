
#include "operators/PelletStackOperator.h"

namespace AMP {
  namespace Operator {

    PelletStackOperator :: PelletStackOperator(const boost::shared_ptr<OperatorParameters> & params)
      : Operator(params) {
        d_totalNumberOfPellets = (params->d_db)->getInteger("TOTAL_NUMBER_OF_PELLETS");
        d_useSerial = (params->d_db)->getBool("USE_SERIAL");
        d_onlyZcorrection = (params->d_db)->getBool("ONLY_Z_CORRECTION");
        if((params->d_db)->keyExists("SCALING_FACTOR")) {
          d_useScaling = true;
          d_scalingFactor = (params->d_db)->getDouble("SCALING_FACTOR");
        } else {
          d_useScaling = false;
        }
        d_currentPellet = 0;
      }

    bool PelletStackOperator :: useSerial() {
      return d_useSerial;
    }

    bool PelletStackOperator :: onlyZcorrection() {
      return d_onlyZcorrection;
    }

    bool PelletStackOperator :: useScaling() {
      return d_useScaling;
    }

    void PelletStackOperator :: setCurrentPellet(unsigned int pellId) {
      d_currentPellet = pellId;
    }

    unsigned int PelletStackOperator :: getTotalNumberOfPellets() {
      return d_totalNumberOfPellets;
    }

    bool PelletStackOperator :: hasPellet(unsigned int pellId) {
      return true;
    }

    void PelletStackOperator :: applyUnscaling(AMP::LinearAlgebra::Vector::shared_ptr f) {
    }

    AMP::LinearAlgebra::Variable::shared_ptr PelletStackOperator :: getOutputVariable() {
      AMP::LinearAlgebra::Variable::shared_ptr emptyPointer;
      return emptyPointer;
    }

    AMP::LinearAlgebra::Variable::shared_ptr PelletStackOperator :: getInputVariable(int varId) {
      AMP::LinearAlgebra::Variable::shared_ptr emptyPointer;
      return emptyPointer;
    }

    void PelletStackOperator :: apply(const AMP::LinearAlgebra::Vector::shared_ptr &f,
        const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr &r,
        const double a, const double b) {
    }

  }
}


