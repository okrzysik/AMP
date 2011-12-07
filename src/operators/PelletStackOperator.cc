
#include "operators/PelletStackOperator.h"

namespace AMP {
  namespace Operator {

    PelletStackOperator :: PelletStackOperator(const boost::shared_ptr<OperatorParameters> & params)
      : Operator(params) {
        d_totalNumberOfPellets = (params->d_db)->getInteger("TOTAL_NUMBER_OF_PELLETS");
        d_useSerial = (params->d_db)->getBool("USE_SERIAL");
        d_onlyZcorrection = (params->d_db)->getBool("ONLY_Z_CORRECTION");
        d_masterId = (params->d_db)->getInteger("MASTER");
        d_slaveId = (params->d_db)->getInteger("SLAVE");
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

    void PelletStackOperator :: setVariables(AMP::LinearAlgebra::Variable::shared_ptr rhs, 
        AMP::LinearAlgebra::Variable::shared_ptr sol) {
      d_rhsVar = boost::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVariable>(rhs);
      d_solVar = boost::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVariable>(sol);
    }

    void PelletStackOperator :: setLocalMeshes(std::vector<AMP::Mesh::MeshManager::Adapter::shared_ptr> inp) {
      d_meshes = inp;
    }

    void PelletStackOperator :: setLocalPelletIds(std::vector<unsigned int> inp) {
      d_pelletIds = inp;
    }

    bool PelletStackOperator :: hasPellet(unsigned int pellId) {
      return ( find(d_pelletIds.begin(), d_pelletIds.end(), pellId) != d_pelletIds.end() );
    }

    void PelletStackOperator :: apply(const AMP::LinearAlgebra::Vector::shared_ptr &f,
        const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr &r,
        const double a, const double b) {
      if(d_useSerial) {
        applySerial(f, u, r);
      } else if(d_onlyZcorrection) {
        applyOnlyZcorrection(r);
      } else {
        applyXYZcorrection(f, u, r);
      }
    }

    void PelletStackOperator :: applyUnscaling(AMP::LinearAlgebra::Vector::shared_ptr f) {
    }

    void PelletStackOperator :: applySerial(const AMP::LinearAlgebra::Vector::shared_ptr &f,
        const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr  &r) {
    }

    void PelletStackOperator :: applyOnlyZcorrection(AMP::LinearAlgebra::Vector::shared_ptr &r) {
    }

    void PelletStackOperator :: applyXYZcorrection(const AMP::LinearAlgebra::Vector::shared_ptr &f,
        const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr  &r) {
    }

  }
}


