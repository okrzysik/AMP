
#include "operators/boundary/TractionBoundaryOperator.h"

namespace AMP {
  namespace Operator {

    TractionBoundaryOperator :: TractionBoundaryOperator(const boost::shared_ptr<TractionBoundaryOperatorParameters> & params)
      : BoundaryOperator(params) {
        AMP_INSIST( params->d_db->keyExists("InputVariable"), "key not found");
        AMP_INSIST( params->d_db->keyExists("OutputVariable"), "key not found");
        std::string inpVarName = params->d_db->getString("InputVariable");
        std::string outVarName = params->d_db->getString("OutputVariable");
        d_inputVar.reset(new AMP::LinearAlgebra::Variable(inpVarName));
        d_outputVar.reset(new AMP::LinearAlgebra::Variable(outVarName));
        d_residualMode = params->d_db->getBool("ResidualMode");
        d_boundaryId = params->d_db->getInteger("BoundaryID");
        setTractionVec(params->d_tractionVec);
      }

    void TractionBoundaryOperator :: addRHScorrection(AMP::LinearAlgebra::Vector::shared_ptr rhs) {
      if(!d_residualMode) {
        AMP::LinearAlgebra::Vector::shared_ptr myRhs = mySubsetVector(rhs, d_outputVar);
        if(d_correction == NULL) {
          d_correction = myRhs->cloneVector();
        }
        computeCorrection();
        myRhs->add(myRhs, d_correction);
      }
    }

    void TractionBoundaryOperator :: apply(const AMP::LinearAlgebra::Vector::shared_ptr &, const AMP::LinearAlgebra::Vector::shared_ptr &,
        AMP::LinearAlgebra::Vector::shared_ptr &r, const double, const double) {
      if(d_residualMode) {
        AMP::LinearAlgebra::Vector::shared_ptr rInternal = mySubsetVector(r, d_outputVar);
        if(d_correction == NULL) {
          d_correction = rInternal->cloneVector();
        }
        computeCorrection();
        rInternal->subtract(rInternal, d_correction);
      }
    }


    void TractionBoundaryOperator :: computeCorrection() {

    }

    AMP::LinearAlgebra::Vector::shared_ptr TractionBoundaryOperator :: mySubsetVector(AMP::LinearAlgebra::Vector::shared_ptr vec, 
        AMP::LinearAlgebra::Variable::shared_ptr var) {
      if(vec != NULL) {
        if(d_Mesh.get() != NULL) {
          AMP::LinearAlgebra::VS_Mesh meshSelector(var->getName(), d_Mesh);
          AMP::LinearAlgebra::Vector::shared_ptr meshSubsetVec = vec->select(meshSelector, var->getName());
          return meshSubsetVec->subsetVectorForVariable(var);
        } else {
          return vec->subsetVectorForVariable(var);
        }
      } else {
        return vec;
      }
    }

  }
}

