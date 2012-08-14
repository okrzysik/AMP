
#include "operators/MeshBasedOperator.h"
#include "vectors/VectorSelector.h"

namespace AMP {
  namespace Operator {

    AMP::LinearAlgebra::Vector::shared_ptr MeshBasedOperator :: subsetOutputVector(
        AMP::LinearAlgebra::Vector::shared_ptr vec) {
      AMP::LinearAlgebra::Vector::shared_ptr varSubsetVec = vec->subsetVectorForVariable(getOutputVariable());
      if(varSubsetVec == NULL) {
        return varSubsetVec;
      } else {
        AMP::LinearAlgebra::VS_Mesh meshSelector("meshSelector", d_Mesh);
        return varSubsetVec->select(meshSelector, getOutputVariable()->getName());
      }
    }

    AMP::LinearAlgebra::Vector::shared_ptr MeshBasedOperator :: subsetInputVector(
        AMP::LinearAlgebra::Vector::shared_ptr vec) {
      AMP::LinearAlgebra::Vector::shared_ptr varSubsetVec = vec->subsetVectorForVariable(getInputVariable());
      if(varSubsetVec == NULL) {
        return varSubsetVec;
      } else {
        AMP::LinearAlgebra::VS_Mesh meshSelector("meshSelector", d_Mesh);
        return varSubsetVec->select(meshSelector, getInputVariable()->getName());
      }
    }

    AMP::LinearAlgebra::Vector::const_shared_ptr  MeshBasedOperator :: subsetOutputVector(
        AMP::LinearAlgebra::Vector::const_shared_ptr vec) {
      AMP::LinearAlgebra::Variable::shared_ptr var = getOutputVariable();
      if(d_Mesh.get() != NULL) {
        AMP::LinearAlgebra::VS_Mesh meshSelector(var->getName(), d_Mesh);
        AMP::LinearAlgebra::Vector::const_shared_ptr meshSubsetVec = vec->constSelect(meshSelector, var->getName());
        AMP::LinearAlgebra::Vector::const_shared_ptr varSubsetVec = meshSubsetVec->constSubsetVectorForVariable(var);
        return varSubsetVec;
      } else {
        return vec->constSubsetVectorForVariable(var);
      }
    }

    AMP::LinearAlgebra::Vector::const_shared_ptr  MeshBasedOperator :: subsetInputVector(
        AMP::LinearAlgebra::Vector::const_shared_ptr vec) {
      AMP::LinearAlgebra::Variable::shared_ptr var = getInputVariable();
      if(d_Mesh.get() != NULL) {
        AMP::LinearAlgebra::VS_Mesh meshSelector(var->getName(), d_Mesh);
        AMP::LinearAlgebra::Vector::const_shared_ptr meshSubsetVec = vec->constSelect(meshSelector, var->getName());
        AMP::LinearAlgebra::Vector::const_shared_ptr varSubsetVec = meshSubsetVec->constSubsetVectorForVariable(var);
        return varSubsetVec;
      } else {
        return vec->constSubsetVectorForVariable(var);
      }
    }

  }
}


