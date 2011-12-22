
#include "vectors/newFrozenVectorDesign/newFrozenVectorDesignHelpers.h"

#include "vectors/MultiVector.h"

namespace AMP {
  namespace LinearAlgebra {

    AMP::LinearAlgebra::Vector::shared_ptr subsetExceptForVariable(AMP::LinearAlgebra::Vector::shared_ptr inVec, 
        AMP::LinearAlgebra::Variable::shared_ptr var) {
      boost::shared_ptr<AMP::LinearAlgebra::MultiVector> outVec = boost::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVector>(
          AMP::LinearAlgebra::MultiVector::create("MultiVariable", AMP_COMM_SELF));
      boost::shared_ptr<AMP::LinearAlgebra::MultiVector> castedInVec = 
        boost::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVector>(inVec);
      AMP::LinearAlgebra::Vector::shared_ptr vec2skip = inVec->subsetVectorForVariable(var);
      AMP::LinearAlgebra::MultiVector::vector_iterator curr = castedInVec->beginVector();
      AMP::LinearAlgebra::MultiVector::vector_iterator end = castedInVec->endVector();
      for(; curr != end; curr++) {
        if((*curr) != vec2skip) {
          outVec->addVector(*curr);
        }
      }

      return outVec;
    }

    AMP::LinearAlgebra::Vector::shared_ptr joinVectors(AMP::LinearAlgebra::Vector::shared_ptr vec1, 
        AMP::LinearAlgebra::Vector::shared_ptr vec2) {
      boost::shared_ptr<AMP::LinearAlgebra::MultiVector> outVec = boost::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVector>(
          AMP::LinearAlgebra::MultiVector::create("MultiVariable", AMP_COMM_SELF));
      outVec->addVector(vec1);
      outVec->addVector(vec2);
      return outVec;
    }

  }
}


