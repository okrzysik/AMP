
#include "LinearFEOperator.h"
#include "utils/Utilities.h"
#include "matrices/MatrixBuilder.h"
#include "vectors/VectorBuilder.h"

namespace AMP {
  namespace Operator {

    void LinearFEOperator :: reset(const boost::shared_ptr<OperatorParameters>& params) 
    {
      AMP_INSIST( (params != NULL), "NULL parameter" );
      AMP_INSIST( ((params->d_db) != NULL), "NULL database" );

      const bool reuse_matrix = (params->d_db)->getBoolWithDefault("reset_reuses_matrix", true);

      boost::shared_ptr<FEOperatorParameters> feOpParams = boost::dynamic_pointer_cast<FEOperatorParameters>(params);

      if( (d_matrix.get() == NULL) || (!reuse_matrix) ) {
        AMP::LinearAlgebra::Vector::shared_ptr inVec = AMP::LinearAlgebra::createVector((feOpParams->d_inDofMap), getInputVariable(), false);
        AMP::LinearAlgebra::Vector::shared_ptr outVec = AMP::LinearAlgebra::createVector((feOpParams->d_outDofMap), getOutputVariable(), false);
        d_matrix = AMP::LinearAlgebra::createMatrix(inVec, outVec);
      }

      AMP::Mesh::MeshIterator  el = d_Mesh->getIterator(AMP::Mesh::Volume, 0);
      AMP::Mesh::MeshIterator  end_el = el.end();

      this->preAssembly(params);

      for( ; el != end_el; ++el) {

        this->preElementOperation(*el);

        d_elemOp->apply();

        this->postElementOperation();

      }//end for el

      this->postAssembly();
    }

  }
}//end namespace


