
#include "LinearFEOperator.h"
#include "utils/Utilities.h"
#include "matrices/MatrixBuilder.h"
#include "vectors/VectorBuilder.h"
#include "utils/ProfilerApp.h"

#include "cell_hex8.h"
#include "node.h"

namespace AMP {
  namespace Operator {

    void LinearFEOperator :: reset(const boost::shared_ptr<OperatorParameters>& params) 
    {
      PROFILE_START("reset");
      AMP_INSIST( (params != NULL), "NULL parameter" );
      AMP_INSIST( ((params->d_db) != NULL), "NULL database" );

      const bool reuse_matrix = (params->d_db)->getBoolWithDefault("reset_reuses_matrix", true);

      if( (d_matrix.get() == NULL) || (!reuse_matrix) ) {
        AMP::LinearAlgebra::Vector::shared_ptr inVec = AMP::LinearAlgebra::createVector(d_inDofMap, getInputVariable(), true);
        AMP::LinearAlgebra::Vector::shared_ptr outVec = AMP::LinearAlgebra::createVector(d_outDofMap, getOutputVariable(), true);
        d_matrix = AMP::LinearAlgebra::createMatrix(inVec, outVec);
        d_matrix->zero();
        d_matrix->makeConsistent();
      }
      AMP_ASSERT((*d_inDofMap)==(*d_matrix->getLeftDOFManager()));
      AMP_ASSERT((*d_inDofMap)==(*d_matrix->getRightDOFManager()));

      AMP::Mesh::MeshIterator el = d_Mesh->getIterator(AMP::Mesh::Volume, 0);
      AMP::Mesh::MeshIterator end_el = el.end();

      this->preAssembly(params);

      for( ; el != end_el; ++el) {

        this->preElementOperation(*el);

        d_elemOp->apply();

        this->postElementOperation();

      }//end for el

      this->postAssembly();
      PROFILE_STOP("reset");
    }

    void LinearFEOperator :: createCurrentLibMeshElement() {
      d_currElemPtr = new ::Hex8;
      for(size_t j = 0; j < d_currNodes.size(); ++j) {
        std::vector<double> pt = d_currNodes[j].coord();
        d_currElemPtr->set_node(j) = new ::Node(pt[0], pt[1], pt[2], j);
      }//end for j
    }

    void LinearFEOperator :: destroyCurrentLibMeshElement() {
      for(size_t j = 0; j < d_currElemPtr->n_nodes(); ++j) {
        delete (d_currElemPtr->get_node(j));
        d_currElemPtr->set_node(j) = NULL;
      }//end for j
      delete d_currElemPtr;
      d_currElemPtr = NULL;
    }


  }
}//end namespace


