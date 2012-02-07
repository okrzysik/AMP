
#include "NonlinearFEOperator.h"
#include "utils/Utilities.h"
#include "cell_hex8.h"
#include "node.h"

namespace AMP {
  namespace Operator {

    void NonlinearFEOperator :: apply(const boost::shared_ptr<AMP::LinearAlgebra::Vector>  &f, 
        const boost::shared_ptr<AMP::LinearAlgebra::Vector>  &u, boost::shared_ptr<AMP::LinearAlgebra::Vector>  &r,
        const double a,  const double b)
    {
      AMP_INSIST( (r != NULL), "NULL Residual/Output Vector" );

      AMP::LinearAlgebra::Vector::shared_ptr rInternal = this->subsetOutputVector(r);

      AMP::Mesh::MeshIterator el = d_Mesh->getIterator(AMP::Mesh::Volume, 0);
      AMP::Mesh::MeshIterator end_el = el.end();

      this->preAssembly(u, rInternal);

      for( ; el != end_el; ++el) {
        this->preElementOperation(*el);

        d_elemOp->apply();

        this->postElementOperation();
      }//end for el

      this->postAssembly();

      if(f == NULL) {
        rInternal->scale(a);
      } else {
        AMP::LinearAlgebra::Vector::shared_ptr fInternal = this->subsetOutputVector(f);
        if(fInternal == NULL) {
          rInternal->scale(a);
        } else {
          rInternal->axpby(b, a, fInternal);
        }
      }

      if(d_iDebugPrintInfoLevel>2)
      {
        AMP::pout << "L2 norm of result of NonlinearFEOperator::apply is: " << rInternal->L2Norm() << std::endl;
      }
      if(d_iDebugPrintInfoLevel>5)
      {
        std::cout << rInternal << std::endl;
      }
    }

    void NonlinearFEOperator :: createCurrentLibMeshElement() {
      d_currElemPtr = new ::Hex8;
      for(size_t j = 0; j < d_currNodes.size(); j++) {
        std::vector<double> pt = d_currNodes[j].coord();
        d_currElemPtr->set_node(j) = new ::Node(pt[0], pt[1], pt[2], j);
      }//end for j
    }

    void NonlinearFEOperator :: destroyCurrentLibMeshElement() {
      for(size_t j = 0; j < d_currElemPtr->n_nodes(); j++) {
        delete (d_currElemPtr->get_node(j));
        d_currElemPtr->set_node(j) = NULL;
      }//end for j
      delete d_currElemPtr;
      d_currElemPtr = NULL;
    }

  }
}//end namespace


