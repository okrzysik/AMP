
#include "NonlinearFEOperator.h"
#include "utils/Utilities.h"

namespace AMP {
  namespace Operator {

    void NonlinearFEOperator :: apply(const boost::shared_ptr<AMP::LinearAlgebra::Vector>  &f, 
        const boost::shared_ptr<AMP::LinearAlgebra::Vector>  &u, boost::shared_ptr<AMP::LinearAlgebra::Vector>  &r,
        const double a,  const double b)
    {
      AMP_INSIST( (r != NULL), "NULL Residual/Output Vector" );

      AMP::LinearAlgebra::Vector::shared_ptr rInternal = r->subsetVectorForVariable( this->getOutputVariable() );

      AMP::Mesh::MeshIterator el = d_Mesh->getIterator(AMP::Mesh::Volume, 0);
      AMP::Mesh::MeshIterator end_el = el.end();

      this->preAssembly(u, rInternal);

      for( ; el != end_el; ++el) {
        this->preElementOperation(*el);

        d_elemOp->apply();

        this->postElementOperation();
      }//end for el

      this->postAssembly();

      if(f.get() == NULL) {
        rInternal->scale(a);
      } else {
        AMP::LinearAlgebra::Vector::shared_ptr fInternal = f->subsetVectorForVariable( this->getOutputVariable() );
        if(fInternal.get() == NULL) {
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

  }
}//end namespace


