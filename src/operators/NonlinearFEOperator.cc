
#include "NonlinearFEOperator.h"
#include "utils/Utilities.h"

namespace AMP {
  namespace Operator {

    void NonlinearFEOperator :: apply(const boost::shared_ptr<AMP::LinearAlgebra::Vector>  &f, 
        const boost::shared_ptr<AMP::LinearAlgebra::Vector>  &u, boost::shared_ptr<AMP::LinearAlgebra::Vector>  &r,
        const double a,  const double b)
    {
AMP_ERROR("NonlinearFEOperator is not converted yet");
/*

      d_applyCount++;

      AMP_INSIST( (r != NULL), "NULL Residual/Output Vector" );

      unsigned int numDOFMaps = this->numberOfDOFMaps();
      std::vector<AMP::Mesh::DOFMap::shared_ptr> dof_maps(numDOFMaps);

      for(unsigned int i = 0; i < numDOFMaps; i++) {
        dof_maps[i] = d_MeshAdapter->getDOFMap( this->getVariableForDOFMap(i) );
      }

      AMP::Mesh::MeshManager::Adapter::ElementIterator  el = d_MeshAdapter->beginElement();
      AMP::Mesh::MeshManager::Adapter::ElementIterator  end_el = d_MeshAdapter->endElement();

      AMP::LinearAlgebra::Vector::shared_ptr rInternal = r->subsetVectorForVariable( this->getOutputVariable() );

      this->preAssembly(u, rInternal);

      for( ; el != end_el; ++el) {
        this->preElementOperation(*el, dof_maps);

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
*/
    }

  }
}//end namespace


