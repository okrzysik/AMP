
#include "LinearFEOperator.h"
#include "utils/Utilities.h"

namespace AMP {
  namespace Operator {

    void LinearFEOperator :: reset(const boost::shared_ptr<OperatorParameters>& params) 
    {
      AMP_ERROR("LinearFEOperator is not converted yet");
      /*
         AMP_INSIST( ((params.get()) != NULL), "NULL parameter" );
         AMP_INSIST( (((params->d_db).get()) != NULL), "NULL database" );

         const bool reuse_matrix = (params->d_db)->getBoolWithDefault("reset_reuses_matrix", true);

         if( (d_matrix.get() == NULL) || (!reuse_matrix) ) {
         d_matrix = d_MeshAdapter->createMatrix ( this->getInputVariable(), this->getOutputVariable() );
         }

         unsigned int numDOFMaps = this->numberOfDOFMaps();
         std::vector<AMP::Mesh::DOFMap::shared_ptr> dof_maps(numDOFMaps);

         for(unsigned int i = 0; i < numDOFMaps; i++) {
         dof_maps[i] = d_MeshAdapter->getDOFMap( this->getVariableForDOFMap(i) );
         }

         AMP::Mesh::MeshManager::Adapter::ElementIterator  el = d_MeshAdapter->beginElement();
         AMP::Mesh::MeshManager::Adapter::ElementIterator  end_el = d_MeshAdapter->endElement();

         this->preAssembly(params);

         for( ; el != end_el; ++el) {

         this->preElementOperation(*el, dof_maps);

         d_elemOp->apply();

         this->postElementOperation();

         }//end for el

         this->postAssembly();
         */
    }

  }
}//end namespace


