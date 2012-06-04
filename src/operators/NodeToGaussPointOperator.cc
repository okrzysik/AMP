#include "operators/NodeToGaussPointOperator.h"

namespace AMP {
namespace Operator {


    void NodeToGaussPointOperator :: apply(const AMP::LinearAlgebra::Vector::shared_ptr &, const AMP::LinearAlgebra::Vector::shared_ptr &u,
        AMP::LinearAlgebra::Vector::shared_ptr  &r, const double , const double )
    { 

      AMP::LinearAlgebra::Vector::shared_ptr nodalVec = u->subsetVectorForVariable(d_NodalVariable);
      AMP::LinearAlgebra::Vector::shared_ptr gaussPtVec = r->subsetVectorForVariable(d_GaussPtVariable);

      AMP::Discretization::DOFManager::shared_ptr gaussPt_dof_maps = gaussPtVec->getDOFManager();

      libMeshEnums::Order feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>("FIRST");
      libMeshEnums::FEFamily feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>("LAGRANGE");

      boost::shared_ptr < ::FEType > d_feType ( new ::FEType(feTypeOrder, feFamily) );
      boost::shared_ptr < ::FEBase > d_fe ( (::FEBase::build(2, (*d_feType))).release() );

      libMeshEnums::Order qruleOrder = Utility::string_to_enum<libMeshEnums::Order>("SECOND");
      boost::shared_ptr < ::QBase > d_qrule ( (::QBase::build("QGAUSS", 2, qruleOrder)).release() );

      d_fe->attach_quadrature_rule( d_qrule.get() );

      AMP::Mesh::MeshIterator el = d_Mesh->getIterator(AMP::Mesh::Volume, 0);
      AMP::Mesh::MeshIterator end_el = el.end();

      AMP::Discretization::DOFManager::shared_ptr dof_map = TemperatureVec->getDOFManager();

      for( ; el != end_el; ++el) {
        std::vector<AMP::Mesh::MeshElement> d_currNodes = el->getElements(AMP::Mesh::Vertex);  

        ::Elem* d_currElemPtr = new ::Hex8;
        for(size_t j = 0; j < d_currNodes.size(); j++) {
          std::vector<double> pt = d_currNodes[j].coord();
          d_currElemPtr->set_node(j) = new ::Node(pt[0], pt[1], pt[2], j);
        }//end for j

        std::vector<size_t> bndGlobalIds;
        std::vector<AMP::Mesh::MeshElementID> globalIds(d_currNodes.size()) ;
        for (size_t i = 0; i<d_currNodes.size(); i++)
        {
          globalIds[i] = d_currNodes[i].globalID();
        }
        dof_map->getDOFs(globalIds, bndGlobalIds);

        std::vector<size_t> d_gaussPtIndices; 
        gaussPt_dof_maps->getDOFs (el->globalID(), d_gaussPtIndices);
        d_fe->reinit(d_currElemPtr);

        const std::vector<Real> & JxW = d_fe->get_JxW();
        const std::vector<std::vector<Real> > & phi = d_fe->get_phi();

        std::vector<double>  computedAtGauss(d_qrule->n_points(), 0.0);
        for (unsigned int qp = 0; qp < d_qrule->n_points(); qp++){
          for (unsigned int j = 0; j < bndGlobalIds.size(); j++){
            double computedAtNode = nodalVec->getValueByGlobalID(bndGlobalIds[j]); 
            computedAtGauss[qp]     += computedAtNode * phi[j][qp];
          }//end for j
        }//end for qp

        for (unsigned int qp = 0; qp < d_gaussPtIndices.size(); qp++) {
          gaussPtVec->setValueByGlobalID( d_gaussPtIndices[qp] , computedAtGauss[qp]);
        }
        for(size_t j = 0; j < d_currElemPtr->n_nodes(); j++) {
          delete (d_currElemPtr->get_node(j));
          d_currElemPtr->set_node(j) = NULL;
        }//end for j
        delete d_currElemPtr;
        d_currElemPtr = NULL;

      }

    }// end apply

}
}

