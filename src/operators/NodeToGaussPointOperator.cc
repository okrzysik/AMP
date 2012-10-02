#include "operators/NodeToGaussPointOperator.h"
#include "utils/ProfilerApp.h"

/* Libmesh files */
#include "fe_type.h"
#include "fe_base.h"
#include "elem.h"
#include "quadrature.h"

#include "enum_order.h"
#include "enum_fe_family.h"
#include "enum_quadrature_type.h"
#include "auto_ptr.h"
#include "string_to_enum.h"

#include "face_quad4.h"
#include "cell_hex8.h"
#include "node.h"

namespace AMP {
namespace Operator {


void NodeToGaussPointOperator :: apply(AMP::LinearAlgebra::Vector::const_shared_ptr,
        AMP::LinearAlgebra::Vector::const_shared_ptr u,
        AMP::LinearAlgebra::Vector::shared_ptr r, const double , const double ) 
{ 
    PROFILE_START("apply");

    AMP::LinearAlgebra::Vector::const_shared_ptr nodalVec = subsetInputVector(u);
    AMP::LinearAlgebra::Vector::shared_ptr gaussPtVec = subsetOutputVector(r);

    AMP_ASSERT(nodalVec->getUpdateStatus()==AMP::LinearAlgebra::Vector::UNCHANGED);

    AMP::Discretization::DOFManager::shared_ptr dof_map = nodalVec->getDOFManager();
    AMP::Discretization::DOFManager::shared_ptr gaussPt_dof_map = gaussPtVec->getDOFManager();

    libMeshEnums::Order feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>("FIRST");
    libMeshEnums::FEFamily feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>("LAGRANGE");

    boost::shared_ptr < ::FEType > d_feType ( new ::FEType(feTypeOrder, feFamily) );

    int dim;
    if(d_UseSurfaceElements) {
        dim = 2;
    } else {
        dim = 3;
    }

    libMeshEnums::Order qruleOrder = Utility::string_to_enum<libMeshEnums::Order>("SECOND");
    boost::shared_ptr < ::QBase > d_qrule ( (::QBase::build("QGAUSS", dim, qruleOrder)).release() );

    AMP::Mesh::MeshIterator el;
    if(d_UseSurfaceElements) {
        el = d_Mesh->getIterator(AMP::Mesh::Face, 0);
    } else {
        el = d_Mesh->getIterator(AMP::Mesh::Volume, 0);
    }

    AMP::Mesh::MeshIterator end_el = el.end();

    for( ; el != end_el; ++el) {
        boost::shared_ptr < ::FEBase > d_fe ( (::FEBase::build(dim, (*d_feType))).release() );
        d_fe->attach_quadrature_rule( d_qrule.get() );

        std::vector<AMP::Mesh::MeshElement> d_currNodes = el->getElements(AMP::Mesh::Vertex);  

        ::Elem* d_currElemPtr;
        if(d_UseSurfaceElements) {
          d_currElemPtr = new ::Quad4;
        } else {
          d_currElemPtr = new ::Hex8;
        }
        for(size_t j = 0; j < d_currNodes.size(); ++j) {
          std::vector<double> pt = d_currNodes[j].coord();
          d_currElemPtr->set_node(j) = new ::Node(pt[0], pt[1], pt[2], j);
        }//end for j

        std::vector<size_t> bndGlobalIds;
        std::vector<AMP::Mesh::MeshElementID> globalIds(d_currNodes.size()) ;
        for (size_t i = 0; i < d_currNodes.size(); ++i) {
          globalIds[i] = d_currNodes[i].globalID();
        }
        dof_map->getDOFs(globalIds, bndGlobalIds);

        std::vector<size_t> d_gaussPtIndices; 
        gaussPt_dof_map->getDOFs (el->globalID(), d_gaussPtIndices);
        d_fe->reinit(d_currElemPtr);

        const std::vector<std::vector<Real> > & phi = d_fe->get_phi();

        std::vector<double>  computedAtGauss(d_qrule->n_points(), 0.0);
        for (unsigned int j = 0; j < bndGlobalIds.size(); ++j){
          double computedAtNode = nodalVec->getValueByGlobalID(bndGlobalIds[j]); 
          for (unsigned int qp = 0; qp < d_qrule->n_points(); ++qp){
            computedAtGauss[qp] += ( computedAtNode * phi[j][qp] );
          }//end for qp
        }//end for j

        for (unsigned int qp = 0; qp < d_gaussPtIndices.size(); ++qp) {
          gaussPtVec->setLocalValueByGlobalID( d_gaussPtIndices[qp] , computedAtGauss[qp]);
        }
        for(size_t j = 0; j < d_currElemPtr->n_nodes(); ++j) {
          delete (d_currElemPtr->get_node(j));
          d_currElemPtr->set_node(j) = NULL;
        }//end for j
        delete d_currElemPtr;
        d_currElemPtr = NULL;
    }//end for
    PROFILE_STOP("apply");
}// end apply


}
}



