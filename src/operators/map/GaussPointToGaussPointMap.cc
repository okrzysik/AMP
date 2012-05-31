
#include "operators/map/GaussPointToGaussPointMap.h"
#include "vectors/VectorBuilder.h"
#include "ampmesh/MultiMesh.h"
#include "discretization/simpleDOF_Manager.h"

#include "elem.h"
#include "fe_type.h"
#include "fe_base.h"
#include "face_quad4.h"
#include "quadrature.h"

#include "enum_order.h"
#include "enum_fe_family.h"
#include "enum_quadrature_type.h"
#include "auto_ptr.h"
#include "string_to_enum.h"

#include <string>

namespace AMP {
  namespace Operator {

    void GaussPointToGaussPointMap :: correctLocalOrdering() {
      AMP::Discretization::DOFManager::shared_ptr dofMap = d_OutputVector->getDOFManager();
      std::vector<size_t> localDofs(DofsPerObj);
      for(size_t i = 0; i < d_recvList.size(); ++i) {
        dofMap->getDOFs( d_recvList[i], localDofs );
        std::vector<double> vals(DofsPerObj);
        for(int j = 0; j < DofsPerObj; ++j) {
          vals[j] = d_OutputVector->getLocalValueByGlobalID(localDofs[j]);
        }//end j
        int DofsPerGaussPt = DofsPerObj/(d_idxMap[i].size());
        for(size_t j = 0; j < d_idxMap[i].size(); ++j) {
          for(int k = 0; k < DofsPerGaussPt; ++k) {
            d_OutputVector->setLocalValueByGlobalID(localDofs[(j*DofsPerGaussPt) + k],
                vals[((d_idxMap[i][j])*DofsPerGaussPt) + k]);
          }//end k
        }//end j
      }//end i
    }

    void GaussPointToGaussPointMap :: createIdxMap(boost::shared_ptr<AMP::Operator::OperatorParameters> params) {
      boost::shared_ptr<AMP::Database> db = params->d_db;
      std::string feTypeOrderName = db->getStringWithDefault("FE_ORDER", "FIRST");
      libMeshEnums::Order feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>(feTypeOrderName);

      std::string feFamilyName = db->getStringWithDefault("FE_FAMILY", "LAGRANGE");
      libMeshEnums::FEFamily feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>(feFamilyName);

      std::string qruleTypeName = db->getStringWithDefault("QRULE_TYPE", "QGAUSS");
      libMeshEnums::QuadratureType qruleType = Utility::string_to_enum<libMeshEnums::QuadratureType>(qruleTypeName);

      std::string qruleOrderName = db->getStringWithDefault("QRULE_ORDER", "DEFAULT");

      int faceDim = db->getIntegerWithDefault("DIMENSION", 2);

      boost::shared_ptr < ::FEType > feType(new ::FEType(feTypeOrder, feFamily) ); 

      libMeshEnums::Order qruleOrder;

      if(qruleOrderName == "DEFAULT") {
        qruleOrder = feType->default_quadrature_order();
      } else {
        qruleOrder = Utility::string_to_enum<libMeshEnums::Order>(qruleOrderName);
      }

      boost::shared_ptr < ::QBase > qrule( (::QBase::build(qruleType, faceDim, qruleOrder)).release() ); 

      int numGaussPtsPerElem = qrule->n_points();

      int dofsPerElem = (dim*numGaussPtsPerElem);

      AMP::LinearAlgebra::Variable::shared_ptr variable(new AMP::LinearAlgebra::Variable("GaussPoints"));

      std::vector<AMP::Mesh::Mesh::shared_ptr> meshesForMap(2);
      meshesForMap[0] = d_mesh1;
      meshesForMap[1] = d_mesh2;
      AMP::Mesh::Mesh::shared_ptr multiMesh(new AMP::Mesh::MultiMesh(d_MapComm, meshesForMap));

      AMP::Discretization::DOFManager::shared_ptr dofMap = AMP::Discretization::simpleDOFManager::create(multiMesh,
          AMP::Mesh::Face, 0, dofsPerElem);

      AMP::LinearAlgebra::Vector::shared_ptr inVec = AMP::LinearAlgebra::createVector(dofMap, variable);

      AMP::LinearAlgebra::Vector::shared_ptr outVec = inVec->cloneVector();

      std::vector<size_t> localDofs(dofsPerElem);
      for(size_t i = 0; i < d_sendList.size(); ++i) {
        AMP::Mesh::MeshElement el = multiMesh->getElement(d_sendList[i]);

        dofMap->getDOFs( d_sendList[i], localDofs );

        std::vector<AMP::Mesh::MeshElement> currNodes = el.getElements(AMP::Mesh::Vertex);

        ::Elem *elem = new ::Quad4; 
        for(size_t j = 0; j < currNodes.size(); ++j) {
          std::vector<double> pt = currNodes[j].coord();
          elem->set_node(j) = new ::Node(pt[0], pt[1], pt[2], j);
        }//end for j

        boost::shared_ptr < ::FEBase > fe( (::FEBase::build(faceDim, (*feType))).release() ); 
        fe->attach_quadrature_rule( qrule.get() );
        fe->reinit(elem);

        const std::vector< ::Point > & xyz = fe->get_xyz();
        for(size_t j = 0; j < numGaussPtsPerElem; ++j) {
          for(int k = 0; k < dim; ++k) {
            inVec->setLocalValueByGlobalID(localDofs[(j*dim) + k], xyz[j](k));
          }//end for k
        }//end for j

        for(size_t j = 0; j < elem->n_nodes(); ++j) {
          delete (elem->get_node(j));
          elem->set_node(j) = NULL;
        }//end for j
        delete elem;
        elem = NULL;
      }//end i

      db->putInteger("DOFsPerObject", dofsPerElem);
      db->putString("VariableName", "GaussPoints");
      boost::shared_ptr<AMP::Operator::NodeToNodeMap> n2nMap(new AMP::Operator::NodeToNodeMap(params));
      n2nMap->setVector(outVec);

      AMP::LinearAlgebra::Vector::shared_ptr nullVec;
      n2nMap->apply(nullVec, inVec, nullVec, 1, 0);

    }

  }
}


