
#include "operators/boundary/PressureBoundaryOperator.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"
#include "vectors/VectorBuilder.h"
#include "discretization/simpleDOF_Manager.h"

#include "enum_order.h"
#include "enum_fe_family.h"
#include "enum_quadrature_type.h"
#include "string_to_enum.h"
#include "auto_ptr.h"
#include "fe_type.h"
#include "fe_base.h"
#include "quadrature.h"
#include "face_quad4.h"

namespace AMP {
  namespace Operator {

    PressureBoundaryOperator :: PressureBoundaryOperator(const boost::shared_ptr<OperatorParameters> & params)
      : BoundaryOperator(params) {
        short int bndId = (params->d_db)->getInteger("BoundaryID");
        AMP::Mesh::MeshIterator bnd = d_Mesh->getBoundaryIDIterator(AMP::Mesh::Face, bndId, 0);
        AMP::Mesh::MeshIterator end_bnd = bnd.end();

        AMP::Discretization::DOFManager::shared_ptr dofMap = AMP::Discretization::simpleDOFManager::create(bnd, 12); 

        AMP::LinearAlgebra::Variable::shared_ptr pressureVar(new AMP::LinearAlgebra::Variable("Pressure"));
        AMP::LinearAlgebra::Vector::shared_ptr pressure = AMP::LinearAlgebra::createVector(dofMap, pressureVar, true);

        libMeshEnums::Order feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>("FIRST");
        libMeshEnums::FEFamily feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>("LAGRANGE");
        libMeshEnums::QuadratureType qruleType = Utility::string_to_enum<libMeshEnums::QuadratureType>("QGAUSS");
        boost::shared_ptr < ::FEType > feType ( new ::FEType(feTypeOrder, feFamily) );
        libMeshEnums::Order qruleOrder = feType->default_quadrature_order();
        boost::shared_ptr < ::QBase > qrule( (::QBase::build(qruleType, 2, qruleOrder)).release() );

        const double val = (params->d_db)->getDouble("Value");
        for( ; bnd != end_bnd; ++bnd) {
          d_currNodes = bnd->getElements(AMP::Mesh::Vertex);
          size_t numNodesInCurrElem = d_currNodes.size();
          createCurrentLibMeshElement();

          boost::shared_ptr < ::FEBase > fe( (::FEBase::build(2, (*feType))).release() );
          fe->attach_quadrature_rule( qrule.get() );
          fe->reinit( d_currElemPtr );

          const std::vector<Point>& normals = fe->get_normals();

          std::vector<size_t> dofIndices;
          dofMap->getDOFs(bnd->globalID(), dofIndices);
          AMP_ASSERT(dofIndices.size() == 12);

          for(int d = 0; d < 3; ++d) {
            for(size_t qp = 0; qp < qrule->n_points(); ++qp) {
              pressure->setLocalValueByGlobalID((dofIndices[(3*qp) + d]), (val*(normals[qp](d))));
            }//end qp
          }//end d

          destroyCurrentLibMeshElement();
        }//end bnd

        boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));
        tmp_db->putString("InputVariable", "Pressure");
        std::string varName = params->d_db->getString("Variable");
        tmp_db->putString("OutputVariable", varName);
        tmp_db->putInteger("BoundaryID", bndId);
        tmp_db->putBool("ResidualMode", ((params->d_db)->getBool("ResidualMode")));

        boost::shared_ptr<TractionBoundaryOperatorParameters> tracOpParams(new 
            TractionBoundaryOperatorParameters(tmp_db));
        tracOpParams->d_traction = pressure;
        d_tractionOp.reset(new TractionBoundaryOperator(tracOpParams));        
      }

    void PressureBoundaryOperator :: createCurrentLibMeshElement() {
      d_currElemPtr = new ::Quad4;
      for(size_t j = 0; j < d_currNodes.size(); j++) {
        std::vector<double> pt = d_currNodes[j].coord();
        d_currElemPtr->set_node(j) = new ::Node(pt[0], pt[1], pt[2], j);
      }//end for j
    }

    void PressureBoundaryOperator :: destroyCurrentLibMeshElement() {
      for(size_t j = 0; j < d_currElemPtr->n_nodes(); j++) {
        delete (d_currElemPtr->get_node(j));
        d_currElemPtr->set_node(j) = NULL;
      }//end for j
      delete d_currElemPtr;
      d_currElemPtr = NULL;
    }

  }
}


