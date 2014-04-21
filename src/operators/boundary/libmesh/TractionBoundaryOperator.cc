
#include "operators/boundary/libmesh/TractionBoundaryOperator.h"

#include "enum_order.h"
#include "enum_fe_family.h"
#include "enum_quadrature_type.h"
#include "string_to_enum.h"
#include "auto_ptr.h"
#include "fe_type.h"
#include "fe_base.h"
#include "quadrature.h"
#include "elem.h"
#include "cell_hex8.h"
#include "node.h"

namespace AMP {
  namespace Operator {

    TractionBoundaryOperator :: TractionBoundaryOperator(const boost::shared_ptr<TractionBoundaryOperatorParameters> & params)
      : BoundaryOperator(params) {
        AMP_INSIST( params->d_db->keyExists("Variable"), "key not found");
        std::string varName = params->d_db->getString("Variable");
        d_var.reset(new AMP::LinearAlgebra::Variable(varName));
        AMP_INSIST( params->d_db->keyExists("ResidualMode"), "key not found");
        d_residualMode = params->d_db->getBool("ResidualMode");
        d_traction = params->d_traction;
        d_volumeElements = params->d_volumeElements;
        d_sideNumbers = params->d_sideNumbers;
        d_nodeID = params->d_nodeID;
      }

    void TractionBoundaryOperator :: addRHScorrection(AMP::LinearAlgebra::Vector::shared_ptr rhs) {
      if(!d_residualMode) {
        AMP::LinearAlgebra::Vector::shared_ptr myRhs = mySubsetVector(rhs, d_var);
        if(d_correction == NULL) {
          d_correction = myRhs->cloneVector();
        }
        computeCorrection();
        myRhs->add(myRhs, d_correction);
        myRhs->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
      }
    }

    void TractionBoundaryOperator :: apply( AMP::LinearAlgebra::Vector::const_shared_ptr, AMP::LinearAlgebra::Vector::const_shared_ptr,
        AMP::LinearAlgebra::Vector::shared_ptr r, const double, const double) {
      if(d_residualMode) {
        AMP::LinearAlgebra::Vector::shared_ptr rInternal = mySubsetVector(r, d_var);
        if(d_correction == NULL) {
          d_correction = rInternal->cloneVector();
        }
        computeCorrection();
        rInternal->subtract(rInternal, d_correction);
        rInternal->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
      }
    }

    void TractionBoundaryOperator :: computeCorrection() {
      libMeshEnums::Order feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>("FIRST");
      libMeshEnums::FEFamily feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>("LAGRANGE");
      libMeshEnums::QuadratureType qruleType = Utility::string_to_enum<libMeshEnums::QuadratureType>("QGAUSS");
      boost::shared_ptr < ::FEType > feType ( new ::FEType(feTypeOrder, feFamily) );
      libMeshEnums::Order qruleOrder = feType->default_quadrature_order();
      boost::shared_ptr < ::QBase > qrule( (::QBase::build(qruleType, 2, qruleOrder)).release() );

      AMP::Discretization::DOFManager::shared_ptr dofMap = d_correction->getDOFManager();

      d_correction->zero();
      for(size_t b = 0; b < d_sideNumbers.size(); ++b) {
        ::Elem* elem = new ::Hex8;
        for(int j = 0; j < 8; ++j) {
          elem->set_node(j) = new ::Node(d_volumeElements[(24*b) + (3*j) + 0],
              d_volumeElements[(24*b) + (3*j) + 1], d_volumeElements[(24*b) + (3*j) + 2], j);
        }//end j

        boost::shared_ptr < ::FEBase > fe( (::FEBase::build(3, (*feType))).release() );
        fe->attach_quadrature_rule( qrule.get() );
        fe->reinit(elem, d_sideNumbers[b]);

        const std::vector<std::vector<Real> > &phi = fe->get_phi();
        const std::vector<Real> &djxw = fe->get_JxW(); 

        AMP_ASSERT(phi.size() == 8);
        AMP_ASSERT(djxw.size() == 4);
        AMP_ASSERT(qrule->n_points() == 4);

        std::vector<std::vector<size_t> > dofIndices(8);
        for(int i = 0; i < 8; ++i) {
          dofMap->getDOFs(d_nodeID[(8*b) + i], dofIndices[i]);
        }//end i

        for(size_t i = 0; i < 8; ++i) {
          for(int d = 0; d < 3; ++d) {
            double res = 0;
            for(size_t qp = 0; qp < qrule->n_points(); ++qp) {
              double val = d_traction[(12*b) + (3*qp) + d];
              res +=  djxw[qp]*phi[i][qp]*val;
            }//end qp
            d_correction->addValueByGlobalID(dofIndices[i][d], res);
          }//end d
        }//end i

        for(size_t j = 0; j < elem->n_nodes(); ++j) {
          delete (elem->get_node(j));
          elem->set_node(j) = NULL;
        }//end for j
        delete elem;
        elem = NULL;
      }//end b
      d_correction->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_ADD );
    }

    AMP::LinearAlgebra::Vector::shared_ptr TractionBoundaryOperator :: mySubsetVector(AMP::LinearAlgebra::Vector::shared_ptr vec, 
        AMP::LinearAlgebra::Variable::shared_ptr var) {
      if(vec != NULL) {
        if(d_Mesh.get() != NULL) {
          AMP::LinearAlgebra::VS_Mesh meshSelector(d_Mesh);
          AMP::LinearAlgebra::Vector::shared_ptr meshSubsetVec = vec->select(meshSelector, ((vec->getVariable())->getName()));
          return meshSubsetVec->subsetVectorForVariable(var);
        } else {
          return vec->subsetVectorForVariable(var);
        }
      } else {
        return vec;
      }
    }

  }
}

