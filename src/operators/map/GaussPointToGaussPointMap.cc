
#include "operators/map/GaussPointToGaussPointMap.h"

#include "elem.h"
#include "fe_type.h"
#include "fe_base.h"
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

      int dimension = db->getIntegerWithDefault("DIMENSION", 2);

      boost::shared_ptr < ::FEType > feType(new ::FEType(feTypeOrder, feFamily) ); 

      libMeshEnums::Order qruleOrder;

      if(qruleOrderName == "DEFAULT") {
        qruleOrder = feType->default_quadrature_order();
      } else {
        qruleOrder = Utility::string_to_enum<libMeshEnums::Order>(qruleOrderName);
      }

      boost::shared_ptr < ::QBase > qrule( (::QBase::build(qruleType, dimension, qruleOrder)).release() ); 
      ::Elem *elem = NULL; 

      boost::shared_ptr < ::FEBase > fe( (::FEBase::build(dimension, (*feType))).release() ); 
      fe->reinit(elem);

      fe->attach_quadrature_rule( qrule.get() );
    }

  }
}


