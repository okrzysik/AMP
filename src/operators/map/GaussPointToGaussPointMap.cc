
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
    }

    void GaussPointToGaussPointMap :: createIdxMap(boost::shared_ptr<AMP::Database> db) {
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
      const ::Elem *elem; 

      boost::shared_ptr < ::FEBase > fe( (::FEBase::build(dimension, (*feType))).release() ); 

      fe->attach_quadrature_rule( qrule.get() );

    }

  }
}


