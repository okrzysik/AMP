
#include "FlowElement.h"
#include "utils/Utilities.h"

/* Libmesh files */
#include "enum_order.h"
#include "enum_fe_family.h"
#include "enum_quadrature_type.h"
#include "auto_ptr.h"
#include "string_to_enum.h"

#include <string>

namespace AMP {
namespace Operator {

  FlowElement :: FlowElement (const boost::shared_ptr<ElementOperationParameters> & params)
    : ElementOperation(params)
  {
    AMP_INSIST( (params.get() != NULL), "''params'' is NULL");

    AMP_INSIST( (((params->d_db).get()) != NULL), "NULL database" );

    const unsigned int dimension = (params->d_db)->getIntegerWithDefault("DIMENSION", 3);
    // int numApprox = (params->d_db)->getIntegerWithDefault("NUM_APPROX", 2);
    
    std::string U_feTypeOrderName = (params->d_db)->getStringWithDefault("FE_ORDER", "SECOND");
    libMeshEnums::Order feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>(U_feTypeOrderName);

    std::string feFamilyName = (params->d_db)->getStringWithDefault("FE_FAMILY", "LAGRANGE");
    libMeshEnums::FEFamily feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>(feFamilyName);

    std::string qruleTypeName = (params->d_db)->getStringWithDefault("QRULE_TYPE", "QGAUSS");
    libMeshEnums::QuadratureType qruleType = Utility::string_to_enum<libMeshEnums::QuadratureType>(qruleTypeName);

    d_feType.reset( new ::FEType(feTypeOrder, feFamily) );

    d_fe.reset( (::FEBase::build(dimension, (*d_feType))).release() ); 

    std::string qruleOrderName = (params->d_db)->getStringWithDefault("QRULE_ORDER", "DEFAULT");

    libMeshEnums::Order qruleOrder;

    if(qruleOrderName == "DEFAULT") {
      qruleOrder = d_feType->default_quadrature_order();
    } else {
      qruleOrder = Utility::string_to_enum<libMeshEnums::Order>(qruleOrderName);
    }

    d_qrule.reset( (::QBase::build(qruleType, dimension, qruleOrder)).release() );

    d_fe->attach_quadrature_rule( (d_qrule).get() );

  }

}
}

