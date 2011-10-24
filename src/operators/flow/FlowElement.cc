
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
    
    d_feType.resize(2);
    d_fe.resize(2);
    d_qrule.resize(2);

    boost::shared_ptr<AMP::Database> U_db = (params->d_db)->getDatabase("VELOCITY_APPROX");
    boost::shared_ptr<AMP::Database> P_db = (params->d_db)->getDatabase("PRESSURE_APPROX");

    std::string U_feTypeOrderName = U_db->getStringWithDefault("FE_ORDER", "SECOND");
    libMeshEnums::Order U_feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>(U_feTypeOrderName);

    std::string U_feFamilyName = U_db->getStringWithDefault("FE_FAMILY", "LAGRANGE");
    libMeshEnums::FEFamily U_feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>(U_feFamilyName);

    std::string U_qruleTypeName = U_db->getStringWithDefault("QRULE_TYPE", "QGAUSS");
    libMeshEnums::QuadratureType U_qruleType = Utility::string_to_enum<libMeshEnums::QuadratureType>(U_qruleTypeName);

    d_feType[0].reset( new ::FEType(U_feTypeOrder, U_feFamily) );

    d_fe[0].reset( (::FEBase::build(dimension, (*d_feType[0]))).release() ); 

    std::string U_qruleOrderName = (params->d_db)->getStringWithDefault("QRULE_ORDER", "DEFAULT");

    libMeshEnums::Order U_qruleOrder;

    if(U_qruleOrderName == "DEFAULT") {
      U_qruleOrder = d_feType[0]->default_quadrature_order();
    } else {
      U_qruleOrder = Utility::string_to_enum<libMeshEnums::Order>(U_qruleOrderName);
    }

    (d_qrule[0]).reset( (::QBase::build(U_qruleType, dimension, U_qruleOrder)).release() );

    (d_fe[0])->attach_quadrature_rule( (d_qrule[0]).get() );

    /////

    std::string P_feTypeOrderName = P_db->getStringWithDefault("FE_ORDER", "SECOND");
    libMeshEnums::Order P_feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>(P_feTypeOrderName);

    std::string P_feFamilyName = P_db->getStringWithDefault("FE_FAMILY", "LAGRANGE");
    libMeshEnums::FEFamily P_feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>(P_feFamilyName);

    std::string P_qruleTypeName = P_db->getStringWithDefault("QRULE_TYPE", "QGAUSS");
    libMeshEnums::QuadratureType P_qruleType = Utility::string_to_enum<libMeshEnums::QuadratureType>(P_qruleTypeName);

    d_feType[1].reset( new ::FEType(P_feTypeOrder, P_feFamily) );

    d_fe[1].reset( (::FEBase::build(dimension, (*d_feType[1]))).release() ); 

    std::string P_qruleOrderName = P_db->getStringWithDefault("QRULE_ORDER", "DEFAULT");

    libMeshEnums::Order P_qruleOrder;

    if(P_qruleOrderName == "DEFAULT") {
      P_qruleOrder = d_feType[1]->default_quadrature_order();
    } else {
      P_qruleOrder = Utility::string_to_enum<libMeshEnums::Order>(P_qruleOrderName);
    }

    (d_qrule[1]).reset( (::QBase::build(P_qruleType, dimension, P_qruleOrder)).release() );

    (d_fe[1])->attach_quadrature_rule( (d_qrule[0]).get() );


  }

}
}

