
#include "FlowElement.h"
#include "utils/Utilities.h"

/* Libmesh files */
#include "libmesh/auto_ptr.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/string_to_enum.h"

#include <string>

namespace AMP {
namespace Operator {

FlowElement::FlowElement( const AMP::shared_ptr<ElementOperationParameters> &params )
    : ElementOperation( params ), d_elem( nullptr )
{
    AMP_INSIST( ( params.get() != nullptr ), "''params'' is NULL" );

    AMP_INSIST( ( ( ( params->d_db ).get() ) != nullptr ), "NULL database" );

    const unsigned int dimension = ( params->d_db )->getIntegerWithDefault( "DIMENSION", 3 );
    // int numApprox = (params->d_db)->getIntegerWithDefault("NUM_APPROX", 2);

    std::string U_feTypeOrderName = ( params->d_db )->getStringWithDefault( "FE_ORDER", "SECOND" );
    auto feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>( U_feTypeOrderName );

    std::string feFamilyName = ( params->d_db )->getStringWithDefault( "FE_FAMILY", "LAGRANGE" );
    auto feFamily            = Utility::string_to_enum<libMeshEnums::FEFamily>( feFamilyName );

    std::string qruleTypeName = ( params->d_db )->getStringWithDefault( "QRULE_TYPE", "QGAUSS" );
    auto qruleType = Utility::string_to_enum<libMeshEnums::QuadratureType>( qruleTypeName );

    d_feType.reset( new ::FEType( feTypeOrder, feFamily ) );

    d_fe.reset( (::FEBase::build( dimension, ( *d_feType ) ) ).release() );

    std::string qruleOrderName = ( params->d_db )->getStringWithDefault( "QRULE_ORDER", "DEFAULT" );

    libMeshEnums::Order qruleOrder;

    if ( qruleOrderName == "DEFAULT" ) {
        qruleOrder = d_feType->default_quadrature_order();
    } else {
        qruleOrder = Utility::string_to_enum<libMeshEnums::Order>( qruleOrderName );
    }

    d_qrule.reset( (::QBase::build( qruleType, dimension, qruleOrder ) ).release() );

    d_fe->attach_quadrature_rule( ( d_qrule ).get() );
}
} // namespace Operator
} // namespace AMP
