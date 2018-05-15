#include "DiffusionElement.h"
#include "AMP/utils/Utilities.h"

// Libmesh headers
DISABLE_WARNINGS
#include "libmesh/auto_ptr.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/string_to_enum.h"
ENABLE_WARNINGS

#include <string>

namespace AMP {
namespace Operator {


DiffusionElement::DiffusionElement( const AMP::shared_ptr<ElementOperationParameters> &params )
    : ElementOperation( params ),
      d_JxW( nullptr ),
      d_phi( nullptr ),
      d_dphi( nullptr ),
      d_elem( nullptr )
{

    AMP_INSIST( ( params.get() != nullptr ), "''params'' is NULL" );

    AMP_INSIST( ( ( ( params->d_db ).get() ) != nullptr ), "NULL database" );

    std::string feTypeOrderName = ( params->d_db )->getStringWithDefault( "FE_ORDER", "FIRST" );
    auto feTypeOrder            = Utility::string_to_enum<libMeshEnums::Order>( feTypeOrderName );

    std::string feFamilyName = ( params->d_db )->getStringWithDefault( "FE_FAMILY", "LAGRANGE" );
    auto feFamily            = Utility::string_to_enum<libMeshEnums::FEFamily>( feFamilyName );

    std::string qruleTypeName = ( params->d_db )->getStringWithDefault( "QRULE_TYPE", "QGAUSS" );
    auto qruleType = Utility::string_to_enum<libMeshEnums::QuadratureType>( qruleTypeName );

    const unsigned int dimension = 3;

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

    d_fe->attach_quadrature_rule( d_qrule.get() );

    d_JxW = &( d_fe->get_JxW() );

    d_phi = &( d_fe->get_phi() );

    d_dphi = &( d_fe->get_dphi() );
}


void DiffusionElement::initializeForCurrentElement(
    const ::Elem *elem, const AMP::shared_ptr<DiffusionTransportModel> &transportModel )
{
    d_elem           = elem;
    d_transportModel = transportModel;
}
} // namespace Operator
} // namespace AMP
