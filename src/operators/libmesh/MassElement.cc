#include "MassElement.h"
#include "AMP/utils/Utilities.h"

/* Libmesh files */
#include "libmesh/auto_ptr.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/string_to_enum.h"

#include <string>

namespace AMP {
namespace Operator {

MassElement::MassElement( const AMP::shared_ptr<ElementOperationParameters> &params )
    : ElementOperation( params ), d_elem( nullptr )
{

    AMP_INSIST( ( params.get() != nullptr ), "''params'' is NULL" );

    AMP_INSIST( ( ( ( params->d_db ).get() ) != nullptr ), "NULL database" );

    std::string feTypeOrderName = ( params->d_db )->getStringWithDefault( "FE_ORDER", "FIRST" );

    auto feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>( feTypeOrderName );

    std::string feFamilyName = ( params->d_db )->getStringWithDefault( "FE_FAMILY", "LAGRANGE" );

    auto feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>( feFamilyName );

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

    // d_dphi = &(d_fe->get_dphi());
    d_phi = &( d_fe->get_phi() );
}


void MassElement::initializeForCurrentElement(
    const ::Elem *elem, const AMP::shared_ptr<MassDensityModel> &densityModel )
{
    d_elem = elem;

    d_densityModel = densityModel;
}
} // namespace Operator
} // namespace AMP
