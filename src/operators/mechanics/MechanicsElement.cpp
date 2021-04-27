#include "MechanicsElement.h"
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

MechanicsElement::MechanicsElement( std::shared_ptr<const ElementOperationParameters> params )
    : ElementOperation( params ), d_elem( nullptr )
{
    AMP_INSIST( ( params ), "''params'' is NULL" );

    AMP_INSIST( ( ( ( params->d_db ).get() ) != nullptr ), "NULL database" );

    d_useReducedIntegration = params->d_db->getWithDefault( "USE_REDUCED_INTEGRATION", false );

    d_useJaumannRate = params->d_db->getWithDefault( "USE_JAUMANN_RATE", false );

    d_useFlanaganTaylorElem =
        params->d_db->getWithDefault( "USE_FLANAGAN_TAYLOR_ELEMENT_FORMULATION", false );

    if ( d_useFlanaganTaylorElem == true ) {
        AMP_INSIST(
            ( d_useJaumannRate == false ),
            "Flanagan Taylor element formulation can only be used with Green-Naghdi stress rate." );
    }

    std::string feTypeOrderName = params->d_db->getWithDefault<std::string>( "FE_ORDER", "FIRST" );

    auto feTypeOrder = libMesh::Utility::string_to_enum<libMeshEnums::Order>( feTypeOrderName );

    std::string feFamilyName = params->d_db->getWithDefault<std::string>( "FE_FAMILY", "LAGRANGE" );

    auto feFamily = libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( feFamilyName );

    std::string qruleTypeName = params->d_db->getWithDefault<std::string>( "QRULE_TYPE", "QGAUSS" );

    auto qruleType =
        libMesh::Utility::string_to_enum<libMeshEnums::QuadratureType>( qruleTypeName );

    const unsigned int dimension = 3;

    d_feType.reset( new libMesh::FEType( feTypeOrder, feFamily ) );

    d_fe.reset( ( libMesh::FEBase::build( dimension, ( *d_feType ) ) ).release() );

    std::string qruleOrderName =
        params->d_db->getWithDefault<std::string>( "QRULE_ORDER", "DEFAULT" );

    libMeshEnums::Order qruleOrder;

    if ( qruleOrderName == "DEFAULT" ) {
        qruleOrder = d_feType->default_quadrature_order();
    } else {
        qruleOrder = libMesh::Utility::string_to_enum<libMeshEnums::Order>( qruleOrderName );
    }

    d_qrule.reset( ( libMesh::QBase::build( qruleType, dimension, qruleOrder ) ).release() );

    d_fe->attach_quadrature_rule( d_qrule.get() );

    d_iDebugPrintInfoLevel = params->d_db->getWithDefault( "print_info_level", 0 );
}
} // namespace Operator
} // namespace AMP
