#include "ConstructLinearMechanicsRHSVector.h"

// Libmesh headers
DISABLE_WARNINGS
#include "libmesh/auto_ptr.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/elem.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_type.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
ENABLE_WARNINGS


void computeTemperatureRhsVector(
    AMP::Mesh::Mesh::shared_ptr mesh,
    AMP::shared_ptr<AMP::Database> input_db,
    AMP::LinearAlgebra::Variable::shared_ptr,
    AMP::LinearAlgebra::Variable::shared_ptr displacementVar,
    const AMP::shared_ptr<AMP::LinearAlgebra::Vector> &currTemperatureVec,
    const AMP::shared_ptr<AMP::LinearAlgebra::Vector> &prevTemperatureVec,
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec )
{
    currTemperatureVec->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
    prevTemperatureVec->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );

    AMP::LinearAlgebra::Vector::shared_ptr rInternal =
        rhsVec->subsetVectorForVariable( displacementVar );
    rInternal->zero();

    AMP::shared_ptr<AMP::Database> elementRhsDatabase = input_db->getDatabase( "RhsElements" );
    AMP::shared_ptr<AMP::Database> materialModelDatabase =
        input_db->getDatabase( "RhsMaterialModel" );

    AMP::shared_ptr<::FEType> feType;
    AMP::shared_ptr<::FEBase> fe;
    AMP::shared_ptr<::QBase> qrule;

    std::string feTypeOrderName = elementRhsDatabase->getStringWithDefault( "FE_ORDER", "FIRST" );
    auto feTypeOrder            = Utility::string_to_enum<libMeshEnums::Order>( feTypeOrderName );

    std::string feFamilyName = elementRhsDatabase->getStringWithDefault( "FE_FAMILY", "LAGRANGE" );
    auto feFamily            = Utility::string_to_enum<libMeshEnums::FEFamily>( feFamilyName );

    std::string qruleTypeName = elementRhsDatabase->getStringWithDefault( "QRULE_TYPE", "QGAUSS" );
    auto qruleType = Utility::string_to_enum<libMeshEnums::QuadratureType>( qruleTypeName );

    const unsigned int dimension = 3;

    feType.reset( new ::FEType( feTypeOrder, feFamily ) );
    fe.reset( (::FEBase::build( dimension, ( *feType ) ) ).release() );

    std::string qruleOrderName =
        elementRhsDatabase->getStringWithDefault( "QRULE_ORDER", "DEFAULT" );
    libMeshEnums::Order qruleOrder;
    if ( qruleOrderName == "DEFAULT" ) {
        qruleOrder = feType->default_quadrature_order();
    } else {
        qruleOrder = Utility::string_to_enum<libMeshEnums::Order>( qruleOrderName );
    }

    qrule.reset( (::QBase::build( qruleType, dimension, qruleOrder ) ).release() );
    fe->attach_quadrature_rule( qrule.get() );

    const std::vector<Real> &JxW                       = ( fe->get_JxW() );
    const std::vector<std::vector<RealGradient>> &dphi = ( fe->get_dphi() );
    const std::vector<std::vector<Real>> &phi          = ( fe->get_phi() );

    AMP::shared_ptr<AMP::Materials::Material> material;
    double youngsModulus = 1.0e10, poissonsRatio = 0.33, thermalExpansionCoefficient = 2.0e-6;
    double default_BURNUP, default_OXYGEN_CONCENTRATION;

    bool useMaterialsLibrary =
        materialModelDatabase->getBoolWithDefault( "USE_MATERIALS_LIBRARY", false );
    if ( useMaterialsLibrary == true ) {
        AMP_INSIST( ( materialModelDatabase->keyExists( "Material" ) ),
                    "Key ''Material'' is missing!" );
        std::string matname = materialModelDatabase->getString( "Material" );
        material = AMP::voodoo::Factory<AMP::Materials::Material>::instance().create( matname );
    }

    if ( useMaterialsLibrary == false ) {
        AMP_INSIST( materialModelDatabase->keyExists( "THERMAL_EXPANSION_COEFFICIENT" ),
                    "Missing key: THERMAL_EXPANSION_COEFFICIENT" );
        AMP_INSIST( materialModelDatabase->keyExists( "Youngs_Modulus" ),
                    "Missing key: Youngs_Modulus" );
        AMP_INSIST( materialModelDatabase->keyExists( "Poissons_Ratio" ),
                    "Missing key: Poissons_Ratio" );

        thermalExpansionCoefficient =
            materialModelDatabase->getDouble( "THERMAL_EXPANSION_COEFFICIENT" );
        youngsModulus = materialModelDatabase->getDouble( "Youngs_Modulus" );
        poissonsRatio = materialModelDatabase->getDouble( "Poissons_Ratio" );
    }

    default_BURNUP = materialModelDatabase->getDoubleWithDefault( "Default_Burnup", 0.0 );
    default_OXYGEN_CONCENTRATION =
        materialModelDatabase->getDoubleWithDefault( "Default_Oxygen_Concentration", 0.0 );

    AMP::Discretization::DOFManager::shared_ptr dof_map_0 = rInternal->getDOFManager();
    AMP::Discretization::DOFManager::shared_ptr dof_map_1 = currTemperatureVec->getDOFManager();

    AMP::Mesh::MeshIterator el     = mesh->getIterator( AMP::Mesh::GeomType::Volume, 0 );
    AMP::Mesh::MeshIterator end_el = el.end();

    for ( ; el != end_el; ++el ) {
        std::vector<AMP::Mesh::MeshElement> currNodes =
            el->getElements( AMP::Mesh::GeomType::Vertex );
        size_t numNodesInCurrElem = currNodes.size();

        std::vector<std::vector<size_t>> type0DofIndices( currNodes.size() );
        std::vector<std::vector<size_t>> type1DofIndices( currNodes.size() );
        for ( size_t j = 0; j < currNodes.size(); ++j ) {
            dof_map_0->getDOFs( currNodes[j].globalID(), type0DofIndices[j] );
            dof_map_1->getDOFs( currNodes[j].globalID(), type1DofIndices[j] );
        } // end j

        std::vector<double> elementForceVector( ( 3 * numNodesInCurrElem ), 0.0 );
        std::vector<double> currElementTemperatureVector( numNodesInCurrElem );
        std::vector<double> prevElementTemperatureVector( numNodesInCurrElem );

        for ( size_t r = 0; r < numNodesInCurrElem; ++r ) {
            currElementTemperatureVector[r] =
                currTemperatureVec->getValueByGlobalID( type1DofIndices[r][0] );
            prevElementTemperatureVector[r] =
                prevTemperatureVec->getValueByGlobalID( type1DofIndices[r][0] );
        } // end r

        ::Elem *elem = new ::Hex8;
        for ( size_t j = 0; j < currNodes.size(); ++j ) {
            std::vector<double> pt = currNodes[j].coord();
            elem->set_node( j )    = new ::Node( pt[0], pt[1], pt[2], j );
        } // end j

        fe->reinit( elem );

        for ( unsigned int qp = 0; qp < qrule->n_points(); ++qp ) {
            double Bl_np1[6][24];

            for ( auto &elem : Bl_np1 ) {
                for ( size_t j = 0; j < ( 3 * numNodesInCurrElem ); ++j ) {
                    elem[j] = 0.0;
                } // end j
            }     // end i

            for ( size_t i = 0; i < numNodesInCurrElem; ++i ) {
                Bl_np1[0][( 3 * i ) + 0] = dphi[i][qp]( 0 );
                Bl_np1[1][( 3 * i ) + 1] = dphi[i][qp]( 1 );
                Bl_np1[2][( 3 * i ) + 2] = dphi[i][qp]( 2 );
                Bl_np1[3][( 3 * i ) + 1] = dphi[i][qp]( 2 );
                Bl_np1[3][( 3 * i ) + 2] = dphi[i][qp]( 1 );
                Bl_np1[4][( 3 * i ) + 0] = dphi[i][qp]( 2 );
                Bl_np1[4][( 3 * i ) + 2] = dphi[i][qp]( 0 );
                Bl_np1[5][( 3 * i ) + 0] = dphi[i][qp]( 1 );
                Bl_np1[5][( 3 * i ) + 1] = dphi[i][qp]( 0 );
            } // end i

            double currTemperatureAtGaussPoint = 0.0;
            double prevTemperatureAtGaussPoint = 0.0;
            for ( unsigned int k = 0; k < numNodesInCurrElem; k++ ) {
                currTemperatureAtGaussPoint += ( currElementTemperatureVector[k] * phi[k][qp] );
                prevTemperatureAtGaussPoint += ( prevElementTemperatureVector[k] * phi[k][qp] );
            } // end k

            if ( useMaterialsLibrary == true ) {
                std::map<std::string, AMP::shared_ptr<std::vector<double>>> inputMaterialParameters;

                std::string temperatureString = "temperature"; // in the future get from input file
                std::string burnupString      = "burnup";      // in the future get from input file
                std::string oxygenString = "concentration";    // in the future get from input file

                AMP::shared_ptr<std::vector<double>> tempVec( new std::vector<double> );
                AMP::shared_ptr<std::vector<double>> burnupVec( new std::vector<double> );
                AMP::shared_ptr<std::vector<double>> oxygenVec( new std::vector<double> );

                inputMaterialParameters.insert( std::make_pair( temperatureString, tempVec ) );
                inputMaterialParameters.insert( std::make_pair( burnupString, burnupVec ) );
                inputMaterialParameters.insert( std::make_pair( oxygenString, oxygenVec ) );

                tempVec->push_back( currTemperatureAtGaussPoint );
                burnupVec->push_back( default_BURNUP );
                oxygenVec->push_back( default_OXYGEN_CONCENTRATION );

                std::vector<double> YM( 1 );
                std::vector<double> PR( 1 );
                std::vector<double> TEC( 1 );

                std::string ymString  = "YoungsModulus";
                std::string prString  = "PoissonRatio";
                std::string tecString = "ThermalExpansion";

                material->property( ymString )->evalv( YM, inputMaterialParameters );
                material->property( prString )->evalv( PR, inputMaterialParameters );
                material->property( tecString )->evalv( TEC, inputMaterialParameters );

                youngsModulus               = YM[0];
                poissonsRatio               = PR[0];
                thermalExpansionCoefficient = TEC[0];
            }

            double d_thermalStress[6], d_thermalStrain[6], d_constitutiveMatrix[6][6];
            for ( unsigned int i = 0; i < 6; i++ ) {
                d_thermalStrain[i] = 0.0;
                d_thermalStress[i] = 0.0;
                for ( unsigned int j = 0; j < 6; j++ )
                    d_constitutiveMatrix[i][j] = 0.0;
            }

            double E  = youngsModulus;
            double nu = poissonsRatio;
            double K  = E / ( 3.0 * ( 1.0 - ( 2.0 * nu ) ) );
            double G  = E / ( 2.0 * ( 1.0 + nu ) );

            for ( unsigned int i = 0; i < 3; i++ )
                d_constitutiveMatrix[i][i] += ( 2.0 * G );

            for ( unsigned int i = 3; i < 6; i++ )
                d_constitutiveMatrix[i][i] += G;

            for ( unsigned int i = 0; i < 3; i++ ) {
                for ( unsigned int j = 0; j < 3; j++ ) {
                    d_constitutiveMatrix[i][j] += ( K - ( ( 2.0 * G ) / 3.0 ) );
                }
            }

            for ( unsigned int i = 0; i < 3; i++ ) {
                d_thermalStrain[i] = thermalExpansionCoefficient *
                                     ( currTemperatureAtGaussPoint - prevTemperatureAtGaussPoint );
            }

            for ( unsigned int i = 0; i < 6; i++ ) {
                for ( unsigned int j = 0; j < 6; j++ ) {
                    d_thermalStress[i] += ( d_constitutiveMatrix[i][j] * d_thermalStrain[j] );
                }
            }

            for ( unsigned int j = 0; j < numNodesInCurrElem; j++ ) {
                for ( unsigned int d = 0; d < 3; d++ ) {

                    double tmp = 0;
                    for ( unsigned int m = 0; m < 6; m++ ) {
                        tmp += ( Bl_np1[m][( 3 * j ) + d] * d_thermalStress[m] );
                    }

                    elementForceVector[( 3 * j ) + d] += ( JxW[qp] * tmp );
                } // end d
            }     // end j
        }         // end qp

        for ( unsigned int r = 0; r < numNodesInCurrElem; r++ ) {
            for ( unsigned int d = 0; d < 3; d++ ) {
                rInternal->addValueByGlobalID( type0DofIndices[r][d],
                                               elementForceVector[( 3 * r ) + d] );
            } // end d
        }     // end r

        for ( size_t j = 0; j < elem->n_nodes(); ++j ) {
            delete ( elem->get_node( j ) );
            elem->set_node( j ) = nullptr;
        } // end j
        delete elem;
        elem = nullptr;
    } // end el

    rInternal->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_ADD );
}
