#include "SourceNonlinearElement.h"
#include "utils/Utilities.h"
#include "utils/ProfilerApp.h"

/* Libmesh files */
#include "enum_order.h"
#include "enum_fe_family.h"
#include "enum_quadrature_type.h"
#include "auto_ptr.h"
#include "string_to_enum.h"

#include <string>

namespace AMP {
namespace Operator {

SourceNonlinearElement :: SourceNonlinearElement (const boost::shared_ptr<ElementOperationParameters> & params) :
    ElementOperation(params),
    d_elementOutputVector(NULL),
    d_elem(NULL)
{

    AMP_INSIST( (params.get() != NULL), "''params'' is NULL");

    AMP_INSIST( (((params->d_db).get()) != NULL), "NULL database" );

    std::string feTypeOrderName = (params->d_db)->getStringWithDefault("FE_ORDER", "FIRST");
    libMeshEnums::Order feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>(feTypeOrderName);

    std::string feFamilyName = (params->d_db)->getStringWithDefault("FE_FAMILY", "LAGRANGE");
    libMeshEnums::FEFamily feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>(feFamilyName);

    std::string qruleTypeName = (params->d_db)->getStringWithDefault("QRULE_TYPE", "QGAUSS");
    libMeshEnums::QuadratureType qruleType = Utility::string_to_enum<libMeshEnums::QuadratureType>(qruleTypeName);

    d_integrateVolume= (params->d_db)->getBoolWithDefault("INTEGRATEVOLUME", "TRUE");

    const unsigned int dimension = 3;

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

    d_fe->attach_quadrature_rule( d_qrule.get() );

    //        d_isInputType = (params->d_db)->getStringWithDefault("InputVariableType", "NodalScalar");

    d_JxW = &(d_fe->get_JxW());

    d_phi = &(d_fe->get_phi());

    d_dphi = &(d_fe->get_dphi());

}


void SourceNonlinearElement :: initializeForCurrentElement( const ::Elem* elem, const boost::shared_ptr<SourcePhysicsModel> & sourcePhysicsModel)
{
    d_elem = elem;

    if(sourcePhysicsModel.get() != NULL){
        d_sourcePhysicsModel = sourcePhysicsModel;}
}


void SourceNonlinearElement::apply() {

    PROFILE_START("apply",5);

    d_fe->reinit(d_elem);
    const unsigned int n_nodes = d_elem->n_nodes();
    const unsigned int n_points = d_qrule->n_points();
    const std::vector<Real> & JxW = (*d_JxW);
    const std::vector<std::vector<Real> > & phi = (*d_phi);
    // const std::vector<std::vector<RealGradient> > & dphi = (*d_dphi);
    // std::vector<std::vector<double> > & elementInputVector = d_elementInputVector;
    std::vector<double> & elementOutputVector = (*d_elementOutputVector);
    std::vector<double> source_physics(static_cast<int>(n_points));
    std::vector<std::vector<double> > source_vectors(d_elementInputVector.size());
    std::vector<std::vector<double> > auxillary_vectors(d_elementAuxVector.size());
    const std::vector<Point>& coordinates = d_fe->get_xyz();

    if (d_isInputType == "IntegrationPointScalar") {

        for (unsigned int var = 0; var < d_elementInputVector.size(); var++)
        {
            source_vectors[var] = d_elementInputVector[var];
            if(d_elementAuxVector.size()>0){    auxillary_vectors[var] = d_elementAuxVector[var];}
        }

        if(d_sourcePhysicsModel.get() != NULL){
            d_sourcePhysicsModel->getConstitutiveProperty(source_physics, source_vectors, auxillary_vectors, coordinates);}

    }else if(d_isInputType== "NodalScalar"){

        for (unsigned int var = 0; var < d_elementInputVector.size(); var++)
        {
            source_vectors[var].resize(n_points);
            for (unsigned int qp = 0; qp < n_points; qp++)
            {
                source_vectors[var][qp] = 0.0;
                for (unsigned int j = 0; j < n_nodes; j++)
                {
                    source_vectors[var][qp] += d_elementInputVector[var][j] * phi[j][qp];
                }//end for j
            }//end for qp
            if(d_elementAuxVector.size()>0){
                auxillary_vectors[var].resize(n_points);
                for (unsigned int qp = 0; qp < n_points; qp++)
                {
                    auxillary_vectors[var][qp] = 0.0;
                    for (unsigned int j = 0; j < n_nodes; j++)
                    {
                        auxillary_vectors[var][qp]    += d_elementAuxVector[var][j] * phi[j][qp];
                    }//end for j
                }
            }
        }

        if(d_sourcePhysicsModel.get() != NULL){
            d_sourcePhysicsModel->getConstitutiveProperty(source_physics, source_vectors, auxillary_vectors, coordinates);}

    } 
    
    if(d_sourcePhysicsModel.get() == NULL) {
        AMP_INSIST( (source_vectors.size() == 1), "In absence of SourcePhysicsModel Element Operation Expects only one source vector" );
    }
    for (unsigned int j = 0; j < n_nodes; j++) {
        for (unsigned int qp = 0; qp < n_points; qp++) {
            if(d_sourcePhysicsModel.get() != NULL){
                elementOutputVector[j] += (JxW[qp] * source_physics[qp] * phi[j][qp] );
            }else{
                if(d_integrateVolume){
                    elementOutputVector[j] += (JxW[qp] * source_vectors[0][qp] * phi[j][qp] );
                }else{
                    elementOutputVector[j] += ( source_vectors[0][qp] * phi[j][qp] )/8;
                }
            }
        }//end for j
    }//end for qp

    PROFILE_STOP("apply",5);
}


}
}

