#include "NodeToGaussPointOperator.h"
#include "utils/ProfilerApp.h"

/* Libmesh files */
#include "fe_type.h"
#include "fe_base.h"
#include "elem.h"
#include "quadrature.h"
#include "enum_order.h"
#include "enum_fe_family.h"
#include "enum_quadrature_type.h"
#include "auto_ptr.h"
#include "string_to_enum.h"


namespace AMP {
namespace Operator {


// Constructor
NodeToGaussPointOperator::NodeToGaussPointOperator (const boost::shared_ptr<OperatorParameters> & params) : Operator (params)
{
    d_NodalVariable.reset(new AMP::LinearAlgebra::Variable(params->d_db->getString("InputVariable")));
    d_GaussPtVariable.reset(new AMP::LinearAlgebra::Variable(params->d_db->getString("OutputVariable")));
    d_UseSurfaceElements = (params->d_db)->getBoolWithDefault("UseSurfaceElements", true);
    // Create a list of the libmesh elements
    if(d_UseSurfaceElements) {
        d_iterator = d_Mesh->getIterator(AMP::Mesh::Face, 0);
    } else {
        d_iterator = d_Mesh->getIterator(AMP::Mesh::Volume, 0);
    }
    d_libmeshElements.reinit( d_iterator );
}


// Apply operator
void NodeToGaussPointOperator::apply(AMP::LinearAlgebra::Vector::const_shared_ptr,
        AMP::LinearAlgebra::Vector::const_shared_ptr u,
        AMP::LinearAlgebra::Vector::shared_ptr r, const double , const double ) 
{ 
    PROFILE_START("apply");

    AMP::LinearAlgebra::Vector::const_shared_ptr nodalVec = subsetInputVector(u);
    AMP::LinearAlgebra::Vector::shared_ptr gaussPtVec = subsetOutputVector(r);

    AMP_ASSERT(nodalVec->getUpdateStatus()==AMP::LinearAlgebra::Vector::UNCHANGED);

    AMP::Discretization::DOFManager::shared_ptr dof_map = nodalVec->getDOFManager();
    AMP::Discretization::DOFManager::shared_ptr gaussPt_dof_map = gaussPtVec->getDOFManager();

    int dim=0;
    if(d_UseSurfaceElements)
        dim = 2;
    else
        dim = 3;
    libMeshEnums::Order feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>("FIRST");
    libMeshEnums::FEFamily feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>("LAGRANGE");
    libMeshEnums::Order qruleOrder = Utility::string_to_enum<libMeshEnums::Order>("SECOND");
    boost::shared_ptr < ::QBase > d_qrule ( (::QBase::build("QGAUSS", dim, qruleOrder)).release() );
    boost::shared_ptr < ::FEType > d_feType ( new ::FEType(feTypeOrder, feFamily) );

    AMP::Mesh::MeshIterator iterator = d_iterator.begin();
    for (size_t i=0; i<iterator.size(); i++) {

        // Get the dofs for the gauss points
        std::vector<size_t> gaussPtIndices; 
        gaussPt_dof_map->getDOFs( iterator->globalID(), gaussPtIndices );

        // Check if we need to set any gauss points for the current element
        if ( gaussPtIndices.size()==0 ) {
            ++iterator;
            continue;
        }

        // Get the nodes comprising the current element
        std::vector<AMP::Mesh::MeshElement> d_currNodes = iterator->getElements(AMP::Mesh::Vertex);  

        // Get the current libmesh element
        ::Elem* currElem = d_libmeshElements.getElement( iterator->globalID() );
        boost::shared_ptr< ::FEBase > d_fe( (::FEBase::build(dim, (*d_feType))).release() );
        d_fe->attach_quadrature_rule( d_qrule.get() );
        d_fe->reinit(currElem);

        // Get the dofs for the nodes
        std::vector<size_t> bndGlobalIds;
        std::vector<AMP::Mesh::MeshElementID> globalIds(d_currNodes.size()) ;
        for (size_t j=0; j<d_currNodes.size(); j++)
            globalIds[j] = d_currNodes[j].globalID();
        dof_map->getDOFs(globalIds, bndGlobalIds);

        // Get the values at the gauss points
        const std::vector<std::vector<Real> >& phi = d_fe->get_phi();
        AMP_ASSERT(gaussPtIndices.size()==d_qrule->n_points());
        std::vector<double>  computedAtGauss(d_qrule->n_points(), 0.0);
        for (unsigned int j = 0; j < bndGlobalIds.size(); ++j){
            double computedAtNode = nodalVec->getValueByGlobalID(bndGlobalIds[j]); 
            for (unsigned int qp = 0; qp < d_qrule->n_points(); ++qp){
                computedAtGauss[qp] += ( computedAtNode * phi[j][qp] );
            }//end for qp
        }//end for j

        gaussPtVec->setLocalValuesByGlobalID( gaussPtIndices.size(),&gaussPtIndices[0], &computedAtGauss[0]);

        ++iterator;
    }//end for
    PROFILE_STOP("apply");
}// end apply


}
}



