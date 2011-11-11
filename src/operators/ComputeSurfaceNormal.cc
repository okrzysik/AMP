
#include "ComputeSurfaceNormal.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"

namespace AMP {
namespace Operator {

  ComputeSurfaceNormal :: ComputeSurfaceNormal(const boost::shared_ptr<OperatorParameters> & params)
  {
    d_MeshAdapter = params->d_MeshAdapter;

    std::string feTypeOrderName = (params->d_db)->getStringWithDefault("FE_ORDER", "FIRST");

    libMeshEnums::Order feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>(feTypeOrderName);

    std::string feFamilyName = (params->d_db)->getStringWithDefault("FE_FAMILY", "LAGRANGE");

    libMeshEnums::FEFamily feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>(feFamilyName);

    std::string qruleTypeName = (params->d_db)->getStringWithDefault("QRULE_TYPE", "QGAUSS");

    libMeshEnums::QuadratureType qruleType = Utility::string_to_enum<libMeshEnums::QuadratureType>(qruleTypeName);

    d_dimension = (params->d_db)->getIntegerWithDefault("DIMENSION", 3);

    d_feType.reset( new ::FEType(feTypeOrder, feFamily) );

    d_fe_face.reset( (::FEBase::build((d_dimension - 1), (*d_feType))).release() );

    d_fe.reset( (::FEBase::build(d_dimension, (*d_feType))).release() );

    std::string qruleOrderName = (params->d_db)->getStringWithDefault("QRULE_ORDER", "DEFAULT");

    libMeshEnums::Order qruleOrder;

    if(qruleOrderName == "DEFAULT") {
      qruleOrder = d_feType->default_quadrature_order();
    } else {
      qruleOrder = Utility::string_to_enum<libMeshEnums::Order>(qruleOrderName);
    }

    d_qrule.reset( (::QBase::build(qruleType, (d_dimension - 1), qruleOrder)).release() );

    d_fe_face->attach_quadrature_rule( d_qrule.get() );

    d_fe->attach_quadrature_rule( d_qrule.get() );

    d_JxW = &(d_fe->get_JxW());

    d_boundaryIds.resize(d_numIds);

    char key[100];
    for(int j = 0; j < d_numIds; j++) {
      sprintf(key, "id_%d", j);
      AMP_INSIST( (params->d_db)->keyExists(key), "Key is missing!" );
      d_boundaryIds[j] = (params->d_db)->getInteger(key);
    }//end for j
  }

  boost::shared_ptr<AMP::LinearAlgebra::Vector> ComputeSurfaceNormal :: getNormals(const AMP::LinearAlgebra::Vector::shared_ptr &u)
  {

    d_inputVec = u->subsetVectorForVariable(d_variable);
   
    d_outputVec = d_inputVec->cloneVector();

    d_inputVec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );

    AMP::Mesh::DOFMap::shared_ptr dof_map = d_MeshAdapter->getDOFMap(d_variable);

    AMP::LinearAlgebra::Vector::shared_ptr rInternal = d_MeshAdapter->createVector(d_variable);

    rInternal->zero();

    unsigned int numIds = d_boundaryIds.size();

    for(unsigned int j = 0; j < numIds; j++) {
      AMP::Mesh::MeshManager::Adapter::BoundarySideIterator bnd = d_MeshAdapter->beginSideBoundary( d_boundaryIds[j] );
      AMP::Mesh::MeshManager::Adapter::BoundarySideIterator end_bnd = d_MeshAdapter->endSideBoundary( d_boundaryIds[j] );

      for( ; bnd != end_bnd; ++bnd) {

        AMP::Mesh::MeshManager::Adapter::Element cur_elem = d_MeshAdapter->getElementFromSide(*bnd);
        AMP::Mesh::MeshManager::Adapter::Element cur_side = *bnd;

        d_phi = &(d_fe_face->get_phi());
        d_JxW = &(d_fe_face->get_JxW());
        d_normal = &(d_fe->get_normals());

        d_face = &cur_side.getElem();
        d_elem = &cur_elem.getElem();

        unsigned int k = 0;
        for(unsigned int s = 0; s < d_elem->n_sides(); s++) {
          AutoPtr< ::Elem> side (d_elem->build_side(s));
          if(*(side.get()) == *d_face) {
            k = s;
          }
        }

        d_fe_face->reinit ( d_face );
        d_fe->reinit( d_elem, k );

        const std::vector<std::vector<Real> > & phi = *d_phi;
        const std::vector<Real> & djxw = *d_JxW;
        const std::vector<Point> & normal = *d_normal;

        for(int dofId = 0; dofId < 3; dofId++) {

          std::vector<unsigned int> bndGlobalIds;
          dof_map->getDOFs(cur_side, bndGlobalIds, dofId);

          std::vector<double> outputValue(bndGlobalIds.size(), 0.0);

          for(unsigned int i = 0; i < bndGlobalIds.size(); i++) {

            double inputValue = d_inputVec->getValueByGlobalID(bndGlobalIds[i]); 

            for(unsigned int qp = 0; qp < d_qrule->n_points(); qp++) 
            {
                outputValue[i] +=  (djxw[qp]*phi[i][qp]*inputValue*normal[qp](dofId));
            }

          }//end for i

          rInternal->addValuesByGlobalID((int)bndGlobalIds.size() , (int *)&(bndGlobalIds[0]), &(outputValue[0]));

        }//end for dofId

      }//end for bnd

    }//end for j

    d_outputVec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_ADD );

    return d_outputVec;
  }

}
}

