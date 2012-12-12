
#include "PressureBoundaryVectorCorrection.h"
#include "PressureBoundaryVectorCorrectionParameters.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"

/* Libmesh files */

#include "enum_order.h"
#include "enum_fe_family.h"
#include "enum_quadrature_type.h"
#include "auto_ptr.h"
#include "string_to_enum.h"

#include "face_quad4.h"
#include "cell_hex8.h"
#include "node.h"

#include <string>


namespace AMP {
  namespace Operator {

    PressureBoundaryVectorCorrection::PressureBoundaryVectorCorrection(const boost::shared_ptr<PressureBoundaryVectorCorrectionParameters> & params)
      : BoundaryOperator (params)
    {
      d_params = params;

      std::string feTypeOrderName = (params->d_db)->getStringWithDefault("FE_ORDER", "FIRST");

      libMeshEnums::Order feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>(feTypeOrderName);

      std::string feFamilyName = (params->d_db)->getStringWithDefault("FE_FAMILY", "LAGRANGE");

      libMeshEnums::FEFamily feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>(feFamilyName);

      std::string qruleTypeName = (params->d_db)->getStringWithDefault("QRULE_TYPE", "QGAUSS");

      libMeshEnums::QuadratureType qruleType = Utility::string_to_enum<libMeshEnums::QuadratureType>(qruleTypeName);

      d_feType.reset( new ::FEType(feTypeOrder, feFamily) );

      d_fe.reset( (::FEBase::build(2, (*d_feType))).release() );

      d_fe_3d.reset( (::FEBase::build(3, (*d_feType))).release() );

      std::string qruleOrderName = (params->d_db)->getStringWithDefault("QRULE_ORDER", "DEFAULT");

      libMeshEnums::Order qruleOrder;

      if(qruleOrderName == "DEFAULT") {
        qruleOrder = d_feType->default_quadrature_order();
      } else {
        qruleOrder = Utility::string_to_enum<libMeshEnums::Order>(qruleOrderName);
      }

      d_qrule.reset( (::QBase::build(qruleType, 2, qruleOrder)).release() );

      d_fe->attach_quadrature_rule( d_qrule.get() );

      d_fe_3d->attach_quadrature_rule( d_qrule.get() );

      d_JxW = &(d_fe->get_JxW());

      d_variable = params->d_variable;

      reset(params);
    }

    void PressureBoundaryVectorCorrection :: reset(const boost::shared_ptr<OperatorParameters>& params)
    {
      boost::shared_ptr<PressureBoundaryVectorCorrectionParameters> myparams = 
        boost::dynamic_pointer_cast<PressureBoundaryVectorCorrectionParameters>(params);

      AMP_INSIST( ((myparams.get()) != NULL), "NULL parameters" );
      AMP_INSIST( (((myparams->d_db).get()) != NULL), "NULL database" );

      AMP_INSIST( (myparams->d_db)->keyExists("number_of_ids"), "Key ''number_of_ids'' is missing!" );
      d_numIds = (myparams->d_db)->getInteger("number_of_ids");

      d_isConstantPressure = (myparams->d_db)->getBoolWithDefault("constant_pressure", true);
      AMP_INSIST((d_isConstantPressure ==  true), "Variable pressure has not been implemented yet.");

      d_boundaryIds.resize(d_numIds);
      d_pressureValues.resize(d_numIds);

      char key[100];
      for(int j = 0; j < d_numIds; j++) {
        sprintf(key, "id_%d", j);
        AMP_INSIST( (myparams->d_db)->keyExists(key), "Key is missing!" );
        d_boundaryIds[j] = (myparams->d_db)->getInteger(key);
        sprintf(key, "value_%d", j);
        AMP_INSIST( (myparams->d_db)->keyExists(key), "Key is missing!" );
        d_pressureValues[j] = (myparams->d_db)->getDouble(key);
      }//end for j
    }

    void PressureBoundaryVectorCorrection :: addRHScorrection(AMP::LinearAlgebra::Vector::shared_ptr rhsCorrection)
    {
      if(d_iDebugPrintInfoLevel > 8) {
        AMP::pout << "Entered the computeRHScorrection function." << std::endl;
      }

      AMP::LinearAlgebra::Vector::shared_ptr myRhs = mySubsetVector(rhsCorrection, d_variable);
      d_dofManager = myRhs->getDOFManager();

      unsigned int numIds = d_boundaryIds.size();

      for(unsigned int j = 0; j < numIds; ++j) {
        if(d_iDebugPrintInfoLevel > 8) {
          AMP::pout << "Entered the numIds loop. d_boundaryIds[" << j << "] = " << d_boundaryIds[j] << std::endl;
        }

        AMP::Mesh::MeshIterator belem     = d_Mesh->getBoundaryIDIterator( AMP::Mesh::Volume, d_boundaryIds[j], 0 );
        AMP::Mesh::MeshIterator end_belem = belem.end();

        for( ; belem  != end_belem; ++belem) {
          if(d_iDebugPrintInfoLevel > 8) {
            AMP::pout << "Entered the bnd loop." << std::endl;
          }

          d_currFaces = belem->getElements(AMP::Mesh::Face);
          d_currNodes = belem->getElements(AMP::Mesh::Vertex);
          createCurrentLibMeshElement();

          for(size_t bb = 0; bb < d_currFaces.size(); ++bb) {
            if(d_currFaces[bb].isOnBoundary(d_boundaryIds[j])) {
              if(d_iDebugPrintInfoLevel > 8) {
                AMP::pout << "The side that matched is = " << bb << std::endl;
              }

              d_currNodes = d_currFaces[bb].getElements(AMP::Mesh::Vertex);
              unsigned int numNodesInCurrSide = d_currNodes.size();

              createCurrentLibMeshSide();

              getDofIndicesForCurrentSide();

              d_phi = &(d_fe->get_phi());
              d_JxW = &(d_fe->get_JxW());
              d_normal = &(d_fe_3d->get_normals());

              d_fe->reinit ( d_currSidePtr );
              d_fe_3d->reinit( d_currElemPtr , bb );

              const std::vector<std::vector<Real> > & phi = *d_phi;
              const std::vector<Real> & djxw = *d_JxW;
              const std::vector<Point> & normal = *d_normal;

              for(int dofId = 0; dofId < 3; ++dofId) {
                if(d_iDebugPrintInfoLevel > 8) {
                  AMP::pout << "Entered the dofId loop." << std::endl;
                  AMP::pout << " d_qrule->n_points() =  " << d_qrule->n_points()<< std::endl; 
                  AMP::pout << " d_normal->size() = " << d_normal->size() << std::endl;
                  AMP::pout << "d_JxW->size() = " << d_JxW->size() << " d_phi->size() = " << d_phi->size() << std::endl;
                }

                std::vector<double> pressure(numNodesInCurrSide , 0.0);

                for(unsigned int i = 0; i < numNodesInCurrSide ; ++i) {
                  for(unsigned int qp = 0; qp < d_qrule->n_points(); ++qp) {
                    //std::vector<std::vector<double> > temp(1) ;
                    if(d_iDebugPrintInfoLevel > 8) {
                      AMP::pout << "The weights are as follows : "  << djxw[qp] << " shape function " <<
                        phi[i][qp] << " pressure " << d_pressureValues[j] << std::endl;
                      AMP::pout << "Direction of the normals are as follows : " << normal[qp](0) << " , " <<
                        normal[qp](1) << " , " << normal[qp](2) << std::endl;
                    }
                    if(dofId == 0) pressure[i] +=  (djxw[qp]*phi[i][qp]*d_pressureValues[j]*normal[qp](0));
                    if(dofId == 1) pressure[i] +=  (djxw[qp]*phi[i][qp]*d_pressureValues[j]*normal[qp](1));
                    if(dofId == 2) pressure[i] +=  (djxw[qp]*phi[i][qp]*d_pressureValues[j]*normal[qp](2));
                    if(d_iDebugPrintInfoLevel > 8) {
                      AMP::pout << "dofID = " << dofId << " pressure[" << i << "] = " << pressure[i] << std::endl;
                    }
                  }//end qp
                }//end i
                myRhs->addValuesByGlobalID(d_dofIndices.size(), &(d_dofIndices[0]), &(pressure[0]));
              }//end for dofId
              destroyCurrentLibMeshSide();
            }
          }//end bb
          destroyCurrentLibMeshElement();
        }//end for belem
      }//end for j
      ////////////////////////////////////////////////////////////////

      myRhs->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_ADD );
    }

    void PressureBoundaryVectorCorrection :: apply(AMP::LinearAlgebra::Vector::const_shared_ptr f,
        AMP::LinearAlgebra::Vector::const_shared_ptr u, AMP::LinearAlgebra::Vector::shared_ptr r, 
        const double a , const double b ) {
      AMP::LinearAlgebra::Vector::shared_ptr rInternal = mySubsetVector(r, d_variable);

      if(d_iDebugPrintInfoLevel>3)
      {
        AMP::pout << "L2 Norm of rInternal entering PressureBoundaryVectorCorrection::apply is : " << rInternal->L2Norm() << std::endl;
        AMP::pout << "The scaling value is : " << a << std::endl;
      }

      //      computeRHScorrection();
      //      d_rhsCorrectionAdd->scale(a);
      addRHScorrection(rInternal);

      if(d_iDebugPrintInfoLevel>3)
      {
        AMP::pout << "L2 Norm of rInternal exiting PressureBoundaryVectorCorrection::apply is : " << rInternal->L2Norm() << std::endl;
      }
    }

    boost::shared_ptr<OperatorParameters> PressureBoundaryVectorCorrection :: 
      getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& ) {
        boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));

        tmp_db->putString("FE_ORDER","FIRST");
        tmp_db->putString("FE_FAMILY","LAGRANGE");
        tmp_db->putString("QRULE_TYPE","QGAUSS");
        tmp_db->putInteger("DIMENSION",2);
        tmp_db->putString("QRULE_ORDER","DEFAULT");
        tmp_db->putInteger("number_of_ids",d_numIds);
        tmp_db->putBool("constant_pressure", d_isConstantPressure);

        char key[100];
        for(int j = 0; j < d_numIds; j++) {
          sprintf(key, "id_%d", j);
          tmp_db->putInteger(key,d_boundaryIds[j]);
          sprintf(key, "value_%d", j);
          tmp_db->putInteger(key,d_pressureValues[j]);
        }

        tmp_db->putBool("skip_params", true);

        boost::shared_ptr<PressureBoundaryVectorCorrectionParameters> outParams(new PressureBoundaryVectorCorrectionParameters(tmp_db));

        return outParams;
      }

    void PressureBoundaryVectorCorrection :: createCurrentLibMeshSide() {
      d_currSidePtr = new ::Quad4;
      AMP_ASSERT(d_currNodes.size()==4);
      for(size_t j = 0; j < d_currNodes.size(); j++) {
        std::vector<double> pt = d_currNodes[j].coord();
        d_currSidePtr->set_node(j) = new ::Node(pt[0], pt[1], pt[2], j);
      }//end for j
    }

    void PressureBoundaryVectorCorrection :: createCurrentLibMeshElement() {
      d_currElemPtr = new ::Hex8;
      AMP_ASSERT(d_currNodes.size()==8);
      for(size_t j = 0; j < d_currNodes.size(); j++) {
        std::vector<double> pt = d_currNodes[j].coord();
        d_currElemPtr->set_node(j) = new ::Node(pt[0], pt[1], pt[2], j);
      }//end for j
    }

    void PressureBoundaryVectorCorrection :: destroyCurrentLibMeshSide() {
      for(size_t j = 0; j < d_currSidePtr->n_nodes(); j++) {
        delete (d_currSidePtr->get_node(j));
        d_currSidePtr->set_node(j) = NULL;
      }//end for j
      delete d_currSidePtr;
      d_currSidePtr = NULL;
    }

    void PressureBoundaryVectorCorrection :: destroyCurrentLibMeshElement() {
      for(size_t j = 0; j < d_currElemPtr->n_nodes(); j++) {
        delete (d_currElemPtr->get_node(j));
        d_currElemPtr->set_node(j) = NULL;
      }//end for j
      delete d_currElemPtr;
      d_currElemPtr = NULL;
    }

    void PressureBoundaryVectorCorrection :: getDofIndicesForCurrentSide() {
      d_dofIndices.resize(d_currNodes.size());
      std::vector<AMP::Mesh::MeshElementID> globalIDs(d_currNodes.size()); 
      for(unsigned int j = 0; j < d_currNodes.size(); j++) {
        globalIDs[j] = d_currNodes[j].globalID();
      } // end of j
      d_dofManager->getDOFs(globalIDs, d_dofIndices);
    }

    AMP::LinearAlgebra::Vector::shared_ptr PressureBoundaryVectorCorrection :: mySubsetVector(
        AMP::LinearAlgebra::Vector::shared_ptr vec, 
        AMP::LinearAlgebra::Variable::shared_ptr var) {
      if(d_Mesh.get() != NULL) {
        AMP::LinearAlgebra::VS_Mesh meshSelector(d_Mesh);
        AMP::LinearAlgebra::Vector::shared_ptr meshSubsetVec = vec->select(meshSelector, (vec->getVariable())->getName());
        AMP::LinearAlgebra::Vector::shared_ptr varSubsetVec = meshSubsetVec->subsetVectorForVariable(var);
        return varSubsetVec;
      } else {
        return vec->subsetVectorForVariable(var);
      }
    }

    AMP::LinearAlgebra::Vector::const_shared_ptr PressureBoundaryVectorCorrection :: mySubsetVector(
        AMP::LinearAlgebra::Vector::const_shared_ptr vec, 
        AMP::LinearAlgebra::Variable::shared_ptr var) {
      if(d_Mesh.get() != NULL) {
        AMP::LinearAlgebra::VS_Mesh meshSelector(d_Mesh);
        AMP::LinearAlgebra::Vector::const_shared_ptr meshSubsetVec = vec->constSelect(meshSelector, (vec->getVariable())->getName());
        AMP::LinearAlgebra::Vector::const_shared_ptr varSubsetVec = meshSubsetVec->constSubsetVectorForVariable(var);
        return varSubsetVec;
      } else {
        return vec->constSubsetVectorForVariable(var);
      }
    }

  }
}


