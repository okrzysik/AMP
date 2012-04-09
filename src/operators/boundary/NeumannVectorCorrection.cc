
#include "NeumannVectorCorrection.h"
#include "NeumannVectorCorrectionParameters.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"

/* Libmesh files */

#include "enum_order.h"
#include "enum_fe_family.h"
#include "enum_quadrature_type.h"
#include "auto_ptr.h"
#include "string_to_enum.h"

#include "face_quad4.h"
#include "node.h"

#include <string>

namespace AMP {
namespace Operator {


    // Constructor
    NeumannVectorCorrection::NeumannVectorCorrection(const boost::shared_ptr<NeumannVectorCorrectionParameters> & params)
        : BoundaryOperator (params)
      {
	      d_params = params;

              std::string feTypeOrderName = (params->d_db)->getStringWithDefault("FE_ORDER", "FIRST");

              libMeshEnums::Order feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>(feTypeOrderName);

              std::string feFamilyName = (params->d_db)->getStringWithDefault("FE_FAMILY", "LAGRANGE");

              libMeshEnums::FEFamily feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>(feFamilyName);

              std::string qruleTypeName = (params->d_db)->getStringWithDefault("QRULE_TYPE", "QGAUSS");

              libMeshEnums::QuadratureType qruleType = Utility::string_to_enum<libMeshEnums::QuadratureType>(qruleTypeName);

              const unsigned int dimension = 2;

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

              d_JxW = &(d_fe->get_JxW());

              d_variable = params->d_variable;

              reset(params);
      }


    void NeumannVectorCorrection :: reset(const boost::shared_ptr<OperatorParameters>& params)
    {
      boost::shared_ptr<NeumannVectorCorrectionParameters> myparams = 
        boost::dynamic_pointer_cast<NeumannVectorCorrectionParameters>(params);

      AMP_INSIST( ((myparams.get()) != NULL), "NULL parameters" );
      AMP_INSIST( (((myparams->d_db).get()) != NULL), "NULL database" );

      AMP_INSIST( (myparams->d_db)->keyExists("number_of_ids"), "Key ''number_of_ids'' is missing!" );
      d_numIds = (myparams->d_db)->getInteger("number_of_ids");

      d_isConstantFlux = (myparams->d_db)->getBoolWithDefault("constant_flux", true);

      d_boundaryIds.resize(d_numIds);
      d_dofIds.resize(d_numIds);
      d_neumannValues.resize(d_numIds);
      d_IsCoupledBoundary.resize(d_numIds);
      d_numDofIds.resize(d_numIds);

      char key[100];
      for(int j = 0; j < d_numIds; j++) {
        sprintf(key, "id_%d", j);
        AMP_INSIST( (myparams->d_db)->keyExists(key), "Key is missing!" );
        d_boundaryIds[j] = (myparams->d_db)->getInteger(key);

        sprintf(key, "number_of_dofs_%d", j);
        AMP_INSIST( (myparams->d_db)->keyExists(key), "Key is missing!" );
        d_numDofIds[j] = (myparams->d_db)->getInteger(key);

        sprintf(key, "IsCoupledBoundary_%d", j);
        d_IsCoupledBoundary[j] = (params->d_db)->getBoolWithDefault(key, false);

        d_dofIds[j].resize(d_numDofIds[j]);
        d_neumannValues[j].resize(d_numDofIds[j]);
        for(int i = 0; i < d_numDofIds[j]; i++) {
          sprintf(key, "dof_%d_%d", j, i);
          AMP_INSIST( (myparams->d_db)->keyExists(key), "Key is missing!" );
          d_dofIds[j][i] = (myparams->d_db)->getInteger(key);

          if(d_isConstantFlux){
            sprintf(key, "value_%d_%d", j, i);
            AMP_INSIST( (myparams->d_db)->keyExists(key), "Key is missing!" );
            d_neumannValues[j][i] = (myparams->d_db)->getDouble(key);
          }else{
            d_variableFlux = myparams->d_variableFlux;      
          }
        }//end for i
      }//end for j

      if(myparams->d_robinPhysicsModel) {
        d_robinPhysicsModel = myparams->d_robinPhysicsModel;
      }
    }

    void
      NeumannVectorCorrection :: addRHScorrection(AMP::LinearAlgebra::Vector::shared_ptr rhsCorrection)
      {

        AMP::LinearAlgebra::Vector::shared_ptr myRhs = rhsCorrection->subsetVectorForVariable(d_variable);
        d_gamma.resize(1);
        d_gamma[0] = 1.0;

        if(!d_isConstantFlux)
        {
          d_variableFlux->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
        }

        if((d_params->d_db)->keyExists("gamma"))
        {
          d_gamma[0] = (d_params->d_db)->getDouble("gamma");
        }

        AMP::LinearAlgebra::Vector::shared_ptr rInternal = myRhs->cloneVector();
        AMP::Discretization::DOFManager::shared_ptr dofManager = rInternal->getDOFManager();
        rInternal->zero();

        unsigned int numIds = d_boundaryIds.size();
        std::vector<size_t> dofIndices, dofs;
        for(unsigned int j = 0; j < numIds; j++)
        {
          if(!d_IsCoupledBoundary[j])
          {
            unsigned int numDofIds = d_dofIds[j].size();

            for(unsigned int k = 0; k < numDofIds; k++)
            {

              unsigned int dofId;
              dofId = d_dofIds[j][k];

              AMP::Mesh::MeshIterator bnd     = d_Mesh->getIDsetIterator( AMP::Mesh::Face, d_boundaryIds[j], 0 );
              AMP::Mesh::MeshIterator end_bnd = bnd.end();

              int count =0;
              for( ; bnd != end_bnd; ++bnd)
              {
                count++;

                d_currNodes = bnd->getElements(AMP::Mesh::Vertex);
                unsigned int numNodesInCurrElem = d_currNodes.size();

                dofIndices.resize(numNodesInCurrElem);
                for(unsigned int i = 0; i < numNodesInCurrElem ; i++) {
                    dofManager->getDOFs(d_currNodes[i].globalID(), dofs);
                    AMP_ASSERT(dofs.size()==1);
                    dofIndices[i] = dofs[0];
                }

                createCurrentLibMeshElement();

                getDofIndicesForCurrentElement();

                d_phi = &(d_fe->get_phi());
                d_JxW = &(d_fe->get_JxW());

                d_fe->reinit ( d_currElemPtr );

                std::vector<double> flux( numNodesInCurrElem , 0.0);
                const std::vector<std::vector<Real> > &phi = *d_phi;
                const std::vector<Real> &djxw = *d_JxW;

                for(unsigned int i = 0; i < numNodesInCurrElem ; i++)
                {
                  for(unsigned int qp = 0; qp < d_qrule->n_points(); qp++) 
                  {
                    std::vector<std::vector<double> > temp(1) ;

                    // there must be a better way to write this!!
                    if(d_isConstantFlux)
                    {
                      temp[0].push_back(d_neumannValues[j][k]);            
                    }
                    else
                    {
                      temp[0].push_back(d_variableFlux->getValueByGlobalID(dofIndices[i]));
                    }

                    if(d_robinPhysicsModel)
                    {
                      d_robinPhysicsModel->getConductance(d_gamma, d_gamma, temp); 
                    }

                    flux[i] +=  (d_gamma[0])*djxw[qp]*phi[i][qp]*temp[0][0];

                  }

                }//end for i

                rInternal->addValuesByGlobalID((int)dofIndices.size() , (size_t *)&(dofIndices[0]), &(flux[0]));

                destroyCurrentLibMeshElement();
              }//end for bnd

            }//end for k
          }//coupled
        }//end for j

        rInternal->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_ADD );
        myRhs->add(myRhs, rInternal);

      }

    void NeumannVectorCorrection :: apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr  &r, const double a , const double b )
    {
      (void) f; (void) u; (void) r; (void) a; (void) b; 
      //Do Nothing
    }

    boost::shared_ptr<OperatorParameters> NeumannVectorCorrection :: 
      getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& ) {
        boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));

        tmp_db->putString("FE_ORDER","FIRST");
        tmp_db->putString("FE_FAMILY","LAGRANGE");
        tmp_db->putString("QRULE_TYPE","QGAUSS");
        tmp_db->putInteger("DIMENSION",2);
        tmp_db->putString("QRULE_ORDER","DEFAULT");
        tmp_db->putInteger("number_of_ids",d_numIds);
        tmp_db->putBool("constant_flux", d_isConstantFlux);

        char key[100];
        for(int j = 0; j < d_numIds; j++) {
          sprintf(key, "id_%d", j);
          tmp_db->putInteger(key,d_boundaryIds[j]);
          sprintf(key, "number_of_dofs_%d", j);
          tmp_db->putInteger(key,d_numDofIds[j]);
          sprintf(key, "IsCoupledBoundary_%d", j);
          tmp_db->putBool(key, d_IsCoupledBoundary[j]);

          for(int i = 0; i < d_numDofIds[j]; i++) {
            sprintf(key, "dof_%d_%d", j, i);
            tmp_db->putInteger(key,d_dofIds[j][i]);
            if(d_isConstantFlux){
              sprintf(key, "value_%d_%d", j, i);
              tmp_db->putInteger(key,d_neumannValues[j][i]);
            }else{
              //d_variableFlux ??
            }
          }
        }

        tmp_db->putBool("skip_params", true);

        boost::shared_ptr<NeumannVectorCorrectionParameters> outParams(new NeumannVectorCorrectionParameters(tmp_db));

        return outParams;
      }

    void NeumannVectorCorrection :: createCurrentLibMeshElement() {
      d_currElemPtr = new ::Quad4;
      for(size_t j = 0; j < d_currNodes.size(); j++) {
        std::vector<double> pt = d_currNodes[j].coord();
        d_currElemPtr->set_node(j) = new ::Node(pt[0], pt[1], pt[2], j);
      }//end for j
    }

    void NeumannVectorCorrection :: destroyCurrentLibMeshElement() {
      for(size_t j = 0; j < d_currElemPtr->n_nodes(); j++) {
        delete (d_currElemPtr->get_node(j));
        d_currElemPtr->set_node(j) = NULL;
      }//end for j
      delete d_currElemPtr;
      d_currElemPtr = NULL;
    }

    void NeumannVectorCorrection :: getDofIndicesForCurrentElement() {
      std::vector<AMP::Mesh::MeshElementID> globalIDs(d_currNodes.size()); 
      for(unsigned int j = 0; j < d_currNodes.size(); j++) {
        globalIDs[j] = d_currNodes[j].globalID();
      } // end of j
    }


     void NeumannVectorCorrection :: setFrozenVector ( AMP::LinearAlgebra::Vector::shared_ptr f ) {
         if ( d_Frozen )
         {
         d_Frozen->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( f );
         }
         else
         {
         d_Frozen = AMP::LinearAlgebra::MultiVector::view ( f );
         }

         AMP::LinearAlgebra::VS_Mesh meshSelector("meshSelector", d_Mesh);
         d_Frozen = d_Frozen->select(meshSelector, d_Frozen->getVariable()->getName());

    }


    void NeumannVectorCorrection :: setVariableFlux(const AMP::LinearAlgebra::Vector::shared_ptr &flux) {
        if(d_Mesh.get() != NULL) {
            AMP::LinearAlgebra::VS_Mesh meshSelector( d_variable->getName(), d_Mesh);
            AMP::LinearAlgebra::Vector::shared_ptr meshSubsetVec = flux->select(meshSelector, d_variable->getName());
            d_variableFlux = meshSubsetVec->subsetVectorForVariable( d_variable );
        } else {
            d_variableFlux = flux->subsetVectorForVariable ( d_variable );
        }
    }


}
}

