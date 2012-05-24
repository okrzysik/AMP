
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
              std::string feFamilyName = (params->d_db)->getStringWithDefault("FE_FAMILY", "LAGRANGE");
              std::string qruleTypeName = (params->d_db)->getStringWithDefault("QRULE_TYPE", "QGAUSS");
              d_qruleOrderName = (params->d_db)->getStringWithDefault("QRULE_ORDER", "DEFAULT");

              d_feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>(feTypeOrderName);
              d_feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>(feFamilyName);
              d_qruleType = Utility::string_to_enum<libMeshEnums::QuadratureType>(qruleTypeName);

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
      d_numBndIds = (myparams->d_db)->getInteger("number_of_ids");

      d_isConstantFlux = (myparams->d_db)->getBoolWithDefault("constant_flux", true);

      d_boundaryIds.resize(d_numBndIds);
      d_dofIds.resize(d_numBndIds);
      d_neumannValues.resize(d_numBndIds);
      d_IsCoupledBoundary.resize(d_numBndIds);
      d_numDofIds.resize(d_numBndIds);

      char key[100];
      for(int j = 0; j < d_numBndIds; j++) {
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

      // Create the libmesh elements
      AMP::Mesh::MeshIterator iterator;
      for(unsigned int j = 0; j < d_boundaryIds.size() ; j++) {
        AMP::Mesh::MeshIterator iterator2 = d_Mesh->getBoundaryIDIterator( AMP::Mesh::Face, d_boundaryIds[j], 0 );
        iterator = AMP::Mesh::Mesh::getIterator( AMP::Mesh::Union, iterator, iterator2 );
      }
      libmeshElements.reinit( iterator );
    }

    void
      NeumannVectorCorrection :: addRHScorrection(AMP::LinearAlgebra::Vector::shared_ptr rhsCorrection)
      {

        AMP::LinearAlgebra::Vector::shared_ptr myRhs = subsetInputVector( rhsCorrection );
        std::vector<double> gamma(1,1.0);

        if(!d_isConstantFlux)
        {
          d_variableFlux->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
        }

        if((d_params->d_db)->keyExists("gamma"))
        {
          gamma[0] = (d_params->d_db)->getDouble("gamma");
        }

        AMP::LinearAlgebra::Vector::shared_ptr rInternal = myRhs->cloneVector();
        AMP::Discretization::DOFManager::shared_ptr dofManager = rInternal->getDOFManager();
        rInternal->zero();

        unsigned int numBndIds = d_boundaryIds.size();
        std::vector<size_t> dofs;
        std::vector<std::vector<size_t> > dofIndices;

        for(unsigned int j = 0; j < numBndIds ; j++)
        {
          if(!d_IsCoupledBoundary[j])
          {
            unsigned int numDofIds = d_dofIds[j].size();

            for(unsigned int k = 0; k < numDofIds; k++)
            {

              AMP::Mesh::MeshIterator bnd     = d_Mesh->getBoundaryIDIterator( AMP::Mesh::Face, d_boundaryIds[j], 0 );
              AMP::Mesh::MeshIterator end_bnd = bnd.end();

              int count =0;
              for( ; bnd != end_bnd; ++bnd)
              {
                count++;

                const unsigned int dimension = 2;

                boost::shared_ptr < ::FEType > d_feType ( new ::FEType(d_feTypeOrder, d_feFamily) );
                boost::shared_ptr < ::FEBase > d_fe ( (::FEBase::build(dimension, (*d_feType))).release() );

                if(d_qruleOrderName == "DEFAULT") {
                  d_qruleOrder = d_feType->default_quadrature_order();
                } else {
                  d_qruleOrder = Utility::string_to_enum<libMeshEnums::Order>(d_qruleOrderName);
                }
                boost::shared_ptr < ::QBase > d_qrule( (::QBase::build(d_qruleType, dimension, d_qruleOrder)).release() );


                d_currNodes = bnd->getElements(AMP::Mesh::Vertex);
                unsigned int numNodesInCurrElem = d_currNodes.size();

                dofIndices.resize(numNodesInCurrElem);
                for(unsigned int i = 0; i < numNodesInCurrElem ; i++) {
                   dofManager->getDOFs(d_currNodes[i].globalID(), dofIndices[i]);
                }

                // Get the libmesh element
                d_currElemPtr = libmeshElements.getElement( bnd->globalID() );

                d_fe->attach_quadrature_rule( d_qrule.get() );

                d_phi = &(d_fe->get_phi());
                d_JxW = &(d_fe->get_JxW());

                d_fe->reinit ( d_currElemPtr );

                const std::vector<std::vector<Real> > &phi = *d_phi;
                const std::vector<Real> &djxw = *d_JxW;

                dofs.resize(dofIndices.size());
                for(unsigned int i = 0; i < dofIndices.size() ; i++)
                  dofs[i] = dofIndices[i][k];

                std::vector<double> flux( dofIndices.size(), 0.0);

                for(unsigned int i = 0; i < dofIndices.size() ; i++)    // Loop over nodes
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
                      temp[0].push_back(d_variableFlux->getValueByGlobalID(dofs[i]));
                    }

                    if(d_robinPhysicsModel)
                    {
                      d_robinPhysicsModel->getConductance(gamma, gamma, temp); 
                    }

                    flux[i] +=  (gamma[0])*djxw[qp]*phi[i][qp]*temp[0][0];

                  }//end for qp

                }//end for i

                rInternal->addValuesByGlobalID((int)dofs.size() , (size_t *)&(dofs[0]), &(flux[0]));

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
        tmp_db->putInteger("number_of_ids",d_numBndIds);
        tmp_db->putBool("constant_flux", d_isConstantFlux);

        char key[100];
        for(int j = 0; j < d_numBndIds; j++) {
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

