
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

#include <string>

namespace AMP {
namespace Operator {

  void PressureBoundaryVectorCorrection :: reset(const boost::shared_ptr<OperatorParameters>& params)
  {
    boost::shared_ptr<PressureBoundaryVectorCorrectionParameters> myparams = 
      boost::dynamic_pointer_cast<PressureBoundaryVectorCorrectionParameters>(params);

    AMP_INSIST( ((myparams.get()) != NULL), "NULL parameters" );
    AMP_INSIST( (((myparams->d_db).get()) != NULL), "NULL database" );

    AMP_INSIST( (myparams->d_db)->keyExists("number_of_ids"), "Key ''number_of_ids'' is missing!" );
    d_numIds = (myparams->d_db)->getInteger("number_of_ids");

    d_isConstantPressure = (myparams->d_db)->getBoolWithDefault("constant_pressure", true);

    d_boundaryIds.resize(d_numIds);
    d_pressureValues.resize(d_numIds);

    char key[100];
    for(int j = 0; j < d_numIds; j++) {
      sprintf(key, "id_%d", j);
      AMP_INSIST( (myparams->d_db)->keyExists(key), "Key is missing!" );
      d_boundaryIds[j] = (myparams->d_db)->getInteger(key);

      if(d_isConstantPressure){
        sprintf(key, "value_%d", j);
        AMP_INSIST( (myparams->d_db)->keyExists(key), "Key is missing!" );
        d_pressureValues[j] = (myparams->d_db)->getDouble(key);
      }else{
        AMP_INSIST((d_isConstantPressure ==  true), "Variable pressure has not been implemented yet.");
        d_variablePressure = myparams->d_variablePressure;      
      }
    }//end for j
  }

  std::vector<std::vector<std::vector<Point > > > PressureBoundaryVectorCorrection :: getNormals()
  {
    if(!d_isConstantPressure){
      d_variablePressure->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    }

    if(d_iDebugPrintInfoLevel>8)
    {
      AMP::pout << "Entered the getNormals function." << std::endl;
    }

    unsigned int numIds = d_boundaryIds.size();

    std::vector<std::vector<std::vector<Point> > > allBoundaryIdsNormals;

    for(unsigned int j = 0; j < numIds; j++) {
      if(d_iDebugPrintInfoLevel>8)
      {
        AMP::pout << "Entered the numIds loop. d_boundaryIds[" << j << "] = " << d_boundaryIds[j] << std::endl;
      }
      std::vector<std::vector<Point> > thisBoundaryIdNormals;

      //This could be changed to BoundaryOwnNodeIterator
      AMP::Mesh::MeshManager::Adapter::BoundarySideIterator bnd = d_MeshAdapter->beginSideBoundary( d_boundaryIds[j] );
      AMP::Mesh::MeshManager::Adapter::BoundarySideIterator end_bnd = d_MeshAdapter->endSideBoundary( d_boundaryIds[j] );

      for( ; bnd != end_bnd; ++bnd) {
  
        if(d_iDebugPrintInfoLevel>8)
        {
          AMP::pout << "Entered the bnd loop." << std::endl;
        }

        AMP::Mesh::MeshManager::Adapter::Element cur_elem = d_MeshAdapter->getElementFromSide(*bnd);
        AMP::Mesh::MeshManager::Adapter::Element cur_side = *bnd;

        d_normal = &(d_fe_3d->get_normals());

        d_elem = &cur_side.getElem();
        e_elem = &cur_elem.getElem();

        if(d_iDebugPrintInfoLevel>8)
        {
          AMP::pout << "Number of sides in d_elem = " << d_elem->n_sides() << " in e_elem =  " << e_elem->n_sides() << std::endl;
        }

        unsigned int k = 0;
        for(unsigned int s = 0; s < e_elem->n_sides(); s++) {
          AutoPtr< ::Elem> side (e_elem->build_side(s));
          if(*(side.get()) == *d_elem) {
            if(d_iDebugPrintInfoLevel>8)
            {
              AMP::pout << "The side that matched is = " << s << std::endl;
            }
            k = s;
          }
        }
        if(d_iDebugPrintInfoLevel>8)
        {
          AMP::pout << "The side that matched is = " << k << std::endl;
        }

        d_fe->reinit ( d_elem );
        d_fe_3d->reinit( e_elem, k );

        const std::vector<Point> & normal = *d_normal;
        thisBoundaryIdNormals.push_back(normal);

      }//end for bnd

      allBoundaryIdsNormals.push_back(thisBoundaryIdNormals);

    }//end for j

    return allBoundaryIdsNormals;
  }
  void PressureBoundaryVectorCorrection :: computeRHScorrection()
  {
    if(!d_isConstantPressure){
      d_variablePressure->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    }

    if(d_iDebugPrintInfoLevel>8)
    {
      AMP::pout << "Entered the computeRHScorrection function." << std::endl;
    }

    AMP::Mesh::DOFMap::shared_ptr dof_map = d_MeshAdapter->getDOFMap(d_variable);

    AMP::LinearAlgebra::Vector::shared_ptr rInternal = d_MeshAdapter->createVector(d_variable);
    rInternal->zero();

    unsigned int numIds = d_boundaryIds.size();

    for(unsigned int j = 0; j < numIds; j++) {
      if(d_iDebugPrintInfoLevel>8)
      {
        AMP::pout << "Entered the numIds loop. d_boundaryIds[" << j << "] = " << d_boundaryIds[j] << std::endl;
      }
      //This could be changed to BoundaryOwnNodeIterator
      AMP::Mesh::MeshManager::Adapter::BoundarySideIterator bnd = d_MeshAdapter->beginSideBoundary( d_boundaryIds[j] );
      AMP::Mesh::MeshManager::Adapter::BoundarySideIterator end_bnd = d_MeshAdapter->endSideBoundary( d_boundaryIds[j] );

      for( ; bnd != end_bnd; ++bnd) {
  
        if(d_iDebugPrintInfoLevel>8)
        {
          AMP::pout << "Entered the bnd loop." << std::endl;
        }

        AMP::Mesh::MeshManager::Adapter::Element cur_elem = d_MeshAdapter->getElementFromSide(*bnd);
        AMP::Mesh::MeshManager::Adapter::Element cur_side = *bnd;

        d_phi = &(d_fe->get_phi());
        d_JxW = &(d_fe->get_JxW());
        d_normal = &(d_fe_3d->get_normals());

        d_elem = &cur_side.getElem();
        e_elem = &cur_elem.getElem();

        if(d_iDebugPrintInfoLevel>8)
        {
          AMP::pout << "Number of sides in d_elem = " << d_elem->n_sides() << " in e_elem =  " << e_elem->n_sides() << std::endl;
        }

        unsigned int k = 0;
        for(unsigned int s = 0; s < e_elem->n_sides(); s++) {
          AutoPtr< ::Elem> side (e_elem->build_side(s));
          if(*(side.get()) == *d_elem) {
            if(d_iDebugPrintInfoLevel>8)
            {
              AMP::pout << "The side that matched is = " << s << std::endl;
            }
            k = s;
          }
        }
        if(d_iDebugPrintInfoLevel>8)
        {
          AMP::pout << "The side that matched is = " << k << std::endl;
        }

        d_fe->reinit ( d_elem );
        d_fe_3d->reinit( e_elem, k );

        const std::vector<std::vector<Real> > & phi = *d_phi;
        const std::vector<Real> & djxw = *d_JxW;
        const std::vector<Point> & normal = *d_normal;

        for(int dofId = 0; dofId < 3; dofId++) {
          std::vector<unsigned int> bndGlobalIds;
          dof_map->getDOFs(cur_side, bndGlobalIds, dofId);

          if(d_iDebugPrintInfoLevel>8)
          {
            AMP::pout << "Entered the dofId loop." << std::endl;
            AMP::pout << "bndGlobalIds.size() = " << bndGlobalIds.size() << " d_qrule->n_points() =  " << d_qrule->n_points() << " d_normal->size() = " << d_normal->size() << std::endl;
            AMP::pout << "d_JxW->size() = " << d_JxW->size() << " d_phi->size() = " << d_phi->size() << std::endl;
          }

          std::vector<double> pressure(bndGlobalIds.size(), 0.0);

          for(unsigned int i = 0; i < bndGlobalIds.size(); i++) {
            for(unsigned int qp = 0; qp < d_qrule->n_points(); qp++) 
            {
              //std::vector<std::vector<double> > temp(1) ;
              if(d_iDebugPrintInfoLevel>8)
              {
                AMP::pout << "The weights are as follows : "  << djxw[qp] << " shape function " << phi[i][qp] << " pressure " << d_pressureValues[j] << std::endl;
                AMP::pout << "Direction of the normals are as follows : " << normal[qp](0) << " , " << normal[qp](1) << " , " << normal[qp](2) << std::endl;
                AMP::pout << " bndGlobalIds[" << i << "] = " << bndGlobalIds[i] << std::endl;
              }

              if(d_isConstantPressure){
                if(dofId == 0) pressure[i] +=  (djxw[qp]*phi[i][qp]*d_pressureValues[j]*normal[qp](0));
                if(dofId == 1) pressure[i] +=  (djxw[qp]*phi[i][qp]*d_pressureValues[j]*normal[qp](1));
                if(dofId == 2) pressure[i] +=  (djxw[qp]*phi[i][qp]*d_pressureValues[j]*normal[qp](2));
                if(d_iDebugPrintInfoLevel>8)
                {
                  AMP::pout << "dofID = " << dofId << " pressure[" << i << "] = " << pressure[i] << std::endl;
                }
              }else{
                AMP_INSIST((d_isConstantPressure == true), "Variable pressure has not been implemented yet.");
              }
            }

          }//end for i
          rInternal->addValuesByGlobalID((int)bndGlobalIds.size() , (int *)&(bndGlobalIds[0]), &(pressure[0]));
        }//end for dofId

      }//end for bnd

    }//end for j
    ////////////////////////////////////////////////////////////////

    rInternal->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_ADD );

    d_rhsCorrectionAdd = rInternal->cloneVector();
    d_rhsCorrectionAdd->copyVector(rInternal);
  }

  void PressureBoundaryVectorCorrection :: apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr  &r, const double a , const double b )
  {
    AMP::LinearAlgebra::Vector::shared_ptr rInternal = r->subsetVectorForVariable(d_variable);

    if(d_iDebugPrintInfoLevel>3)
    {
      AMP::pout << "L2 Norm of rInternal entering PressureBoundaryVectorCorrection::apply is : " << rInternal->L2Norm() << std::endl;
      AMP::pout << "The scaling value is : " << a << std::endl;
    }

    computeRHScorrection();
    d_rhsCorrectionAdd->scale(a);
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
        if(d_isConstantPressure){
          sprintf(key, "value_%d", j);
          tmp_db->putInteger(key,d_pressureValues[j]);
        }else{
          //d_variableFlux ??
        }
      }

      tmp_db->putBool("skip_params", true);

      boost::shared_ptr<PressureBoundaryVectorCorrectionParameters> outParams(new PressureBoundaryVectorCorrectionParameters(tmp_db));

      return outParams;
    }

}
}

