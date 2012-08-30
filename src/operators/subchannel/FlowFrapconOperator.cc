
#include "FlowFrapconOperator.h"
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

void FlowFrapconOperator :: reset(const boost::shared_ptr<OperatorParameters>& params)
{
      boost::shared_ptr<FlowFrapconOperatorParameters> myparams = 
        boost::dynamic_pointer_cast<FlowFrapconOperatorParameters>(params);

      AMP_INSIST( ((myparams.get()) != NULL), "NULL parameters" );
      AMP_INSIST( (((myparams->d_db).get()) != NULL), "NULL database" );

      AMP_INSIST( (myparams->d_db)->keyExists("numpoints"), "Key ''numpoints'' is missing!" );
      d_numpoints = (myparams->d_db)->getInteger("numpoints");

      AMP_INSIST( (myparams->d_db)->keyExists("Channel_Diameter"), "Missing key: Channel_Dia" );
      d_De = (myparams->d_db)->getDouble("Channel_Diameter");

      AMP_INSIST( (myparams->d_db)->keyExists("Heat_Capacity"), "Missing key: Heat_Capacity" );
      Cp = (myparams->d_db)->getDouble("Heat_Capacity");

      AMP_INSIST( (myparams->d_db)->keyExists("Mass_Flux"), "Missing key: Mass_Flux" );
      d_G  = (myparams->d_db)->getDouble("Mass_Flux");

      AMP_INSIST( (myparams->d_db)->keyExists("Temp_Inlet"), "Missing key: Temp_In" );
      d_Tin = (myparams->d_db)->getDoubleWithDefault("Temp_Inlet", 300.);

      AMP_INSIST( (myparams->d_db)->keyExists("Conductivity"), "Missing key: Kconductivity" );
      d_K   = (myparams->d_db)->getDouble("Conductivity");

      AMP_INSIST( (myparams->d_db)->keyExists("Reynolds"), "Missing key: Reynolds" );
      d_Re  = (myparams->d_db)->getDouble("Reynolds");

      AMP_INSIST( (myparams->d_db)->keyExists("Prandtl"), "Missing key: Prandtl" );
      d_Pr  = (myparams->d_db)->getDouble("Prandtl");

}


// This is an in-place apply
void FlowFrapconOperator :: apply(AMP::LinearAlgebra::Vector::const_shared_ptr f, AMP::LinearAlgebra::Vector::const_shared_ptr u,
    AMP::LinearAlgebra::Vector::shared_ptr r, const double a, const double b)
{

      // AMP::Mesh::DOFMap::shared_ptr dof_map = d_MeshAdapter->getDOFMap(d_inpVariable);

      AMP_INSIST( ((r.get()) != NULL), "NULL Residual Vector" );
      AMP_INSIST( ((u.get()) != NULL), "NULL Solution Vector" );

      if( !zPoints.empty() ){d_numpoints = zPoints.size();}

      std::vector<double> box = d_Mesh->getBoundingBox();
      const double min_z = box[4];
      const double max_z = box[5];
      const double del_z = (max_z-min_z)/d_numpoints ; 

      // std::cout << "Extreme Min Point in z = " << min_z << std::endl;
      // std::cout << "Extreme Max Point in z = " << max_z << std::endl;

      // Subset the vectors
      AMP::LinearAlgebra::Vector::const_shared_ptr flowInputVec = subsetInputVector( u );
      AMP::LinearAlgebra::Vector::shared_ptr outputVec = subsetOutputVector( r );

      // AMP::LinearAlgebra::Variable::shared_ptr localVar ( new AMP::LinearAlgebra::Variable(d_cladVec->getVariable()->getName() ) ); 
      // d_localCladVec = AMP::LinearAlgebra::SimpleVector::create( d_numpoints, localVar ); 
      // 
      // boost::shared_ptr<AMP::InputDatabase> map3to1_db (new AMP::InputDatabase("Dummy"));
      // map3to1_db->putInteger("BoundaryId",4);
      // map3to1_db->putString("InputVariable",localVar->getName());
      // map3to1_db->putString("OutputVariable",localVar->getName());
      // boost::shared_ptr<AMP::Operator::MapOperatorParameters> map3to1Params (new AMP::Operator::MapOperatorParameters( map3to1_db ));
      // map3to1Params->d_MeshAdapter = d_MeshAdapter;
      // boost::shared_ptr<AMP::Operator::Map3Dto1D> mapCladto1 (new AMP::Operator::Map3Dto1D( map3to1Params ));
      // 
      // mapCladto1->setZLocations( d_Map1to3->getZLocations());
      // 
      // mapCladto1->setVector(d_localCladVec); 
      // mapCladto1->apply(nullVec, d_cladVec, nullVec);

      // should never use resize as it is assumed that the vector is created using the right size !! 
      //outputVec->castTo<SimpleVector>().resize (d_numpoints);

      zPoints.resize(d_numpoints);

      // set the inlet flow temperature value
      double T1 = flowInputVec->getValueByLocalID(0);
      outputVec->setValueByLocalID(0, T1-d_Tin);

      zPoints[0] = min_z;
      for( int j=1; j<d_numpoints; j++) {
        zPoints[j] = zPoints[j-1] + del_z;
      } 

      // Iterate through the flow boundary
      for( int i=1; i<d_numpoints; i++) {

        double cur_node, next_node;

        cur_node  = zPoints[i-1];
        next_node = zPoints[i];

        double Heff, he_z, T_b_i, T_b_im1, T_c_i;
        double R_b =0;

        T_c_i = d_cladVec->getValueByLocalID(i);
        T_b_i = flowInputVec->getValueByLocalID(i);
        T_b_im1 = flowInputVec->getValueByLocalID(i-1);

        Heff = (0.023*d_K/d_De)*pow(d_Re,0.8)*pow(d_Pr,0.4);
        //       Cp   = getHeatCapacity(T_b_i);
        he_z = next_node - cur_node;

        R_b = T_b_i - T_b_im1 - ((4*Heff*( T_c_i - T_b_i))/(Cp*d_G*d_De))*he_z;

        outputVec->setValueByLocalID(i, R_b);

      }//end for i

      if(f.get() == NULL) {
        outputVec->scale(a);
      } else {
        AMP::LinearAlgebra::Vector::const_shared_ptr fInternal = subsetInputVector( f );
        if(fInternal.get() == NULL) {
          outputVec->scale(a);
        } else {
          outputVec->axpby(b, a, fInternal);
        }
      }

}


boost::shared_ptr<OperatorParameters> FlowFrapconOperator :: 
    getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& u) 
{
        boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));

        tmp_db->putString("name","FlowFrapconOperator");
        tmp_db->putInteger("numpoints",d_numpoints);
        tmp_db->putDouble("Channel_Diameter",d_De);
        tmp_db->putDouble("Mass_Flux",d_G);
        tmp_db->putDouble("Heat_Capacity",Cp);
        tmp_db->putDouble("Temp_Inlet",d_Tin);
        tmp_db->putDouble("Conductivity",d_K);
        tmp_db->putDouble("Reynolds",d_Re);
        tmp_db->putDouble("Prandtl",d_Pr);

        boost::shared_ptr<FlowFrapconJacobianParameters> outParams(new FlowFrapconJacobianParameters(tmp_db));
        outParams->d_frozenSolution = subsetInputVector(u); 
        return outParams;
}


AMP::LinearAlgebra::Vector::shared_ptr  FlowFrapconOperator::subsetOutputVector(AMP::LinearAlgebra::Vector::shared_ptr vec)
{
    AMP::LinearAlgebra::Variable::shared_ptr var = getInputVariable();
    // Subset the vectors, they are simple vectors and we need to subset for the current comm instead of the mesh
    if(d_Mesh.get() != NULL) {
        AMP::LinearAlgebra::VS_Comm commSelector( d_Mesh->getComm() );
        AMP::LinearAlgebra::Vector::shared_ptr commVec = vec->select(commSelector, var->getName());
        return commVec->subsetVectorForVariable(var);
    } else {
        return vec->subsetVectorForVariable(var);
    }
}


AMP::LinearAlgebra::Vector::shared_ptr  FlowFrapconOperator::subsetInputVector(AMP::LinearAlgebra::Vector::shared_ptr vec)
{
    AMP::LinearAlgebra::Variable::shared_ptr var = getInputVariable();
    // Subset the vectors, they are simple vectors and we need to subset for the current comm instead of the mesh
    if(d_Mesh.get() != NULL) {
        AMP::LinearAlgebra::VS_Comm commSelector( d_Mesh->getComm() );
        AMP::LinearAlgebra::Vector::shared_ptr commVec = vec->select(commSelector, var->getName());
        return commVec->subsetVectorForVariable(var);
    } else {
        return vec->subsetVectorForVariable(var);
    }
}


AMP::LinearAlgebra::Vector::const_shared_ptr  FlowFrapconOperator::subsetOutputVector(AMP::LinearAlgebra::Vector::const_shared_ptr vec)
{
    AMP::LinearAlgebra::Variable::shared_ptr var = getInputVariable();
    // Subset the vectors, they are simple vectors and we need to subset for the current comm instead of the mesh
    if(d_Mesh.get() != NULL) {
        AMP::LinearAlgebra::VS_Comm commSelector( d_Mesh->getComm() );
        AMP::LinearAlgebra::Vector::const_shared_ptr commVec = vec->constSelect(commSelector, var->getName());
        return commVec->constSubsetVectorForVariable(var);
    } else {
        return vec->constSubsetVectorForVariable(var);
    }
}


AMP::LinearAlgebra::Vector::const_shared_ptr  FlowFrapconOperator::subsetInputVector(AMP::LinearAlgebra::Vector::const_shared_ptr vec)
{
    AMP::LinearAlgebra::Variable::shared_ptr var = getInputVariable();
    // Subset the vectors, they are simple vectors and we need to subset for the current comm instead of the mesh
    if(d_Mesh.get() != NULL) {
        AMP::LinearAlgebra::VS_Comm commSelector( d_Mesh->getComm() );
        AMP::LinearAlgebra::Vector::const_shared_ptr commVec = vec->constSelect(commSelector, var->getName());
        return commVec->constSubsetVectorForVariable(var);
    } else {
        return vec->constSubsetVectorForVariable(var);
    }
}


}
}



