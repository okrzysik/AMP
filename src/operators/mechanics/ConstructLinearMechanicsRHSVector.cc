
#include "ConstructLinearMechanicsRHSVector.h"

   //This file has not been converted!
  /*
void computeTemperatureRhsVector( AMP::Mesh::MeshManager::Adapter::shared_ptr mesh,
    boost::shared_ptr<AMP::Database> input_db,  
    AMP::LinearAlgebra::Variable::shared_ptr temperatureVar, 
    AMP::LinearAlgebra::Variable::shared_ptr displacementVar,
    const boost::shared_ptr<AMP::LinearAlgebra::Vector> &currTemperatureVec,
    const boost::shared_ptr<AMP::LinearAlgebra::Vector> &prevTemperatureVec,
    boost::shared_ptr<AMP::LinearAlgebra::Vector> &rhsVec) 
{
  currTemperatureVec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
  prevTemperatureVec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );

  AMP::LinearAlgebra::Vector::shared_ptr rInternal = rhsVec->subsetVectorForVariable(displacementVar);
  rInternal->zero();

  boost::shared_ptr<AMP::Database> elementRhsDatabase = input_db->getDatabase("RhsElements");
  boost::shared_ptr<AMP::Database> materialModelDatabase = input_db->getDatabase("RhsMaterialModel");

  boost::shared_ptr < ::FEType > d_feType;
  boost::shared_ptr < ::FEBase > d_fe;
  boost::shared_ptr < ::QBase > d_qrule;
  ::Elem *d_elem;

  std::string feTypeOrderName = elementRhsDatabase->getStringWithDefault("FE_ORDER", "FIRST");
  libMeshEnums::Order feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>(feTypeOrderName);

  std::string feFamilyName = elementRhsDatabase->getStringWithDefault("FE_FAMILY", "LAGRANGE");
  libMeshEnums::FEFamily feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>(feFamilyName);

  std::string qruleTypeName = elementRhsDatabase->getStringWithDefault("QRULE_TYPE", "QGAUSS");
  libMeshEnums::QuadratureType qruleType = Utility::string_to_enum<libMeshEnums::QuadratureType>(qruleTypeName);

  const unsigned int dimension = 3;

  d_feType.reset( new ::FEType(feTypeOrder, feFamily) );
  d_fe.reset( (::FEBase::build(dimension, (*d_feType))).release() );

  std::string qruleOrderName = elementRhsDatabase->getStringWithDefault("QRULE_ORDER", "DEFAULT");
  libMeshEnums::Order qruleOrder;
  if(qruleOrderName == "DEFAULT") {
    qruleOrder = d_feType->default_quadrature_order();
  } else {
    qruleOrder = Utility::string_to_enum<libMeshEnums::Order>(qruleOrderName);
  }

  d_qrule.reset( (::QBase::build(qruleType, dimension, qruleOrder)).release() );
  d_fe->attach_quadrature_rule( d_qrule.get() );

  int d_iDebugPrintInfoLevel = elementRhsDatabase->getIntegerWithDefault("print_info_level",0);

  const std::vector<Real> & JxW = (d_fe->get_JxW());
  const std::vector<std::vector<RealGradient> > & dphi = (d_fe->get_dphi());
  const std::vector<std::vector<Real> > & phi = (d_fe->get_phi());

  boost::shared_ptr<AMP::Materials::Material> d_material;
  double youngsModulus = 1.0e10, poissonsRatio = 0.33, thermalExpansionCoefficient = 2.0e-6;
  //double default_TEMPERATURE, default_BURNUP, default_OXYGEN_CONCENTRATION;
  double default_BURNUP, default_OXYGEN_CONCENTRATION;

  bool d_useMaterialsLibrary = materialModelDatabase->getBoolWithDefault("USE_MATERIALS_LIBRARY",false);
  if(d_useMaterialsLibrary == true) {
    AMP_INSIST( (materialModelDatabase->keyExists("Material")), "Key ''Material'' is missing!" );
    std::string matname = materialModelDatabase->getString("Material");
    d_material = AMP::voodoo::Factory<AMP::Materials::Material>::instance().create(matname);
  }

  if(d_useMaterialsLibrary == false) {
    AMP_INSIST( materialModelDatabase->keyExists("THERMAL_EXPANSION_COEFFICIENT"), "Missing key: THERMAL_EXPANSION_COEFFICIENT" );
    AMP_INSIST( materialModelDatabase->keyExists("Youngs_Modulus"), "Missing key: Youngs_Modulus" );
    AMP_INSIST( materialModelDatabase->keyExists("Poissons_Ratio"), "Missing key: Poissons_Ratio" );

    thermalExpansionCoefficient = materialModelDatabase->getDouble("THERMAL_EXPANSION_COEFFICIENT");
    youngsModulus = materialModelDatabase->getDouble("Youngs_Modulus");
    poissonsRatio = materialModelDatabase->getDouble("Poissons_Ratio");
  }

  //default_TEMPERATURE = materialModelDatabase->getDoubleWithDefault("Default_Temperature",310.0);
  default_BURNUP = materialModelDatabase->getDoubleWithDefault("Default_Burnup",0.0);
  default_OXYGEN_CONCENTRATION = materialModelDatabase->getDoubleWithDefault("Default_Oxygen_Concentration",0.0);

  AMP::Mesh::DOFMap::shared_ptr dof_map_0;
  dof_map_0 = mesh->getDOFMap( displacementVar );

  AMP::Mesh::DOFMap::shared_ptr dof_map_1;
  dof_map_1 = mesh->getDOFMap( temperatureVar );

  AMP::Mesh::MeshManager::Adapter::ElementIterator  el = mesh->beginElement();
  AMP::Mesh::MeshManager::Adapter::ElementIterator  end_el = mesh->endElement();

  for( ; el != end_el; ++el) {
    AMP::Mesh::MeshManager::Adapter::Element elem = *el;

    unsigned int num_local_type0Dofs = 0;
    std::vector<unsigned int> d_type0DofIndices[3];

    for(unsigned int i = 0; i < 3; i++) {
      dof_map_0->getDOFs(elem, d_type0DofIndices[i], i);
      num_local_type0Dofs += d_type0DofIndices[i].size();
    }

    std::vector<double> elementForceVector;
    elementForceVector.resize(num_local_type0Dofs);
    for(unsigned int i = 0; i < num_local_type0Dofs; i++)
      elementForceVector[i] = 0.0;

    unsigned int num_local_type1Dofs = 0;
    std::vector<unsigned int> d_type1DofIndices;

    dof_map_1->getDOFs(elem, d_type1DofIndices);
    num_local_type1Dofs = d_type1DofIndices.size();

    std::vector<double> currElementTemperatureVector;
    currElementTemperatureVector.resize(num_local_type1Dofs);

    std::vector<double> prevElementTemperatureVector;
    prevElementTemperatureVector.resize(num_local_type1Dofs);

    unsigned int d_numNodesForCurrentElement = elem.numNodes();
    for(unsigned int r = 0; r < d_numNodesForCurrentElement; r++) {
      currElementTemperatureVector[r] = currTemperatureVec->getValueByGlobalID( d_type1DofIndices[r] );
      prevElementTemperatureVector[r] = prevTemperatureVec->getValueByGlobalID( d_type1DofIndices[r] );
    }

    d_elem = &(elem.getElem());

    d_fe->reinit(d_elem);

    const unsigned int num_nodes = d_elem->n_nodes();

    for(unsigned int qp = 0; qp < d_qrule->n_points(); qp++) {
      double Bl_np1[6][24];

      for(int i = 0; i < 6; i++) {
        for(unsigned int j = 0; j < (3 * num_nodes); j++) {
          Bl_np1[i][j] = 0.0;
        }
      }

      for(unsigned int i = 0; i < num_nodes; i++) {
        Bl_np1[0][(3 * i) + 0] = dphi[i][qp](0);
        Bl_np1[1][(3 * i) + 1] = dphi[i][qp](1);
        Bl_np1[2][(3 * i) + 2] = dphi[i][qp](2);
        Bl_np1[3][(3 * i) + 1] = dphi[i][qp](2);
        Bl_np1[3][(3 * i) + 2] = dphi[i][qp](1);
        Bl_np1[4][(3 * i) + 0] = dphi[i][qp](2);
        Bl_np1[4][(3 * i) + 2] = dphi[i][qp](0);
        Bl_np1[5][(3 * i) + 0] = dphi[i][qp](1);
        Bl_np1[5][(3 * i) + 1] = dphi[i][qp](0);
      }

      double currTemperatureAtGaussPoint = 0.0, prevTemperatureAtGaussPoint = 0.0;
      for(unsigned int k = 0; k < num_nodes; k++) {
        currTemperatureAtGaussPoint += (currElementTemperatureVector[k] * phi[k][qp]);
        prevTemperatureAtGaussPoint += (prevElementTemperatureVector[k] * phi[k][qp]);
      }

      if(d_useMaterialsLibrary == true) {
        std::map<std::string, boost::shared_ptr<std::vector<double> > > inputMaterialParameters;

        std::string temperatureString = "temperature"; // in the future get from input file
        std::string burnupString = "burnup"; // in the future get from input file
        std::string oxygenString = "concentration"; // in the future get from input file

        boost::shared_ptr<std::vector<double> > tempVec(new std::vector<double> );      
        boost::shared_ptr<std::vector<double> > burnupVec(new std::vector<double> );      
        boost::shared_ptr<std::vector<double> > oxygenVec(new std::vector<double> );      

        inputMaterialParameters.insert( std::make_pair( temperatureString, tempVec) ); 
        inputMaterialParameters.insert( std::make_pair( burnupString, burnupVec) );
        inputMaterialParameters.insert( std::make_pair( oxygenString, oxygenVec) );

        tempVec->push_back(currTemperatureAtGaussPoint);
        burnupVec->push_back(default_BURNUP);
        oxygenVec->push_back(default_OXYGEN_CONCENTRATION);

        std::vector<double> YM(1);
        std::vector<double> PR(1);
        std::vector<double> TEC(1);

        std::string ymString = "YoungsModulus";
        std::string prString = "PoissonRatio";
        std::string tecString = "ThermalExpansion";

        d_material->property(ymString)->evalv(YM, inputMaterialParameters);
        d_material->property(prString)->evalv(PR, inputMaterialParameters);
        d_material->property(tecString)->evalv(TEC, inputMaterialParameters);

        youngsModulus = YM[0];
        poissonsRatio = PR[0];
        thermalExpansionCoefficient = TEC[0];
      }

      double d_thermalStress[6], d_thermalStrain[6], d_constitutiveMatrix[6][6];
      for(unsigned int i = 0; i < 6; i++) {
        d_thermalStrain[i] = 0.0;
        d_thermalStress[i] = 0.0;
        for(unsigned int j = 0; j < 6; j++)
          d_constitutiveMatrix[i][j] = 0.0;
      }

      double E = youngsModulus;
      double nu = poissonsRatio;
      double K = E / (3.0 * (1.0 - (2.0 * nu)));
      double G = E / (2.0 * (1.0 + nu));

      for(unsigned int i = 0; i < 3; i++)
        d_constitutiveMatrix[i][i] += (2.0 * G);

      for(unsigned int i = 3; i < 6; i++)
        d_constitutiveMatrix[i][i] += G;

      for(unsigned int i = 0; i < 3; i++) {
        for(unsigned int j = 0; j < 3; j++) {
          d_constitutiveMatrix[i][j] += (K - ((2.0 * G) / 3.0));
        }
      }

      for(unsigned int i = 0; i < 3; i++) {
        d_thermalStrain[i] = thermalExpansionCoefficient * (currTemperatureAtGaussPoint - prevTemperatureAtGaussPoint);
      }

      for(unsigned int i = 0; i < 6; i++) {
        for(unsigned int j = 0; j < 6; j++) {
          d_thermalStress[i] += (d_constitutiveMatrix[i][j] * d_thermalStrain[j]);
        }
      }

      if(d_iDebugPrintInfoLevel > 11) {
        for(unsigned int i = 0; i < 6; i++) {
          std::cout<<"d_thermalStrain["<<i<<"] = "<<d_thermalStrain[i]<<std::endl;
          std::cout<<"d_thermalStress["<<i<<"] = "<<d_thermalStress[i]<<std::endl;
        }
      }

      for(unsigned int j = 0; j < num_nodes; j++) {
        for(unsigned int d = 0; d < 3; d++) {

          double tmp = 0;
          for(unsigned int m = 0; m < 6; m++) {
            tmp += (Bl_np1[m][(3 * j) + d] * d_thermalStress[m]);
          }

          elementForceVector[(3*j) + d] += (JxW[qp] * tmp);
          if(d_iDebugPrintInfoLevel > 8) {
            AMP::pout<<"Force vector component: "<<elementForceVector[(3*j) + d]<<std::endl;
          }
        } // end of d
      } //end of j
    } // end of qp

    for(unsigned int r = 0; r < num_nodes; r++) {
      for(unsigned int d = 0; d < 3; d++) {
        if(d_iDebugPrintInfoLevel > 11) {
          std::cout<<"d_type0DofIndices["<<d<<"]["<<r<<"] = "<<d_type0DofIndices[d][r]<<" elementForceVector["<<(3*r)+d<<"] = "<<elementForceVector[(3*r)+d]<<std::endl;
        }
        rInternal->addValueByGlobalID( d_type0DofIndices[d][r], elementForceVector[(3*r) + d] );
      }
    }

  } // end of el

  rInternal->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_ADD );

}

*/






















