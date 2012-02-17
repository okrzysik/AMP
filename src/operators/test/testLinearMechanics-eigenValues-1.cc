#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>

#include "boost/shared_ptr.hpp"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"

#include "ampmesh/Mesh.h"

#include "mesh.h"
#include "libmesh.h"
#include "mesh_generation.h"

#include "operators/mechanics/IsotropicElasticModel.h"
#include "operators/mechanics/MechanicsLinearElement.h"
#include "operators/mechanics/MechanicsLinearFEOperator.h"

#include "operators/boundary/DirichletMatrixCorrection.h"
#include "operators/boundary/DirichletVectorCorrection.h"


void myTest(AMP::UnitTest *ut)
{
  std::string exeName("testLinearMechanics-eigenValues-1");
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
  int npes = globalComm.getSize();

  if(npes == 1) {
    boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
    AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
    input_db->printClassData(AMP::plog);

    AMP_INSIST(input_db->keyExists("OutputFileName"), "Key ''OutputFileName'' is missing!");
    std::string outFileName = input_db->getString("OutputFileName");

    FILE* fp; 
    fp = fopen(outFileName.c_str(), "w");
    fprintf(fp, "clc; \n clear; \n A = zeros(24, 24); \n \n");

    AMP_INSIST(input_db->keyExists("DISTORT_ELEMENT"), "Key ''DISTORT_ELEMENT'' is missing!");
    bool distortElement = input_db->getBool("DISTORT_ELEMENT");

    boost::shared_ptr< ::Mesh > mesh(new ::Mesh(3));
    MeshTools::Generation::build_cube((*(mesh.get())), 1, 1, 1, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, HEX8, false);

    if(distortElement) {
      ::Elem* elemPtr = mesh->elem(0);

      (elemPtr->point(0))(0) -= 0.1;
      (elemPtr->point(0))(1) -= 0.2;
      (elemPtr->point(0))(2) -= 0.3;

      (elemPtr->point(3))(0) -= 0.2;
      (elemPtr->point(3))(1) += 0.2;
      (elemPtr->point(3))(2) -= 0.1;

      (elemPtr->point(6))(0) += 0.3;
      (elemPtr->point(6))(1) += 0.2;
      (elemPtr->point(6))(2) += 0.1;
    }

    AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::MeshManager::Adapter::shared_ptr (
        new AMP::Mesh::MeshManager::Adapter (mesh) );

    AMP_INSIST(input_db->keyExists("Isotropic_Model"), "Key ''Isotropic_Model'' is missing!");
    boost::shared_ptr<AMP::Database> matModel_db = input_db->getDatabase("Isotropic_Model");
    boost::shared_ptr<AMP::Operator::MechanicsMaterialModelParameters> matModelParams(new
        AMP::Operator::MechanicsMaterialModelParameters( matModel_db ) );
    boost::shared_ptr<AMP::Operator::IsotropicElasticModel> isotropicModel (new AMP::Operator::IsotropicElasticModel( matModelParams));

    AMP_INSIST(input_db->keyExists("Mechanics_Linear_Element"), "Key ''Mechanics_Linear_Element'' is missing!");
    boost::shared_ptr<AMP::Database> elemOp_db = input_db->getDatabase("Mechanics_Linear_Element");
    boost::shared_ptr<AMP::Operator::ElementOperationParameters> elemOpParams (new AMP::Operator::ElementOperationParameters( elemOp_db ));
    boost::shared_ptr<AMP::Operator::MechanicsLinearElement> mechLinElem (new AMP::Operator::MechanicsLinearElement( elemOpParams ));

    AMP_INSIST(input_db->keyExists("Mechanics_Assembly"), "Key ''Mechanics_Assembly'' is missing!");
    boost::shared_ptr<AMP::Database> mechAssembly_db = input_db->getDatabase("Mechanics_Assembly");
    boost::shared_ptr<AMP::Operator::MechanicsLinearFEOperatorParameters> mechOpParams(new
        AMP::Operator::MechanicsLinearFEOperatorParameters( mechAssembly_db ));
    mechOpParams->d_materialModel = isotropicModel;
    mechOpParams->d_elemOp = mechLinElem;
    mechOpParams->d_MeshAdapter = meshAdapter;
    boost::shared_ptr<AMP::Operator::MechanicsLinearFEOperator> mechOp (new AMP::Operator::MechanicsLinearFEOperator( mechOpParams ));

    boost::shared_ptr<AMP::LinearAlgebra::Matrix> mechMat = mechOp->getMatrix();

    for(int i = 0; i < 24; i++) {
      std::vector<unsigned int> matCols;
      std::vector<double> matVals;
      mechMat->getRowByGlobalID(i, matCols, matVals);
      for(unsigned int j = 0; j < matCols.size(); j++) {
        fprintf(fp, "A(%d, %d) = %.15lf ; \n", (i + 1), (int)(matCols[j] + 1), matVals[j]);
      }//end for j
      fprintf(fp, "\n");
    }//end for i
   
    AMP::Mesh::DOFMap::shared_ptr dof_map = meshAdapter->getDOFMap(mechOp->getOutputVariable());
    std::vector<unsigned int> dofs(3);
    dofs[0] = 0;
    dofs[1] = 1;
    dofs[2] = 2;
    for(int i = 0 ; i < 8; i++) {
        AMP::Mesh::LibMeshNode nd = meshAdapter->getNode( i );
        fprintf(fp, "nd = %d, x = %.15lf, y = %.15lf, z = %.15lf \n", i, nd.x(), nd.y(), nd.z());
        std::vector<unsigned int> globalIds;
        dof_map->getDOFs(nd, globalIds, dofs);
        fprintf(fp, "nd = %d, d0 = %d, d1 = %d, d2 = %d \n", i, (int)globalIds[0], (int)globalIds[1], (int)globalIds[2]);
    }

    fprintf(fp, "format long e; \n\n");
    fprintf(fp, "sort(eig(A)) \n\n");
    fclose(fp);
  }

  ut->passes(exeName);

}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    try {
        myTest(&ut);
    } catch (std::exception &err) {
        std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
        ut.failure("ERROR: While testing");
    } catch( ... ) {
        std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
        ut.failure("ERROR: While testing");
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}   


