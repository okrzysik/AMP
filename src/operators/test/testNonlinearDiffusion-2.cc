#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "boost/shared_ptr.hpp"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"

#include "ampmesh/MeshManager.h"
#include "ampmesh/MeshVariable.h"

#include "libmesh.h"

#include "operators/diffusion/DiffusionTransportModel.h"
#include "operators/diffusion/DiffusionConstants.h"
#include "operators/diffusion/DiffusionLinearElement.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionLinearFEOperatorParameters.h"
#include "operators/diffusion/DiffusionNonlinearElement.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperatorParameters.h"
#include "../ElementPhysicsModelParameters.h"
#include "../ElementPhysicsModelFactory.h"
#include "../OperatorBuilder.h"

#include "materials/Material.h"

#include "patchfunctions.h"


/**
 * This test is patch test for the diffusion operator.
 */
void nonlinearTest(AMP::UnitTest *ut, std::string exeName,
        double function(const double, const double, const double))
{
  // Initialization
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm = AMP::AMP_MPI(AMP_COMM_WORLD);

  std::cout << "testing with input file " << input_file << std::endl;
  std::cout.flush();

  // Test create
  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  std::string mesh_file = input_db->getString("Mesh");

  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = AMP::Mesh::MeshManager::Adapter::shared_ptr ( new AMP::Mesh::MeshManager::Adapter () );
  meshAdapter->readExodusIIFile ( mesh_file.c_str() );

  boost::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> diffOp;
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel;
  boost::shared_ptr<AMP::InputDatabase> diffFEOp_db =
    boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("NonlinearDiffusionOp"));
  boost::shared_ptr<AMP::Operator::Operator> nonlinearOperator = AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
														"NonlinearDiffusionOp",
														input_db,
														elementModel);
  diffOp = boost::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(nonlinearOperator);

  ut->passes(exeName+": create");
  std::cout.flush();

  // set up defaults for materials arguments and create transport model
  boost::shared_ptr<AMP::Database> transportModel_db;
  if (input_db->keyExists("DiffusionTransportModel"))
	  transportModel_db = input_db->getDatabase("DiffusionTransportModel");
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel =
          AMP::Operator::ElementPhysicsModelFactory::createElementPhysicsModel(transportModel_db);
  boost::shared_ptr<AMP::Operator::DiffusionTransportModel> transportModel =
          boost::dynamic_pointer_cast<AMP::Operator::DiffusionTransportModel>(elementPhysicsModel);

  double defTemp = transportModel_db->getDoubleWithDefault("Default_Temperature", 400.0);
  double defConc = transportModel_db->getDoubleWithDefault("Default_Concentration", .33);
  double defBurn = transportModel_db->getDoubleWithDefault("Default_Burnup", .5);

  // create parameters
  boost::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperatorParameters> diffOpParams(new
  AMP::Operator::DiffusionNonlinearFEOperatorParameters( diffFEOp_db ));

  // nullify vectors in parameters
  diffOpParams->d_temperature .reset();
  diffOpParams->d_concentration.reset();
  diffOpParams->d_burnup.reset();

  // create vectors for parameters
  boost::shared_ptr<AMP::Database> active_db = diffFEOp_db->getDatabase("ActiveInputVariables");
  AMP::LinearAlgebra::Variable::shared_ptr tVar(new AMP::Mesh::NodalScalarVariable(
          active_db->getStringWithDefault("Temperature","not_specified")));
  AMP::LinearAlgebra::Variable::shared_ptr cVar(new AMP::Mesh::NodalScalarVariable(
          active_db->getStringWithDefault("Concentration","not_specified")));
  AMP::LinearAlgebra::Variable::shared_ptr bVar(new AMP::Mesh::NodalScalarVariable(
          active_db->getStringWithDefault("Burnup","not_specified")));
  AMP::LinearAlgebra::Vector::shared_ptr tVec = meshAdapter->createVector( tVar );
  AMP::LinearAlgebra::Vector::shared_ptr cVec = meshAdapter->createVector( cVar );
  AMP::LinearAlgebra::Vector::shared_ptr bVec = meshAdapter->createVector( bVar );
  tVec->setToScalar(defTemp);
  cVec->setToScalar(defConc);
  bVec->setToScalar(defBurn);

  // set principal variable vector
  if (diffOp->getPrincipalVariableId() == AMP::Operator::Diffusion::TEMPERATURE)
      diffOpParams->d_temperature = tVec;
  if (diffOp->getPrincipalVariableId() == AMP::Operator::Diffusion::CONCENTRATION)
      diffOpParams->d_concentration = cVec;
  if (diffOp->getPrincipalVariableId() == AMP::Operator::Diffusion::BURNUP)
      diffOpParams->d_burnup = bVec;

  // set frozen vectors in parameters
  if (diffFEOp_db->getBoolWithDefault("FreezeTemperature",false))
    diffOpParams->d_FrozenTemperature = tVec;
  if (diffFEOp_db->getBoolWithDefault("FreezeConcentration",false))
    diffOpParams->d_FrozenConcentration = cVec;
  if (diffFEOp_db->getBoolWithDefault("FreezeBurnup",false))
    diffOpParams->d_FrozenBurnup = bVec;

  // set transport model
  diffOpParams->d_transportModel = transportModel;

  // Reset with new parameters
  diffOp->reset(diffOpParams);

  // set  up variables for apply tests
  AMP::LinearAlgebra::Variable::shared_ptr diffSolVar = diffOp->getInputVariable(diffOp->getPrincipalVariableId());
  AMP::LinearAlgebra::Variable::shared_ptr diffRhsVar = diffOp->getOutputVariable();
  AMP::LinearAlgebra::Variable::shared_ptr diffResVar = diffOp->getOutputVariable();
  AMP::LinearAlgebra::Variable::shared_ptr workVar(new AMP::Mesh::NodalScalarVariable("work"));
  std::vector<unsigned int> nonPrincIds = diffOp->getNonPrincipalVariableIds();
  unsigned int numNonPrincIds = nonPrincIds.size();
  std::vector<AMP::LinearAlgebra::Variable::shared_ptr> nonPrincVars(numNonPrincIds);
  for (size_t i=0; i<numNonPrincIds; i++) {
      nonPrincVars[i] = diffOp->getInputVariable(nonPrincIds[i]);
  }

  // set up vectors for apply tests
  //std::string msgPrefix=exeName+": apply";
  AMP::LinearAlgebra::Vector::shared_ptr diffSolVec = meshAdapter->createVector( diffSolVar );
  AMP::LinearAlgebra::Vector::shared_ptr diffRhsVec = meshAdapter->createVector( diffRhsVar );
  AMP::LinearAlgebra::Vector::shared_ptr diffResVec = meshAdapter->createVector( diffResVar );
  std::vector<AMP::LinearAlgebra::Vector::shared_ptr> nonPrincVecs(numNonPrincIds);
  for (unsigned int i=0; i<numNonPrincIds; i++) {
      nonPrincVecs[i] = meshAdapter->createVector( nonPrincVars[i] );
      if (nonPrincIds[i] == AMP::Operator::Diffusion::TEMPERATURE) nonPrincVecs[i]->setToScalar(defTemp);
      if (nonPrincIds[i] == AMP::Operator::Diffusion::CONCENTRATION) nonPrincVecs[i]->setToScalar(defConc);
      if (nonPrincIds[i] == AMP::Operator::Diffusion::BURNUP) nonPrincVecs[i]->setToScalar(defBurn);
  }
  diffRhsVec->setToScalar(0.0);

  for (AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator curNode =
              meshAdapter->beginOwnedNode();
          curNode != meshAdapter->endOwnedNode(); curNode++)
  {
      double x = curNode->x();
      double y = curNode->y();
      double z = curNode->z();
      size_t i = curNode->globalID();
      double fval = function(x,y,z);
      diffSolVec->setValueByGlobalID(i, fval);
  }

  // Compute finite element operator
  diffOp->apply(diffRhsVec, diffSolVec, diffResVec, 1.,0.);

  // Check that interior values are zero.
  double totalBnd = 0.;
  for (size_t face = 1; face<64; face ++) {
      for (AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator curNode =
                  meshAdapter->beginOwnedBoundary(face);
              curNode != meshAdapter->endOwnedBoundary(face); curNode++ )
      {
          size_t i = curNode->globalID();
          double fval = diffResVec->getValueByGlobalID(i);
          totalBnd += fabs(fval);
      }
  }
  double totalAll = 0.;
  for (AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator curNode =
              meshAdapter->beginOwnedNode();
          curNode != meshAdapter->endOwnedNode(); curNode++)
  {
      size_t i = curNode->globalID();
      double fval = diffResVec->getValueByGlobalID(i);
      totalAll += fabs(fval);
  }
  totalBnd = globalComm.sumReduce(totalBnd);
  totalAll = globalComm.sumReduce(totalAll);
  int rank = globalComm.getRank();
  if (rank == 0) {
      std::cout << "******** All = " << totalAll <<", Bnd = " << totalBnd << ", Err = "
          << totalAll-totalBnd << "********" << std::endl;
  }

  // write values in mathematica form
  int nranks = globalComm.getSize();
  if (nranks == 1) {
      size_t  nnodes = meshAdapter->numLocalNodes();
      int proc = globalComm.getRank();
      int nproc = globalComm.getSize();
      std::string filename = "values-"+exeName;
      std::ofstream file(filename.c_str());
      if (proc == 0) {
          file << "values={"<<"\n";
      }
      for (size_t i=0; i<nnodes; i++)
      {
          AMP::Mesh::MeshManager::Adapter::Node node = meshAdapter->getNode(i);
          double x = node.x();
          double y = node.y();
          double z = node.z();

          int ii = i;
          double rval = diffResVec->getValueByLocalID(ii);
          double fval = function(x,y,z);
          file << "{"<<x<<","<<y<<","<<z<<","<<rval<<","<<fval<<"}";
          if (i<nnodes-1) file <<",\n";
      }
      if (proc < nproc-1) {
          file <<",\n";
      } else {
          file <<"}\n";
      }

      /* view with the following Mathematica commands:
       *
(* sed -e 's/e\([+-]\)/10.*^\1/g' file_name > values2 *)
dir = "W:\\amp\\code43\\trunk\\build\\debug\\src\\operators\\test";
SetDirectory[dir];
ReadList["values2"];
tval = Transpose[values];
pts = Point /@ Transpose[Take[tval, {1, 3}]];
vals = tval[[4]];
funs = tval[[5]];
svals = (vals - Min[vals]) / (Max[vals] - Min[vals]);
sfuns = (funs - Min[funs]) / (Max[funs] - Min[funs]);
hvals = Hue[#, 1., 1.] & /@ svals;
hfuns = Hue[#, 1., 1.] & /@ sfuns;
gvals = Graphics3D@Flatten[Transpose[{hvals, pts}]];
gfuns = Graphics3D@Flatten[Transpose[{hfuns, pts}]];
valuesbnd = Select[values, Abs[#[[4]]] > .00000001 &]; tvalbnd =
Transpose[Take[ Transpose[valuesbnd], {1, 3}]];

Show[gvals, Axes -> True, AxesLabel -> {"x", "y", "z"}]

Show[gfuns, Axes -> True, AxesLabel -> {"x", "y", "z"}]

Show[Graphics3D[Point /@ tvalbnd], AspectRatio -> 1]
       */
  }

  ut->passes("values-"+exeName);
}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;
     const int NUMFILES=1;
    std::string files[NUMFILES] = {
        "Diffusion-TUI-Thermal-1"
    };

    try {
        for (int i=0; i<NUMFILES; i++)
            nonlinearTest(&ut, files[i], x_linear);
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



