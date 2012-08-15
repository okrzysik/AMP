#include <iostream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include "mpi.h"

#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "ampmesh/Mesh.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/Variable.h"
#include "vectors/VectorBuilder.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/VectorSelector.h"
#include "operators/boundary/NeumannVectorCorrection.h"

void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  boost::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
  boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
  params->setComm(globalComm);
  boost::shared_ptr<AMP::Mesh::Mesh> mesh = AMP::Mesh::Mesh::buildMesh(params);

  AMP::LinearAlgebra::Variable::shared_ptr var(new  AMP::LinearAlgebra::Variable("myVar")); 
  AMP::Discretization::DOFManager::shared_ptr nodalScalarDOF = AMP::Discretization::simpleDOFManager::create(mesh,AMP::Mesh::Vertex,1,1,true);  
  AMP::Discretization::DOFManager::shared_ptr nodalVectorDOF = AMP::Discretization::simpleDOFManager::create(mesh,AMP::Mesh::Vertex,1,2,true);  

  AMP::LinearAlgebra::Vector::shared_ptr scalarVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, var, true );
  AMP::LinearAlgebra::Vector::shared_ptr vectorVec = AMP::LinearAlgebra::createVector( nodalVectorDOF, var, true );

  scalarVec->zero();
  vectorVec->zero();

  boost::shared_ptr<AMP::Database>  bnd_db = input_db->getDatabase("NeumannVectorCorrection1");
  boost::shared_ptr<AMP::Operator::NeumannVectorCorrectionParameters> vectorCorrectionParameters (new AMP::Operator::NeumannVectorCorrectionParameters( bnd_db ) );
  vectorCorrectionParameters->d_variable = var;
  vectorCorrectionParameters->d_Mesh = mesh;
  boost::shared_ptr<AMP::Operator::NeumannVectorCorrection> neumannBndOp (new AMP::Operator::NeumannVectorCorrection(vectorCorrectionParameters));

  neumannBndOp->addRHScorrection(scalarVec);
  AMP::pout << scalarVec << std::endl;

  bnd_db = input_db->getDatabase("NeumannVectorCorrection2");
  vectorCorrectionParameters.reset(new AMP::Operator::NeumannVectorCorrectionParameters( bnd_db ) );
  vectorCorrectionParameters->d_variable = var;
  vectorCorrectionParameters->d_Mesh = mesh;
  neumannBndOp.reset(new AMP::Operator::NeumannVectorCorrection(vectorCorrectionParameters));

  neumannBndOp->addRHScorrection(vectorVec);
  //AMP::LinearAlgebra::Vector::shared_ptr vectorVec1 = vectorVec->select( AMP::LinearAlgebra::VS_Stride("U", 0, 2) , "V1" );
  AMP::LinearAlgebra::Vector::shared_ptr vectorVec2 = vectorVec->select( AMP::LinearAlgebra::VS_Stride("V", 1, 2) , "V2" );

  AMP::pout << vectorVec2 << std::endl;

}



int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::string exeName = "testNeumannVectorCorrection";

  try {
    myTest(&ut, exeName);
    ut.passes(exeName);
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
