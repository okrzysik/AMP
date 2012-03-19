#include <string>
#include "utils/UnitTest.h"
#include "utils/AMPManager.h"
#include "utils/InputManager.h"

#include "ampmesh/SiloIO.h"
#include "ampmesh/Mesh.h"
#include "vectors/VectorBuilder.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"

#include <vector>
#include <fstream>
#include <algorithm>
#include <iostream>
#include "utils/Utilities.h"

#include "boost/shared_ptr.hpp"
#include "operators/VolumeIntegralOperator.h"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/AMP_MPI.h"
#include "utils/PIO.h"

#include "operators/BVPOperatorParameters.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/DiffusionNonlinearFEOperator.h"
#include "operators/ColumnOperator.h"
#include "operators/CoupledOperator.h"
#include "operators/ComputeSurfaceNormal.h"
#include "operators/OperatorBuilder.h"
#include "operators/NeumannVectorCorrection.h"
#include "operators/ColumnBoundaryOperator.h"
#include "operators/boundary/RobinVectorCorrection.h"

#include "operators/map/ScalarZAxisMap.h"
#include "operators/map/AsyncMapColumnOperator.h"

#include "solvers/ColumnSolver.h"
#include "solvers/PetscKrylovSolverParameters.h"
#include "solvers/PetscKrylovSolver.h"
#include "solvers/PetscSNESSolverParameters.h"
#include "solvers/PetscSNESSolver.h"
#include "solvers/TrilinosMLSolver.h"

/* Libmesh files */
#include "fe_type.h"
#include "fe_base.h"
#include "elem.h"
#include "quadrature.h"

#include "enum_order.h"
#include "enum_fe_family.h"
#include "enum_quadrature_type.h"
#include "auto_ptr.h"
#include "string_to_enum.h"

#include "cell_hex8.h"
#include "node.h"

#define _PI_ 3.14159265
#define __INIT_FN__(x,y) (800+((0.00004-20*pow(x,2)-20*pow(y,2))*pow(10,6))) // Manufactured Solution
#define __FsnK__() (80000000)


void calculateManufacturedSolution(AMP::Mesh::Mesh::shared_ptr meshAdapter,
    AMP::LinearAlgebra::Vector::shared_ptr manufacturedSolution)
{
  //------------------------------------------
  //  CALCULATE THE MANFACTURED SOLUTION //
  //------------------------------------------
                
    AMP::Discretization::DOFManager::shared_ptr dof_map = manufacturedSolution->getDOFManager();

    AMP::Mesh::MeshIterator el = meshAdapter->getIterator(AMP::Mesh::Volume, 0);
    AMP::Mesh::MeshIterator end_el = el.end();

    for( ; el != end_el; ++el) {
      std::vector<AMP::Mesh::MeshElement> d_currNodes = el->getElements(AMP::Mesh::Vertex);

      std::vector<AMP::Mesh::MeshElementID> globalIDs(d_currNodes.size()); 

      std::vector<size_t> d_dofIndices; 
      for(unsigned int j = 0; j < d_currNodes.size(); j++) {
        globalIDs[j] = d_currNodes[j].globalID();
      }  // end of j
      dof_map->getDOFs(globalIDs, d_dofIndices);

      for(unsigned int j = 0; j < d_currNodes.size(); j++) {
        std::vector<double> pt = d_currNodes[j].coord();
        double val1 = __INIT_FN__(pt[0],pt[1]); 
        manufacturedSolution->setLocalValueByGlobalID(d_dofIndices[j],val1);
      } //end for node
    }
    manufacturedSolution->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );

}

void calculateSources(AMP::Mesh::Mesh::shared_ptr meshAdapter,
    AMP::Discretization::DOFManager::shared_ptr gaussPointDOF,
    boost::shared_ptr<AMP::LinearAlgebra::Vector> manufacturedRHS)
{
  // Compute the source on the gauss point

  AMP::Mesh::MeshIterator el     = meshAdapter->getIterator(AMP::Mesh::Volume, 0);
  AMP::Mesh::MeshIterator end_el = el.end();

  libMeshEnums::Order feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>("FIRST");
  libMeshEnums::FEFamily feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>("LAGRANGE");

  boost::shared_ptr < ::FEType > d_feType ( new ::FEType(feTypeOrder, feFamily) );
  boost::shared_ptr < ::FEBase > d_fe ( (::FEBase::build(3, (*d_feType))).release() );

  libMeshEnums::Order qruleOrder = Utility::string_to_enum<libMeshEnums::Order>("SECOND");
  boost::shared_ptr < ::QBase > d_qrule ( (::QBase::build("QGAUSS", 3, qruleOrder)).release() );

  d_fe->attach_quadrature_rule( d_qrule.get() );

  for( ; el != end_el; ++el) {

    std::vector<size_t> d_gaussPtIndices; 
    gaussPointDOF->getDOFs ( el->globalID(), d_gaussPtIndices);

    std::vector<Point> coordinates = d_fe->get_xyz();
    for (unsigned int qp = 0; qp < d_gaussPtIndices.size(); qp++) {
      double manufacturedAtGauss1;
      manufacturedAtGauss1 = __FsnK__(); 

      manufacturedRHS->setValueByGlobalID( d_gaussPtIndices[qp], manufacturedAtGauss1 );
    }

  }

  manufacturedRHS->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
}

void computeL2Norm(AMP::Mesh::Mesh::shared_ptr meshAdapter,
    AMP::AMP_MPI globalComm,
    AMP::LinearAlgebra::Vector::shared_ptr  TemperatureVec,
    double * discretizationErrorNorm2)
{
  //------------------------------------------
  // CALCULATE THE L2Norm OF (U-Uh)         //
  //------------------------------------------
  AMP::Mesh::MeshIterator el     = meshAdapter->getIterator(AMP::Mesh::Volume, 0);
  AMP::Mesh::MeshIterator end_el = el.end();

  AMP::Discretization::DOFManager::shared_ptr dof_map = TemperatureVec->getDOFManager();

  libMeshEnums::Order feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>("FIRST");
  libMeshEnums::FEFamily feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>("LAGRANGE");

  boost::shared_ptr < ::FEType > d_feType ( new ::FEType(feTypeOrder, feFamily) );
  boost::shared_ptr < ::FEBase > d_fe ( (::FEBase::build(3, (*d_feType))).release() );

  libMeshEnums::Order qruleOrder = Utility::string_to_enum<libMeshEnums::Order>("SECOND");
  boost::shared_ptr < ::QBase > d_qrule ( (::QBase::build("QGAUSS", 3, qruleOrder)).release() );

  d_fe->attach_quadrature_rule( d_qrule.get() );

  for( ; el != end_el; ++el) {

    std::vector<AMP::Mesh::MeshElement> d_currNodes = el->getElements(AMP::Mesh::Vertex);

    std::vector<size_t> bndGlobalIds;
    std::vector<AMP::Mesh::MeshElementID> globalIDs(d_currNodes.size()); 
    for(unsigned int j = 0; j < d_currNodes.size(); j++) {
      globalIDs[j] = d_currNodes[j].globalID();
    }  // end of j
    dof_map->getDOFs(globalIDs, bndGlobalIds);

    ::Elem* d_currElemPtr ( new ::Hex8);
    for(size_t j = 0; j < d_currNodes.size(); j++) {
      std::vector<double> pt = d_currNodes[j].coord();
      d_currElemPtr->set_node(j) = new ::Node(pt[0], pt[1], pt[2], j);
    }//end for j

    d_fe->reinit(d_currElemPtr);

    const std::vector<Real> & JxW = d_fe->get_JxW();
    std::vector<Point> coordinates = d_fe->get_xyz();
    const std::vector<std::vector<Real> > & phi = d_fe->get_phi();

    std::vector<double>  computedAtGauss(d_qrule->n_points(), 0.0);
    for (unsigned int qp = 0; qp < d_qrule->n_points(); qp++){
      for (unsigned int j = 0; j < bndGlobalIds.size(); j++){
        double computedAtNode = TemperatureVec->getValueByGlobalID(bndGlobalIds[j]); 
        computedAtGauss[qp]     += computedAtNode * phi[j][qp];
      }//end for j
    }//end for qp

    for (unsigned int qp = 0; qp < d_qrule->n_points(); qp++) {
      double px = coordinates[qp](0);
      double py = coordinates[qp](1);
      double manufacturedAtGauss = __INIT_FN__(px,py); 

      *discretizationErrorNorm2 += JxW[qp]* pow( (computedAtGauss[qp]
            - manufacturedAtGauss), 2);
    }
  }
  *discretizationErrorNorm2 = globalComm.sumReduce(*discretizationErrorNorm2); 
}

  void
createThermalOperators(boost::shared_ptr<AMP::InputDatabase> global_input_db,
    AMP::Mesh::Mesh::shared_ptr manager,
    boost::shared_ptr<AMP::Operator::ColumnOperator> &nonlinearColumnOperator,
    boost::shared_ptr<AMP::Operator::ColumnOperator> &linearColumnOperator)
{
  AMP::pout << "Entering createThermalOperators" << std::endl;

  boost::shared_ptr<AMP::Operator::OperatorParameters> emptyParams;
  nonlinearColumnOperator.reset(new AMP::Operator::ColumnOperator(emptyParams));
  linearColumnOperator.reset(new AMP::Operator::ColumnOperator(emptyParams));

  AMP::Mesh::Mesh::shared_ptr  meshAdapter1 = manager->Subset( "LeftMesh" );
  AMP::Mesh::Mesh::shared_ptr  meshAdapter2 = manager->Subset( "RightMesh" );
  //-----------------------------------------------
  //   CREATE THE NONLINEAR THERMAL OPERATOR 1 ----
  //-----------------------------------------------
  AMP_INSIST( global_input_db->keyExists("LeftNonlinearThermalOperator"), "key missing!" );
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> thermalNonlinearOperator = boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator( meshAdapter1 , "LeftNonlinearThermalOperator", global_input_db,thermalTransportModel));
  nonlinearColumnOperator->append(thermalNonlinearOperator);

  //-------------------------------------
  //   CREATE THE LINEAR THERMAL OPERATOR 1 ----
  //-------------------------------------
  AMP_INSIST( global_input_db->keyExists("LeftLinearThermalOperator"), "key missing!" );
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> thermalLinearOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter1 , "LeftLinearThermalOperator",global_input_db, thermalTransportModel));
  linearColumnOperator->append(thermalLinearOperator);

  //-----------------------------------------------
  //   CREATE THE NONLINEAR THERMAL OPERATOR 2 ----
  //-----------------------------------------------
  AMP_INSIST( global_input_db->keyExists("LeftNonlinearThermalOperator"), "key missing!" );
  thermalNonlinearOperator = boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator( meshAdapter2 , "RightNonlinearThermalOperator", global_input_db,thermalTransportModel));
  nonlinearColumnOperator->append(thermalNonlinearOperator);

  //-------------------------------------
  //   CREATE THE LINEAR THERMAL OPERATOR 2 ----
  //-------------------------------------
  AMP_INSIST( global_input_db->keyExists("LeftLinearThermalOperator"), "key missing!" );
  thermalLinearOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter2 , "RightLinearThermalOperator",global_input_db, thermalTransportModel));
  linearColumnOperator->append(thermalLinearOperator);

  AMP_ASSERT(nonlinearColumnOperator!=NULL);
  AMP_ASSERT(linearColumnOperator!=NULL);

  AMP::pout << "Leaving createThermalOperators" << std::endl;
}

void
createThermalSolvers(boost::shared_ptr<AMP::InputDatabase> &global_input_db,
		     AMP::LinearAlgebra::Vector::shared_ptr &globalSolVec,
		     boost::shared_ptr<AMP::Operator::Operator> &nonlinearOperator,
		     boost::shared_ptr<AMP::Operator::Operator> &linearOperator,
		     boost::shared_ptr<AMP::Solver::PetscSNESSolver>  &nonlinearSolver,
		     boost::shared_ptr<AMP::Solver::PetscKrylovSolver>  &linearSolver)
{
  //----------------------------------------------------------------//
  // initialize the nonlinear solver
  boost::shared_ptr<AMP::Database> nonlinearSolver_db = global_input_db->getDatabase("NonlinearThermalSolver");
  boost::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(new AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db));

  // change the next line to get the correct communicator out
  nonlinearSolverParams->d_comm          = AMP::AMP_MPI(AMP_COMM_WORLD);
  nonlinearSolverParams->d_pOperator     = nonlinearOperator;
  nonlinearSolverParams->d_pInitialGuess = globalSolVec ;
  nonlinearSolver.reset(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));

  //-------------------------------------------------------------------------//
  // initialize the column preconditioner which is a diagonal block preconditioner

  boost::shared_ptr<AMP::Database> linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver");
  boost::shared_ptr<AMP::Database> columnPreconditioner_db = linearSolver_db->getDatabase("Preconditioner");
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> columnPreconditionerParams(new AMP::Solver::SolverStrategyParameters(columnPreconditioner_db));
  boost::shared_ptr<AMP::Operator::ColumnOperator> linearColumnOperator = boost::dynamic_pointer_cast<AMP::Operator::ColumnOperator>(linearOperator);
  AMP_ASSERT(linearColumnOperator);
  columnPreconditionerParams->d_pOperator = linearColumnOperator;
  boost::shared_ptr<AMP::Solver::ColumnSolver> columnPreconditioner(new AMP::Solver::ColumnSolver(columnPreconditionerParams));

  //-------------------------------------------------------------------------//

  boost::shared_ptr<AMP::Database> trilinosPreconditioner_db = columnPreconditioner_db->getDatabase("TrilinosPreconditioner");
  for(int id = 0; id != linearColumnOperator->getNumberOfOperators(); id++)
  {
    boost::shared_ptr<AMP::Solver::SolverStrategyParameters> trilinosPreconditionerParams(new AMP::Solver::SolverStrategyParameters(trilinosPreconditioner_db));
    trilinosPreconditionerParams->d_pOperator = linearColumnOperator->getOperator(id);
    boost::shared_ptr<AMP::Solver::TrilinosMLSolver> trilinosPreconditioner(new AMP::Solver::TrilinosMLSolver(trilinosPreconditionerParams));
    columnPreconditioner->append(trilinosPreconditioner);
  }

  //--------------------------------------------------------------------//
  // register the preconditioner with the Jacobian free Krylov solver
  linearSolver = nonlinearSolver->getKrylovSolver();
  linearSolver->setPreconditioner(columnPreconditioner);
}

void
createThermalMaps( boost::shared_ptr<AMP::InputDatabase> input_db,
		   AMP::Mesh::Mesh::shared_ptr manager,
                   AMP::LinearAlgebra::Vector::shared_ptr &thermalMapVec ,
		   boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator> &mapsColumn)
{
    boost::shared_ptr<AMP::Database> map_db = input_db->getDatabase( "MeshToMeshMaps" );

    AMP::Discretization::DOFManager::shared_ptr dof_map = thermalMapVec->getDOFManager();

    mapsColumn = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::ScalarZAxisMap> ( manager, dof_map , map_db );
    mapsColumn->setVector ( thermalMapVec );
}

void
registerMapswithThermalOperator( boost::shared_ptr<AMP::InputDatabase> input_db , 
		     boost::shared_ptr<AMP::Operator::ColumnOperator> &nonlinearThermalColumnOperator,
                     AMP::LinearAlgebra::Vector::shared_ptr &thermalMapVec )
{
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator>  curBVPop = boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator> ( nonlinearThermalColumnOperator->getOperator ( 0 ) );
  boost::shared_ptr<AMP::Operator::ColumnBoundaryOperator>  curBCcol = boost::dynamic_pointer_cast<AMP::Operator::ColumnBoundaryOperator> ( curBVPop->getBoundaryOperator () );
  boost::shared_ptr<AMP::Database> operator_db = input_db->getDatabase ( "LeftNonlinearThermalOperator" );
  boost::shared_ptr<AMP::Database>  curBCdb = input_db->getDatabase(operator_db->getString ( "BoundaryOperator" ));
  std::vector<std::string>  opNames = curBCdb->getStringArray ( "boundaryOperators" );
  boost::shared_ptr<AMP::Operator::RobinVectorCorrection>  gapBC = boost::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection> ( curBCcol->getBoundaryOperator ( 0 ) );
  gapBC->setVariableFlux ( thermalMapVec );
  gapBC->reset ( gapBC->getParameters() );

  curBVPop = boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator> ( nonlinearThermalColumnOperator->getOperator ( 0 ) );
  curBCcol = boost::dynamic_pointer_cast<AMP::Operator::ColumnBoundaryOperator> ( curBVPop->getBoundaryOperator () );
  operator_db = input_db->getDatabase ( "RightNonlinearThermalOperator" );
  curBCdb = input_db->getDatabase(operator_db->getString ( "BoundaryOperator" ));
  opNames = curBCdb->getStringArray ( "boundaryOperators" );
  gapBC = boost::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection> ( curBCcol->getBoundaryOperator ( 0 ) );
  gapBC->setVariableFlux ( thermalMapVec );
  gapBC->reset ( gapBC->getParameters() );
}



///////////////////////////////////////////////
//       Main Program     //
///////////////////////////////////////////////

void myTest(AMP::UnitTest *ut, std::string exeName)
{
  std::string input_file = "input_" + exeName;

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
  AMP::InputManager::getManager()->parseInputFile(input_file,input_db);
  input_db->printClassData(AMP::plog);

  AMP_INSIST(input_db->keyExists("NumberOfMeshes"), "Key does not exist");
  int numMeshes = input_db->getInteger("NumberOfMeshes");

  AMP::pout<<"Num meshes = "<<numMeshes<<std::endl;

  AMP::Mesh::MeshParameters::shared_ptr meshmgrParams (new AMP::Mesh::MeshParameters ( input_db ) );
  boost::shared_ptr<AMP::Mesh::Mesh> manager = AMP::Mesh::Mesh::buildMesh(meshmgrParams);

  //------------------------------------------
  //  CREATE THE THERMAL OPERATOR  //
  //------------------------------------------
  // create the nonlinear and linear thermal operators
  boost::shared_ptr<AMP::Operator::ColumnOperator> nonlinearThermalColumnOperator;
  boost::shared_ptr<AMP::Operator::ColumnOperator> linearThermalColumnOperator;

  createThermalOperators(input_db, manager, nonlinearThermalColumnOperator, linearThermalColumnOperator);

  AMP_ASSERT(nonlinearThermalColumnOperator!=NULL);
  AMP_ASSERT(linearThermalColumnOperator!=NULL);

  AMP::LinearAlgebra::Variable::shared_ptr outputVar = nonlinearThermalColumnOperator->getOutputVariable();

  AMP::Discretization::DOFManager::shared_ptr nodalScalarDOF = AMP::Discretization::simpleDOFManager::create(manager,AMP::Mesh::Vertex,1,1,true);
  //  create solution, rhs, and  residual vectors
  AMP::LinearAlgebra::Vector::shared_ptr  TemperatureVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, outputVar, true );
  AMP::LinearAlgebra::Vector::shared_ptr  ResidualVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, outputVar, true );

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;

  AMP::LinearAlgebra::Vector::shared_ptr manufacturedSolution = AMP::LinearAlgebra::createVector( nodalScalarDOF, outputVar, true );
  AMP::LinearAlgebra::Vector::shared_ptr thermalConductivity  = AMP::LinearAlgebra::createVector( nodalScalarDOF, outputVar, true );
  AMP::LinearAlgebra::Vector::shared_ptr manufacturedNormalGradient = AMP::LinearAlgebra::createVector( nodalScalarDOF, outputVar, true );
  AMP::LinearAlgebra::Vector::shared_ptr solutionError        = AMP::LinearAlgebra::createVector( nodalScalarDOF, outputVar, true );
  AMP::LinearAlgebra::Vector::shared_ptr thermMapVec          = AMP::LinearAlgebra::createVector( nodalScalarDOF, outputVar, true );

  boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator> mapsColumn;
  createThermalMaps( input_db,  manager, thermMapVec , mapsColumn);

  registerMapswithThermalOperator( input_db , nonlinearThermalColumnOperator , thermMapVec );

  int DOFsPerElement = 8;
  AMP::Discretization::DOFManager::shared_ptr gaussPointDOF = AMP::Discretization::simpleDOFManager::create(manager,AMP::Mesh::Volume,1,DOFsPerElement,true);

  AMP::pout<<"Creating gauss Vectors "<<std::endl;

  AMP::LinearAlgebra::Variable::shared_ptr manuSourceVar (new AMP::LinearAlgebra::Variable("SpecificPowerInWattsPerGram"));
  AMP::LinearAlgebra::Vector::shared_ptr manufacturedRHS = AMP::LinearAlgebra::createVector( gaussPointDOF, manuSourceVar , true );

  AMP::pout<<"Calculating Manufactured Solution and Sources "<<std::endl;

  calculateManufacturedSolution( manager, manufacturedSolution );

  calculateSources( manager , gaussPointDOF , manufacturedRHS );

  //------------------------------------------
#ifdef USE_SILO
  AMP::Mesh::SiloIO::shared_ptr  siloWriter( new AMP::Mesh::SiloIO);
  siloWriter->registerMesh(manager);

  siloWriter->registerVector( manufacturedSolution , manager, AMP::Mesh::Vertex , "ManufacturedSolution" );
  siloWriter->registerVector( TemperatureVec , manager, AMP::Mesh::Vertex , "ComputedSolution" );
  siloWriter->registerVector( manufacturedRHS , manager, AMP::Mesh::Vertex ,"ManufacturedRhs");
  siloWriter->writeFile( exeName, 0 );
#endif 

  TemperatureVec->copyVector(manufacturedSolution);
  std::cout << "Max value of manufactured solution : "<< manufacturedSolution->max() << " thermal conductivity : "<< thermalConductivity->max() << std::endl;
  std::cout << "Min value of manufactured solution : "<< manufacturedSolution->min() << " thermal conductivity : "<< thermalConductivity->min() << std::endl;
  std::cout << "Max value of manufactured RHS : "<< manufacturedRHS->max() << " normal gradient: "<< manufacturedNormalGradient->max() << std::endl;
  std::cout << "Min value of manufactured RHS : "<< manufacturedRHS->min() << " normal gradient: "<< manufacturedNormalGradient->min() << std::endl;

  //------------------------------------------
  // OPERATOR APPLY TO CALCULATE        //
  // MANUFACTURED RHS                   //
  //------------------------------------------

  AMP::Mesh::Mesh::shared_ptr  meshAdapter1 = manager->Subset( "LeftMesh" );
  AMP::Mesh::Mesh::shared_ptr  meshAdapter2 = manager->Subset( "RightMesh" );

  boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));
  boost::shared_ptr<AMP::Operator::OperatorParameters> columnParams (new AMP::Operator::OperatorParameters(tmp_db));
  boost::shared_ptr<AMP::Operator::ColumnOperator>     volumeIntegralColumnOperator (new AMP::Operator::ColumnOperator(columnParams));

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> sourceModel1;
  boost::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator1 = boost::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(AMP::Operator::OperatorBuilder::createOperator( meshAdapter1,"LeftVolumeIntegralOperator",input_db, sourceModel1));
  volumeIntegralColumnOperator->append(sourceOperator1);

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> sourceModel2;
  boost::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator2 = boost::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(AMP::Operator::OperatorBuilder::createOperator( meshAdapter2,"RightVolumeIntegralOperator",input_db, sourceModel2));
  volumeIntegralColumnOperator->append(sourceOperator2);

  AMP::LinearAlgebra::Variable::shared_ptr rhsVar = volumeIntegralColumnOperator-> getOutputVariable();
  AMP::LinearAlgebra::Vector::shared_ptr integratedRHSVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, rhsVar, true );
  integratedRHSVec->zero();

  volumeIntegralColumnOperator->apply(nullVec, manufacturedRHS, integratedRHSVec, 1.,0);  

  integratedRHSVec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );

#ifdef USE_SILO
  siloWriter->registerVector( integratedRHSVec, manager, AMP::Mesh::Vertex ,"Source");
#endif

  // modify the RHS to take into account boundary conditions
  //  for(int id = 0; id != boost::dynamic_pointer_cast<AMP::Operator::ColumnOperator>(nonlinearThermalColumnOperator)->getNumberOfOperators(); id++)
  for(int i=0; i<2; i++)
  {
    AMP::Operator::NonlinearBVPOperator::shared_ptr nonlinearThermalOperator ;
    nonlinearThermalOperator = (boost::dynamic_pointer_cast<AMP::Operator::ColumnOperator> ( nonlinearThermalColumnOperator ) )->getOperator(i);
    (boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator> (nonlinearThermalOperator))->modifyInitialSolutionVector(TemperatureVec);
    (boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator> (nonlinearThermalOperator))->modifyRHSvector(integratedRHSVec);
  }
  boost::shared_ptr<AMP::InputDatabase> emptyDb;
  boost::shared_ptr<AMP::Operator::CoupledOperatorParameters> thermalCoupledOpParams(new AMP::Operator::CoupledOperatorParameters(emptyDb));
  thermalCoupledOpParams->d_MapOperator = mapsColumn;
  thermalCoupledOpParams->d_BVPOperator = nonlinearThermalColumnOperator;
  boost::shared_ptr<AMP::Operator::Operator> nonlinearThermalCoupledOperator(new AMP::Operator::CoupledOperator(thermalCoupledOpParams));

  nonlinearThermalCoupledOperator->apply(integratedRHSVec, TemperatureVec, ResidualVec, 1.0, -1.0);
  ResidualVec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
  double initialResidualNorm  = ResidualVec->L2Norm();

  AMP::pout<<"Initial Residual Norm: "<<initialResidualNorm<<std::endl;

  boost::shared_ptr<AMP::Solver::PetscSNESSolver>  nonlinearThermalSolver;
  boost::shared_ptr<AMP::Solver::PetscKrylovSolver>  linearThermalSolver;
  boost::shared_ptr<AMP::Operator::Operator> linearThermalOperator  = boost::dynamic_pointer_cast<AMP::Operator::Operator> ( linearThermalColumnOperator ) ;
  createThermalSolvers(input_db, TemperatureVec ,  nonlinearThermalCoupledOperator, linearThermalOperator, nonlinearThermalSolver, linearThermalSolver);

  nonlinearThermalSolver->setZeroInitialGuess(false);

  nonlinearThermalSolver->solve(integratedRHSVec , TemperatureVec );

  solutionError->subtract(TemperatureVec, manufacturedSolution);

  std::cout << "Max of ||U-Uh|| : "<< solutionError->max() << " Min of ||U-Uh|| : "<< solutionError->min()<< std::endl;

  TemperatureVec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );

  double discretizationErrorNorm2; 
  double TotalNorm2 = 0;

  discretizationErrorNorm2 = 0;
  AMP::LinearAlgebra::VS_Mesh meshSelector1("meshSelector", meshAdapter1 );
  computeL2Norm( meshAdapter1 , globalComm, TemperatureVec->select(meshSelector1, "Temperature"), &discretizationErrorNorm2 );
  TotalNorm2 += discretizationErrorNorm2; 
  std::cout << "Discretized error norm ^2 for Mesh  1: "<< discretizationErrorNorm2 << std::endl;
  AMP::LinearAlgebra::VS_Mesh meshSelector2("meshSelector", meshAdapter2 );
  computeL2Norm( meshAdapter2 , globalComm, TemperatureVec->select(meshSelector2, "Temperature"), &discretizationErrorNorm2 );
  TotalNorm2 += discretizationErrorNorm2; 
  std::cout << "Discretized error norm ^2 for Mesh  2: "<< discretizationErrorNorm2 << std::endl;

  std::cout << "Discretized error norm for ||U-Uh|| : "<< sqrt(TotalNorm2) << std::endl;

  std::cout << "Max of U : "<< TemperatureVec->max() << " Min of U : "<< TemperatureVec->min()<< std::endl;


#ifdef USE_SILO
  siloWriter->writeFile( exeName, 1 );
#endif

}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;
  std::vector<std::string> exeNames;
  exeNames.push_back("testMeshRefinementDiffusion-1");

  for(unsigned int i=0; i < exeNames.size(); i++)
  {
    try
    {
      myTest(&ut, exeNames[i]);
    } catch (std::exception &err) {
      std::cout << "ERROR:While testing " << argv[0]<<err.what()<<std::endl;
      ut.failure("ERROR: While testing");
    } catch(...){
      std::cout   << "ERROR: While testing " << argv[0] << "An unknown exception was thrown." << std::endl;
      ut.failure("ERROR: While testing");
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;

  }

}
