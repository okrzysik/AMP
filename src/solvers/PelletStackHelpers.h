#ifndef included_AMP_PelletStackHelpers
#define included_AMP_PelletStackHelpers

#include "operators/OperatorBuilder.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/PelletStackOperator.h"
#include "operators/CoupledOperator.h"
#include "operators/map/NodeToNodeMap.h"
#include "operators/map/AsyncMapColumnOperator.h"
#include "operators/mechanics/MechanicsNonlinearFEOperator.h"

#include "discretization/simpleDOF_Manager.h"

#include "ampmesh/Mesh.h"
#include "vectors/VectorBuilder.h"
#include "matrices/MatrixBuilder.h"

#include "solvers/PetscKrylovSolver.h"
#include "solvers/ColumnSolver.h"
#include "solvers/PelletStackMechanicsSolver.h"
#include "solvers/trilinos/TrilinosMLSolver.h"

void helperCreateStackOperatorForPelletMechanics(AMP::Mesh::Mesh::shared_ptr manager,
    boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator> n2nmaps, 
    boost::shared_ptr<AMP::InputDatabase> global_input_db, 
    boost::shared_ptr<AMP::Operator::PelletStackOperator> & pelletStackOp);


void helperCreateColumnOperatorsForPelletMechanics(std::vector<unsigned int> localPelletIds, 
    std::vector<AMP::Mesh::Mesh::shared_ptr> localMeshes,
    boost::shared_ptr<AMP::InputDatabase> global_input_db,
    boost::shared_ptr<AMP::Operator::ColumnOperator> & nonlinearColumnOperator,
    boost::shared_ptr<AMP::Operator::ColumnOperator> & linearColumnOperator);


void helperCreateCoupledOperatorForPelletMechanics(boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator> n2nmaps, 
    boost::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator,
    boost::shared_ptr<AMP::Operator::CoupledOperator> & coupledOp);


void helperSetFrozenVectorForMapsForPelletMechanics(AMP::Mesh::Mesh::shared_ptr manager, 
    boost::shared_ptr<AMP::Operator::CoupledOperator> coupledOp);


void helperCreateAllOperatorsForPelletMechanics(AMP::Mesh::Mesh::shared_ptr manager,
    AMP::AMP_MPI globalComm, boost::shared_ptr<AMP::InputDatabase> global_input_db,
    boost::shared_ptr<AMP::Operator::CoupledOperator> & coupledOp,
    boost::shared_ptr<AMP::Operator::ColumnOperator> & linearColumnOperator, 
    boost::shared_ptr<AMP::Operator::PelletStackOperator> & pelletStackOp);


void helperCreateVectorsForPelletMechanics(AMP::Mesh::Mesh::shared_ptr manager,
    boost::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,  
    AMP::LinearAlgebra::Vector::shared_ptr & solVec,
    AMP::LinearAlgebra::Vector::shared_ptr & rhsVec,
    AMP::LinearAlgebra::Vector::shared_ptr & scaledRhsVec);


void helperBuildPointLoadRHSForPelletMechanics(boost::shared_ptr<AMP::InputDatabase> global_input_db, 
    boost::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec);


void helperApplyBoundaryCorrectionsForPelletMechanics(boost::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,
    AMP::LinearAlgebra::Vector::shared_ptr solVec, 
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec);


void helperCreateTemperatureVectorsForPelletMechanics(AMP::Mesh::Mesh::shared_ptr manager, 
    AMP::LinearAlgebra::Vector::shared_ptr & initialTemperatureVec,
    AMP::LinearAlgebra::Vector::shared_ptr & finalTemperatureVec);


void helperSetReferenceTemperatureForPelletMechanics(boost::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,
    AMP::LinearAlgebra::Vector::shared_ptr initialTemperatureVec);


void helperSetFinalTemperatureForPelletMechanics(boost::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,
    AMP::LinearAlgebra::Vector::shared_ptr finalTemperatureVec);


void helperBuildColumnSolverForPelletMechanics(boost::shared_ptr<AMP::Database> columnSolver_db,
    boost::shared_ptr<AMP::Operator::ColumnOperator> linearColumnOperator,
    boost::shared_ptr<AMP::Solver::ColumnSolver> & columnSolver);


void helperBuildStackSolverForPelletMechanics(boost::shared_ptr<AMP::Database> pelletStackSolver_db, 
    boost::shared_ptr<AMP::Operator::PelletStackOperator> pelletStackOp, 
    boost::shared_ptr<AMP::Operator::ColumnOperator> linearColumnOperator, 
    boost::shared_ptr<AMP::Solver::SolverStrategy> & pelletStackSolver);


void helperResetNonlinearOperatorForPelletMechanics(boost::shared_ptr<AMP::Operator::CoupledOperator> coupledOp);


#endif




