#ifndef included_AMP_PelletStackHelpers
#define included_AMP_PelletStackHelpers

#include "AMP/ampmesh/Mesh.h"
#include "AMP/operators/CoupledOperator.h"
#include "AMP/operators/libmesh/PelletStackOperator.h"
#include "AMP/operators/map/AsyncMapColumnOperator.h"
#include "AMP/solvers/ColumnSolver.h"
#include "AMP/solvers/libmesh/PelletStackMechanicsSolver.h"
#include "AMP/vectors/Vector.h"


void helperCreateStackOperatorForPelletMechanics(
    AMP::Mesh::Mesh::shared_ptr manager,
    std::shared_ptr<AMP::Operator::AsyncMapColumnOperator> n2nmaps,
    std::shared_ptr<AMP::Database> global_input_db,
    std::shared_ptr<AMP::Operator::PelletStackOperator> &pelletStackOp );


void helperCreateColumnOperatorsForPelletMechanics(
    std::vector<unsigned int> localPelletIds,
    std::vector<AMP::Mesh::Mesh::shared_ptr> localMeshes,
    std::shared_ptr<AMP::Database> global_input_db,
    std::shared_ptr<AMP::Operator::ColumnOperator> &nonlinearColumnOperator,
    std::shared_ptr<AMP::Operator::ColumnOperator> &linearColumnOperator );


void helperCreateCoupledOperatorForPelletMechanics(
    std::shared_ptr<AMP::Operator::AsyncMapColumnOperator> n2nmaps,
    std::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator,
    std::shared_ptr<AMP::Operator::CoupledOperator> &coupledOp );


void helperSetFrozenVectorForMapsForPelletMechanics(
    AMP::Mesh::Mesh::shared_ptr manager,
    std::shared_ptr<AMP::Operator::CoupledOperator> coupledOp );


void helperCreateAllOperatorsForPelletMechanics(
    AMP::Mesh::Mesh::shared_ptr manager,
    AMP::AMP_MPI globalComm,
    std::shared_ptr<AMP::Database> global_input_db,
    std::shared_ptr<AMP::Operator::CoupledOperator> &coupledOp,
    std::shared_ptr<AMP::Operator::ColumnOperator> &linearColumnOperator,
    std::shared_ptr<AMP::Operator::PelletStackOperator> &pelletStackOp );


void helperCreateVectorsForPelletMechanics(
    AMP::Mesh::Mesh::shared_ptr manager,
    std::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,
    AMP::LinearAlgebra::Vector::shared_ptr &solVec,
    AMP::LinearAlgebra::Vector::shared_ptr &rhsVec,
    AMP::LinearAlgebra::Vector::shared_ptr &scaledRhsVec );


void helperBuildPointLoadRHSForPelletMechanics(
    std::shared_ptr<AMP::Database> global_input_db,
    std::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec );


void helperApplyBoundaryCorrectionsForPelletMechanics(
    std::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,
    AMP::LinearAlgebra::Vector::shared_ptr solVec,
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec );


void helperCreateTemperatureVectorsForPelletMechanics(
    AMP::Mesh::Mesh::shared_ptr manager,
    AMP::LinearAlgebra::Vector::shared_ptr &initialTemperatureVec,
    AMP::LinearAlgebra::Vector::shared_ptr &finalTemperatureVec );


void helperSetReferenceTemperatureForPelletMechanics(
    std::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,
    AMP::LinearAlgebra::Vector::shared_ptr initialTemperatureVec );


void helperSetFinalTemperatureForPelletMechanics(
    std::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,
    AMP::LinearAlgebra::Vector::shared_ptr finalTemperatureVec );


void helperBuildColumnSolverForPelletMechanics(
    std::shared_ptr<AMP::Database> columnSolver_db,
    std::shared_ptr<AMP::Operator::ColumnOperator> linearColumnOperator,
    std::shared_ptr<AMP::Solver::ColumnSolver> &columnSolver );


void helperBuildStackSolverForPelletMechanics(
    std::shared_ptr<AMP::Database> pelletStackSolver_db,
    std::shared_ptr<AMP::Operator::PelletStackOperator> pelletStackOp,
    std::shared_ptr<AMP::Operator::ColumnOperator> linearColumnOperator,
    std::shared_ptr<AMP::Solver::SolverStrategy> &pelletStackSolver );


void helperResetNonlinearOperatorForPelletMechanics(
    std::shared_ptr<AMP::Operator::CoupledOperator> coupledOp );


#endif
