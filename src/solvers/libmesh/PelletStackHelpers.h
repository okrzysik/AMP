#ifndef included_AMP_PelletStackHelpers
#define included_AMP_PelletStackHelpers

#include "ampmesh/Mesh.h"
#include "operators/CoupledOperator.h"
#include "operators/libmesh/PelletStackOperator.h"
#include "operators/map/AsyncMapColumnOperator.h"
#include "solvers/ColumnSolver.h"
#include "solvers/libmesh/PelletStackMechanicsSolver.h"
#include "solvers/trilinos/TrilinosMLSolver.h"
#include "vectors/Vector.h"


void helperCreateStackOperatorForPelletMechanics(
    AMP::Mesh::Mesh::shared_ptr manager,
    AMP::shared_ptr<AMP::Operator::AsyncMapColumnOperator>
        n2nmaps,
    AMP::shared_ptr<AMP::Database>
        global_input_db,
    AMP::shared_ptr<AMP::Operator::PelletStackOperator> &pelletStackOp );


void helperCreateColumnOperatorsForPelletMechanics(
    std::vector<unsigned int> localPelletIds,
    std::vector<AMP::Mesh::Mesh::shared_ptr>
        localMeshes,
    AMP::shared_ptr<AMP::Database>
        global_input_db,
    AMP::shared_ptr<AMP::Operator::ColumnOperator> &nonlinearColumnOperator,
    AMP::shared_ptr<AMP::Operator::ColumnOperator> &linearColumnOperator );


void helperCreateCoupledOperatorForPelletMechanics(
    AMP::shared_ptr<AMP::Operator::AsyncMapColumnOperator> n2nmaps,
    AMP::shared_ptr<AMP::Operator::ColumnOperator>
        nonlinearColumnOperator,
    AMP::shared_ptr<AMP::Operator::CoupledOperator> &coupledOp );


void helperSetFrozenVectorForMapsForPelletMechanics( AMP::Mesh::Mesh::shared_ptr manager,
                                                     AMP::shared_ptr<AMP::Operator::CoupledOperator>
                                                         coupledOp );


void helperCreateAllOperatorsForPelletMechanics(
    AMP::Mesh::Mesh::shared_ptr manager,
    AMP::AMP_MPI globalComm,
    AMP::shared_ptr<AMP::Database>
        global_input_db,
    AMP::shared_ptr<AMP::Operator::CoupledOperator> &coupledOp,
    AMP::shared_ptr<AMP::Operator::ColumnOperator> &linearColumnOperator,
    AMP::shared_ptr<AMP::Operator::PelletStackOperator> &pelletStackOp );


void helperCreateVectorsForPelletMechanics( AMP::Mesh::Mesh::shared_ptr manager,
                                            AMP::shared_ptr<AMP::Operator::CoupledOperator>
                                                coupledOp,
                                            AMP::LinearAlgebra::Vector::shared_ptr &solVec,
                                            AMP::LinearAlgebra::Vector::shared_ptr &rhsVec,
                                            AMP::LinearAlgebra::Vector::shared_ptr &scaledRhsVec );


void helperBuildPointLoadRHSForPelletMechanics( AMP::shared_ptr<AMP::Database> global_input_db,
                                                AMP::shared_ptr<AMP::Operator::CoupledOperator>
                                                    coupledOp,
                                                AMP::LinearAlgebra::Vector::shared_ptr rhsVec );


void helperApplyBoundaryCorrectionsForPelletMechanics(
    AMP::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,
    AMP::LinearAlgebra::Vector::shared_ptr solVec,
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec );


void helperCreateTemperatureVectorsForPelletMechanics(
    AMP::Mesh::Mesh::shared_ptr manager,
    AMP::LinearAlgebra::Vector::shared_ptr &initialTemperatureVec,
    AMP::LinearAlgebra::Vector::shared_ptr &finalTemperatureVec );


void helperSetReferenceTemperatureForPelletMechanics(
    AMP::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,
    AMP::LinearAlgebra::Vector::shared_ptr initialTemperatureVec );


void helperSetFinalTemperatureForPelletMechanics(
    AMP::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,
    AMP::LinearAlgebra::Vector::shared_ptr finalTemperatureVec );


void helperBuildColumnSolverForPelletMechanics(
    AMP::shared_ptr<AMP::Database> columnSolver_db,
    AMP::shared_ptr<AMP::Operator::ColumnOperator>
        linearColumnOperator,
    AMP::shared_ptr<AMP::Solver::ColumnSolver> &columnSolver );


void helperBuildStackSolverForPelletMechanics(
    AMP::shared_ptr<AMP::Database> pelletStackSolver_db,
    AMP::shared_ptr<AMP::Operator::PelletStackOperator>
        pelletStackOp,
    AMP::shared_ptr<AMP::Operator::ColumnOperator>
        linearColumnOperator,
    AMP::shared_ptr<AMP::Solver::SolverStrategy> &pelletStackSolver );


void helperResetNonlinearOperatorForPelletMechanics(
    AMP::shared_ptr<AMP::Operator::CoupledOperator> coupledOp );


#endif
