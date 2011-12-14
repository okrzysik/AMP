
#ifndef included_AMP_PelletStackHelpers
#define included_AMP_PelletStackHelpers

void helperCreateVectors(AMP::Mesh::MeshManager::shared_ptr manager, AMP::Operator::Operator::shared_ptr nonlinearColumnOperator,
    AMP::AMP_MPI globalComm, AMP::LinearAlgebra::Vector::shared_ptr & solVec, AMP::LinearAlgebra::Vector::shared_ptr & rhsVec,
    AMP::LinearAlgebra::Vector::shared_ptr & scaledRhsVec, AMP::LinearAlgebra::Vector::shared_ptr & resVec) {
  AMP::LinearAlgebra::Variable::shared_ptr dispVar = nonlinearColumnOperator->getOutputVariable();
  solVec = AMP::LinearAlgebra::CommCollectVector::view (
      manager->createVector ( dispVar ) , globalComm );
  rhsVec = AMP::LinearAlgebra::CommCollectVector::view (
      manager->createVector ( dispVar ) , globalComm );
  scaledRhsVec = AMP::LinearAlgebra::CommCollectVector::view (
      manager->createVector ( dispVar ) , globalComm );
  resVec = AMP::LinearAlgebra::CommCollectVector::view (
      manager->createVector ( dispVar ) , globalComm );
}


#endif


