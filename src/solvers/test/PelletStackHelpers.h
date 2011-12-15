
#ifndef included_AMP_PelletStackHelpers
#define included_AMP_PelletStackHelpers

void helperCreatePelletStackOperator(AMP::Mesh::MeshManager::shared_ptr manager,
    boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator> n2nmaps, AMP::AMP_MPI globalComm,
    boost::shared_ptr<AMP::InputDatabase> global_input_db, 
    boost::shared_ptr<AMP::Operator::PelletStackOperator> & pelletStackOp) {
  boost::shared_ptr<AMP::Database> pelletStackOp_db = global_input_db->getDatabase("PelletStackOperator");
  boost::shared_ptr<AMP::Operator::PelletStackOperatorParameters> pelletStackOpParams(new 
      AMP::Operator::PelletStackOperatorParameters(pelletStackOp_db));
  pelletStackOpParams->d_pelletStackComm = globalComm;
  pelletStackOpParams->d_n2nMaps = n2nmaps;
  pelletStackOpParams->d_meshManager = manager;
  pelletStackOp.reset(new AMP::Operator::PelletStackOperator(pelletStackOpParams));
}

void helperCreateColumnOperators(std::vector<unsigned int> localPelletIds, 
    std::vector<AMP::Mesh::MeshManager::Adapter::shared_ptr> localMeshes,
    boost::shared_ptr<AMP::InputDatabase> global_input_db,
    boost::shared_ptr<AMP::Operator::ColumnOperator> & nonlinearColumnOperator,
    boost::shared_ptr<AMP::Operator::ColumnOperator> & linearColumnOperator) {
  boost::shared_ptr<AMP::Operator::OperatorParameters> emptyParams;
  nonlinearColumnOperator.reset(new AMP::Operator::ColumnOperator(emptyParams));
  linearColumnOperator.reset(new AMP::Operator::ColumnOperator(emptyParams));
  for(unsigned int id = 0; id < localPelletIds.size(); id++) {
    std::string prefix = "";
    if(localPelletIds[id] == 0)
    {
      prefix = "Bottom";
    }

    AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = localMeshes[id];

    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> mechModel;
    boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinOperator =
      boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
            prefix+"NonlinearMechanicsOperator", global_input_db, mechModel));
    nonlinearColumnOperator->append(nonlinOperator);

    boost::shared_ptr<AMP::Operator::LinearBVPOperator> linOperator =
      boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
            prefix+"LinearMechanicsOperator", global_input_db, mechModel));
    linearColumnOperator->append(linOperator);
  }//end for id
} 

void helperCreateCoupledOperator(boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator> n2nmaps, 
    boost::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator,
    boost::shared_ptr<AMP::Operator::CoupledOperator> & coupledOp) {
  boost::shared_ptr<AMP::InputDatabase> emptyDb;
  boost::shared_ptr<AMP::Operator::CoupledOperatorParameters> coupledOpParams(new
      AMP::Operator::CoupledOperatorParameters(emptyDb));
  coupledOpParams->d_MapOperator = n2nmaps;
  coupledOpParams->d_BVPOperator = nonlinearColumnOperator;
  coupledOp.reset(new AMP::Operator::CoupledOperator(coupledOpParams));
}

void helperSetFrozenVectorForMaps(AMP::Mesh::MeshManager::shared_ptr manager, AMP::AMP_MPI globalComm,
    boost::shared_ptr<AMP::Operator::CoupledOperator> coupledOp) {
  boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator> n2nmaps = 
    boost::dynamic_pointer_cast<AMP::Operator::AsyncMapColumnOperator>(coupledOp->getOperator(1)); 
  boost::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::ColumnOperator>(coupledOp->getOperator(2));
  AMP::LinearAlgebra::Variable::shared_ptr dispVar = nonlinearColumnOperator->getOutputVariable();
  AMP::LinearAlgebra::Vector::shared_ptr dirichletValues = AMP::LinearAlgebra::CommCollectVector::view (
      manager->createVector ( dispVar ) , globalComm );
  n2nmaps->setVector(dirichletValues);
  for(int id = 0; id < nonlinearColumnOperator->getNumberOfOperators(); id++) {
    boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletOp =
      boost::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
          boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            nonlinearColumnOperator->getOperator(id))->getBoundaryOperator());
    dirichletOp->setDirichletValues(dirichletValues);
  }//end for id
}

void helperCreateVectors(AMP::Mesh::MeshManager::shared_ptr manager, AMP::Operator::Operator::shared_ptr nonlinearColumnOperator,
    AMP::AMP_MPI globalComm, AMP::LinearAlgebra::Vector::shared_ptr & solVec, AMP::LinearAlgebra::Vector::shared_ptr & rhsVec,
    AMP::LinearAlgebra::Vector::shared_ptr & scaledRhsVec) {
  AMP::LinearAlgebra::Variable::shared_ptr dispVar = nonlinearColumnOperator->getOutputVariable();
  solVec = AMP::LinearAlgebra::CommCollectVector::view (
      manager->createVector ( dispVar ) , globalComm );
  rhsVec = AMP::LinearAlgebra::CommCollectVector::view (
      manager->createVector ( dispVar ) , globalComm );
  scaledRhsVec = AMP::LinearAlgebra::CommCollectVector::view (
      manager->createVector ( dispVar ) , globalComm );
}

void helperBuildPointLoadRHS(boost::shared_ptr<AMP::InputDatabase> global_input_db, 
    boost::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator,
    AMP::LinearAlgebra::Vector::shared_ptr & rhsVec) {
  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  rhsVec->zero();
  for(int id = 0; id < nonlinearColumnOperator->getNumberOfOperators(); id++) {
    AMP::Operator::Operator::shared_ptr currOp = nonlinearColumnOperator->getOperator(id);
    AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = currOp->getMeshAdapter();
    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
    boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> loadOp = 
      boost::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
          AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
            "PointLoad", global_input_db, dummyModel));
    loadOp->setVariable(currOp->getOutputVariable());
    loadOp->apply(nullVec, nullVec, rhsVec, 1.0, 0.0);
  }//end for id
}

void helperApplyBoundaryCorrections(boost::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator,
    AMP::LinearAlgebra::Vector::shared_ptr solVec, 
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec) {
  for(int id = 0; id < nonlinearColumnOperator->getNumberOfOperators(); id++) {
    boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinOperator =
      boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
          nonlinearColumnOperator->getOperator(id));
    nonlinOperator->modifyInitialSolutionVector(solVec);
    nonlinOperator->modifyRHSvector(rhsVec);
  }//end for id
}

#endif


