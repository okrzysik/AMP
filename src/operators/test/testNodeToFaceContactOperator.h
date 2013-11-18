void rotateMesh(AMP::Mesh::Mesh::shared_ptr mesh);


void rotateMesh(AMP::Mesh::Mesh::shared_ptr mesh) {
  AMP::Discretization::DOFManager::shared_ptr dofManager = AMP::Discretization::simpleDOFManager::create(mesh, AMP::Mesh::Vertex, 0, 3, false);
  AMP::LinearAlgebra::Variable::shared_ptr dispVar(new AMP::LinearAlgebra::Variable("disp"));
  AMP::LinearAlgebra::Vector::shared_ptr dispVec = AMP::LinearAlgebra::createVector(dofManager, dispVar, false);
  dispVec->zero();
  AMP::Mesh::MeshIterator meshIterator = mesh->getIterator(AMP::Mesh::Vertex);
  AMP::Mesh::MeshIterator meshIterator_begin = meshIterator.begin();
  AMP::Mesh::MeshIterator meshIterator_end = meshIterator.end();
  for (meshIterator = meshIterator_begin; meshIterator != meshIterator_end; ++meshIterator) {
    std::vector<double> oldVertexCoord = meshIterator->coord();
    std::vector<size_t> dofIndices;
    dofManager->getDOFs(meshIterator->globalID(), dofIndices);
    AMP_ASSERT(dofIndices.size() == 3);
    std::vector<double> newVertexCoord(oldVertexCoord);
    // 15 degrees around the z axis
//    rotate_points(2, M_PI / 12.0, 1, &(newVertexCoord[0]));
    rotate_points(2, M_PI / 2.0, 1, &(newVertexCoord[0]));
    std::vector<double> vertexDisp(3, 0.0);
    make_vector_from_two_points(&(oldVertexCoord[0]), &(newVertexCoord[0]), &(vertexDisp[0]));  
    dispVec->setLocalValuesByGlobalID(3, &(dofIndices[0]), &(vertexDisp[0]));
  } // end for
  mesh->displaceMesh(dispVec);
}

struct dummyClass {
  bool operator() (std::pair<AMP::Mesh::MeshElementID, size_t> const & lhs, 
                   std::pair<AMP::Mesh::MeshElementID, size_t> const & rhs) {
    return lhs.first < rhs.first;
  }
};


void computeStressTensor(AMP::Mesh::Mesh::shared_ptr mesh, AMP::LinearAlgebra::Vector::shared_ptr displacementField, 
    AMP::LinearAlgebra::Vector::shared_ptr sigma_xx, AMP::LinearAlgebra::Vector::shared_ptr sigma_yy, AMP::LinearAlgebra::Vector::shared_ptr sigma_zz, 
    AMP::LinearAlgebra::Vector::shared_ptr sigma_yz, AMP::LinearAlgebra::Vector::shared_ptr sigma_xz, AMP::LinearAlgebra::Vector::shared_ptr sigma_xy, 
    AMP::LinearAlgebra::Vector::shared_ptr sigma_eff, double YoungsModulus, double PoissonsRatio, 
    double referenceTemperature = 273.0, double thermalExpansionCoefficient = 2.0e-6, 
    AMP::LinearAlgebra::Vector::shared_ptr temperatureField = AMP::LinearAlgebra::Vector::shared_ptr()) { 

  sigma_eff->setToScalar(-1.0);

  double constitutiveMatrix[36];
  compute_constitutive_matrix(YoungsModulus, PoissonsRatio, constitutiveMatrix);

  double verticesLocalCoordinates[24] = {
      -1.0, -1.0, -1.0,
      +1.0, -1.0, -1.0,
      +1.0, +1.0, -1.0,
      -1.0, +1.0, -1.0,
      -1.0, -1.0, +1.0,
      +1.0, -1.0, +1.0,
      +1.0, +1.0, +1.0,
      -1.0, +1.0, +1.0
  };

  std::set<std::pair<AMP::Mesh::MeshElementID, size_t>, dummyClass> verticesGlobalIDsAndCount;

  size_t numLocalVertices = mesh->numLocalElements(AMP::Mesh::Vertex);
  size_t countVertices = 0;
  std::vector<size_t> verticesInHowManyVolumeElements(numLocalVertices, 1);
  std::vector<size_t> verticesDOFIndex(numLocalVertices);
  AMP::Mesh::MeshIterator meshIterator = mesh->getIterator(AMP::Mesh::Volume);
  AMP::Mesh::MeshIterator meshIterator_begin = meshIterator.begin();
  AMP::Mesh::MeshIterator meshIterator_end = meshIterator.end();
  for (meshIterator = meshIterator_begin; meshIterator != meshIterator_end; ++meshIterator) {
    std::vector<AMP::Mesh::MeshElement> volumeVertices = meshIterator->getElements(AMP::Mesh::Vertex);
    AMP_ASSERT(volumeVertices.size() == 8);
    std::vector<AMP::Mesh::MeshElementID> volumeVerticesGlobalIDs(8);
    double volumeVerticesCoordinates[24];
    for (size_t v = 0; v < 8; ++v) {
      std::vector<double> vertexCoord = volumeVertices[v].coord();
      std::copy(vertexCoord.begin(), vertexCoord.end(), &(volumeVerticesCoordinates[3*v]));
      volumeVerticesGlobalIDs[v] = volumeVertices[v].globalID();
    } // end for v
    hex8_element_t volumeElement(volumeVerticesCoordinates);

    std::vector<size_t> displacementIndices;
    displacementField->getDOFManager()->getDOFs(volumeVerticesGlobalIDs, displacementIndices);
    AMP_ASSERT(displacementIndices.size() == 24);
    double displacementValues[24];
    displacementField->getValuesByGlobalID(24, &(displacementIndices[0]), &(displacementValues[0]));

    double temperatureValues[8];
    if (temperatureField.get() != NULL) {
      std::vector<size_t> temperatureIndices;
      temperatureField->getDOFManager()->getDOFs(volumeVerticesGlobalIDs, temperatureIndices);
      AMP_ASSERT(temperatureIndices.size() == 8);
      temperatureField->getValuesByGlobalID(8, &(temperatureIndices[0]), &(temperatureValues[0]));
    } // end if

    double strainTensor[6];
    double stressTensor[6];
    for (size_t v = 0; v < 8; ++v) {
      volumeElement.compute_strain_tensor(&(verticesLocalCoordinates[3*v]), displacementValues, strainTensor);
      if (temperatureField.get() != NULL) {
        for (size_t i = 0; i < 3; ++i) {
          strainTensor[i] -= (thermalExpansionCoefficient * (temperatureValues[v] - referenceTemperature));
        } // end for i
      } // end if
      compute_stress_tensor(constitutiveMatrix, strainTensor, stressTensor);
      double vonMisesStress = compute_von_mises_stress(stressTensor);

      std::vector<size_t> indices;
      sigma_xx->getDOFManager()->getDOFs(volumeVerticesGlobalIDs[v], indices);
      AMP_ASSERT(indices.size() == 1);
      // good job buddy
      std::pair<std::set<std::pair<AMP::Mesh::MeshElementID, size_t>, dummyClass>::iterator, bool> dummy = verticesGlobalIDsAndCount.insert(std::pair<AMP::Mesh::MeshElementID, size_t>(volumeVerticesGlobalIDs[v], countVertices));
      if (dummy.second) {
        verticesDOFIndex[countVertices] = indices[0];
        ++countVertices;
        sigma_xx->setLocalValueByGlobalID(indices[0], stressTensor[0]);
        sigma_yy->setLocalValueByGlobalID(indices[0], stressTensor[1]);
        sigma_zz->setLocalValueByGlobalID(indices[0], stressTensor[2]);
        sigma_yz->setLocalValueByGlobalID(indices[0], stressTensor[3]);
        sigma_xz->setLocalValueByGlobalID(indices[0], stressTensor[4]);
        sigma_xy->setLocalValueByGlobalID(indices[0], stressTensor[5]);
        sigma_eff->setLocalValueByGlobalID(indices[0], vonMisesStress);
      } else {
        // sigh...
        ++verticesInHowManyVolumeElements[dummy.first->second];
        sigma_xx->addLocalValueByGlobalID(indices[0], stressTensor[0]);
        sigma_yy->addLocalValueByGlobalID(indices[0], stressTensor[1]);
        sigma_zz->addLocalValueByGlobalID(indices[0], stressTensor[2]);
        sigma_yz->addLocalValueByGlobalID(indices[0], stressTensor[3]);
        sigma_xz->addLocalValueByGlobalID(indices[0], stressTensor[4]);
        sigma_xy->addLocalValueByGlobalID(indices[0], stressTensor[5]);
        sigma_eff->addLocalValueByGlobalID(indices[0], vonMisesStress);

/*        double deltaStressTensor[6];
        deltaStressTensor[0] = sigma_xx->getLocalValueByGlobalID(indices[0]) - stressTensor[0];
        deltaStressTensor[1] = sigma_yy->getLocalValueByGlobalID(indices[0]) - stressTensor[1];
        deltaStressTensor[2] = sigma_zz->getLocalValueByGlobalID(indices[0]) - stressTensor[2];
        deltaStressTensor[3] = sigma_yz->getLocalValueByGlobalID(indices[0]) - stressTensor[3];
        deltaStressTensor[4] = sigma_xz->getLocalValueByGlobalID(indices[0]) - stressTensor[4];
        deltaStressTensor[5] = sigma_xy->getLocalValueByGlobalID(indices[0]) - stressTensor[5];
        double deltaVonMisesStress = sigma_eff->getLocalValueByGlobalID(indices[0]) - vonMisesStress;
        std::cout<<"xx  "<<deltaStressTensor[0]<<"  "<<100.0*deltaStressTensor[0]/std::abs(stressTensor[0])<<" %  |  ";
        std::cout<<"yy  "<<deltaStressTensor[1]<<"  "<<100.0*deltaStressTensor[1]/std::abs(stressTensor[1])<<" %  |  ";
        std::cout<<"zz  "<<deltaStressTensor[2]<<"  "<<100.0*deltaStressTensor[2]/std::abs(stressTensor[2])<<" %  |  ";
        std::cout<<"yz  "<<deltaStressTensor[3]<<"  "<<100.0*deltaStressTensor[3]/std::abs(stressTensor[3])<<" %  |  ";
        std::cout<<"xz  "<<deltaStressTensor[4]<<"  "<<100.0*deltaStressTensor[4]/std::abs(stressTensor[4])<<" %  |  ";
        std::cout<<"xy  "<<deltaStressTensor[5]<<"  "<<100.0*deltaStressTensor[5]/std::abs(stressTensor[5])<<" %  |  ";
        std::cout<<"eff  "<<deltaVonMisesStress<<"  "<<100.0*deltaVonMisesStress/std::abs(vonMisesStress)<<" % \n";
*/
      } // end if
      
    } // end for v
  } // end for
  AMP_ASSERT(verticesGlobalIDsAndCount.size() == numLocalVertices);
  AMP_ASSERT(countVertices == numLocalVertices);
  AMP_ASSERT(find(verticesInHowManyVolumeElements.begin(), verticesInHowManyVolumeElements.end(), 0) == verticesInHowManyVolumeElements.end());

  for (size_t v = 0; v < numLocalVertices; ++v) {
    if (verticesInHowManyVolumeElements[v] > 1) {
//      AMP_ASSERT(verticesInHowManyVolumeElements[v] < 9);
      sigma_xx->setLocalValueByGlobalID(verticesDOFIndex[v], sigma_xx->getLocalValueByGlobalID(verticesDOFIndex[v]) / static_cast<double>(verticesInHowManyVolumeElements[v]));
      sigma_yy->setLocalValueByGlobalID(verticesDOFIndex[v], sigma_yy->getLocalValueByGlobalID(verticesDOFIndex[v]) / static_cast<double>(verticesInHowManyVolumeElements[v]));
      sigma_zz->setLocalValueByGlobalID(verticesDOFIndex[v], sigma_zz->getLocalValueByGlobalID(verticesDOFIndex[v]) / static_cast<double>(verticesInHowManyVolumeElements[v]));
      sigma_yz->setLocalValueByGlobalID(verticesDOFIndex[v], sigma_yz->getLocalValueByGlobalID(verticesDOFIndex[v]) / static_cast<double>(verticesInHowManyVolumeElements[v]));
      sigma_xz->setLocalValueByGlobalID(verticesDOFIndex[v], sigma_xz->getLocalValueByGlobalID(verticesDOFIndex[v]) / static_cast<double>(verticesInHowManyVolumeElements[v]));
      sigma_xy->setLocalValueByGlobalID(verticesDOFIndex[v], sigma_xy->getLocalValueByGlobalID(verticesDOFIndex[v]) / static_cast<double>(verticesInHowManyVolumeElements[v]));
      sigma_eff->setLocalValueByGlobalID(verticesDOFIndex[v], sigma_eff->getLocalValueByGlobalID(verticesDOFIndex[v]) / static_cast<double>(verticesInHowManyVolumeElements[v]));
    } //end if
  } // end for v
  
} 
 

void computeStressTensor(AMP::Mesh::Mesh::shared_ptr mesh, AMP::LinearAlgebra::Vector::shared_ptr displacementField, 
    AMP::LinearAlgebra::Vector::shared_ptr sigmaXX, AMP::LinearAlgebra::Vector::shared_ptr sigmaYY, AMP::LinearAlgebra::Vector::shared_ptr sigmaZZ, 
    AMP::LinearAlgebra::Vector::shared_ptr sigmaYZ, AMP::LinearAlgebra::Vector::shared_ptr sigmaXZ, AMP::LinearAlgebra::Vector::shared_ptr sigmaXY, 
    AMP::LinearAlgebra::Vector::shared_ptr sigmaEff, boost::shared_ptr<AMP::Operator::MechanicsMaterialModel> mechanicsMaterialModel,
    double referenceTemperature = 273.0, double thermalExpansionCoefficient = 2.0e-6, 
    AMP::LinearAlgebra::Vector::shared_ptr temperatureField = AMP::LinearAlgebra::Vector::shared_ptr()) { 

  AMP::LinearAlgebra::VS_Mesh vectorSelector(mesh);
  AMP::LinearAlgebra::Vector::shared_ptr subsetDisplacementField = displacementField->select(vectorSelector, (displacementField->getVariable())->getName());
  AMP::LinearAlgebra::Vector::shared_ptr subsetSigmaXX = sigmaXX->select(vectorSelector, (sigmaXX->getVariable())->getName());
  AMP::LinearAlgebra::Vector::shared_ptr subsetSigmaYY = sigmaYY->select(vectorSelector, (sigmaYY->getVariable())->getName());
  AMP::LinearAlgebra::Vector::shared_ptr subsetSigmaZZ = sigmaZZ->select(vectorSelector, (sigmaZZ->getVariable())->getName());
  AMP::LinearAlgebra::Vector::shared_ptr subsetSigmaYZ = sigmaYZ->select(vectorSelector, (sigmaYZ->getVariable())->getName());
  AMP::LinearAlgebra::Vector::shared_ptr subsetSigmaXZ = sigmaXZ->select(vectorSelector, (sigmaXZ->getVariable())->getName());
  AMP::LinearAlgebra::Vector::shared_ptr subsetSigmaXY = sigmaXY->select(vectorSelector, (sigmaXY->getVariable())->getName());
  AMP::LinearAlgebra::Vector::shared_ptr subsetSigmaEff = sigmaEff->select(vectorSelector, (sigmaEff->getVariable())->getName());
  AMP::LinearAlgebra::Vector::shared_ptr subsetTemperatureField = temperatureField->select(vectorSelector, (temperatureField->getVariable())->getName());

  double youngsModulus = boost::dynamic_pointer_cast<AMP::Operator::IsotropicElasticModel>(mechanicsMaterialModel)->getYoungsModulus();
  double poissonsRatio = boost::dynamic_pointer_cast<AMP::Operator::IsotropicElasticModel>(mechanicsMaterialModel)->getPoissonsRatio();

  computeStressTensor(mesh, subsetDisplacementField,
      subsetSigmaXX, subsetSigmaYY, subsetSigmaZZ, subsetSigmaYZ, subsetSigmaXZ, subsetSigmaXY,
      subsetSigmaEff, youngsModulus, poissonsRatio,
      referenceTemperature, thermalExpansionCoefficient, subsetTemperatureField);
}



void drawVerticesOnBoundaryID(AMP::Mesh::Mesh::shared_ptr meshAdapter, int boundaryID, std::ostream &os, double const * point_of_view, const std::string & option = "") {
  AMP::Mesh::MeshIterator boundaryIterator = meshAdapter->getBoundaryIDIterator(AMP::Mesh::Vertex, boundaryID);
  AMP::Mesh::MeshIterator boundaryIterator_begin = boundaryIterator.begin(), 
      boundaryIterator_end = boundaryIterator.end();
  std::vector<double> vertexCoordinates;

  os<<std::setprecision(6)<<std::fixed;

  for (boundaryIterator = boundaryIterator_begin; boundaryIterator != boundaryIterator_end; ++boundaryIterator) {
    vertexCoordinates = boundaryIterator->coord();
    AMP_ASSERT( vertexCoordinates.size() == 3 );
    draw_point(&(vertexCoordinates[0]), option, os);
  } // end for
}

void drawFacesOnBoundaryID(AMP::Mesh::Mesh::shared_ptr meshAdapter, int boundaryID, std::ostream &os, double const * point_of_view, const std::string & option = "") {
  AMP::Mesh::MeshIterator boundaryIterator = meshAdapter->getBoundaryIDIterator(AMP::Mesh::Face, boundaryID);
  AMP::Mesh::MeshIterator boundaryIterator_begin = boundaryIterator.begin(), 
      boundaryIterator_end = boundaryIterator.end();
  std::vector<AMP::Mesh::MeshElement> faceVertices;
  std::vector<double> faceVertexCoordinates;
  double faceData[12];
  double const * faceDataPtr[4] = { faceData, faceData+3, faceData+6, faceData+9 };

  os<<std::setprecision(6)<<std::fixed;

  for (boundaryIterator = boundaryIterator_begin; boundaryIterator != boundaryIterator_end; ++boundaryIterator) {
    faceVertices = boundaryIterator->getElements(AMP::Mesh::Vertex);
    AMP_ASSERT( faceVertices.size() == 4 );
    for (size_t i = 0; i < 4; ++i) {
      faceVertexCoordinates = faceVertices[i].coord();
      AMP_ASSERT( faceVertexCoordinates.size() == 3 );
      std::copy(faceVertexCoordinates.begin(), faceVertexCoordinates.end(), faceData+3*i);
    } // end for i
    triangle_t t(faceDataPtr[0], faceDataPtr[1], faceDataPtr[2]);

//    if (compute_scalar_product(point_of_view, t.get_normal()) > 0.0) {
    if (true) {
      os<<"\\draw["<<option<<"]\n";
      write_face(faceDataPtr, os);
    } // end if
  } // end for
}

void myPCG(AMP::LinearAlgebra::Vector::shared_ptr rhs, AMP::LinearAlgebra::Vector::shared_ptr sol, 
    AMP::Operator::Operator::shared_ptr op, boost::shared_ptr<AMP::Solver::SolverStrategy> pre,
    size_t maxIters, double relTol, double absTol, bool verbose = false, std::ostream &os = std::cout) {
  AMP::LinearAlgebra::Vector::shared_ptr res = sol->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr dir = sol->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr ext = sol->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr oldSol = sol->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr oldRes = sol->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr oldDir = sol->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr matVec = sol->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr nullVec;

  op->apply(nullVec, sol, matVec, 1.0, 0.0);
  oldRes->subtract(rhs, matVec);
  pre->solve(oldRes, ext);
  oldDir->copyVector(ext);
  oldSol->copyVector(sol);
  double initialResNorm = oldRes->L2Norm();
  AMP::AMP_MPI comm = sol->getComm();
  int rank = comm.getRank();
  verbose = verbose && !rank;
  double tol = absTol + relTol * initialResNorm;
  if (verbose) { os<<std::setprecision(15)<<"  iter=0  itialResNorm="<<initialResNorm<<"\n"; }
  for (size_t iter = 0; iter < maxIters; ++iter) {
    if (verbose) { os<<"  iter="<<iter+1<<"  "; }
    op->apply(nullVec, oldDir, matVec, 1.0, 0.0);
    double extDOToldRes = ext->dot(oldRes);
    double oldDirDOTmatVec = oldDir->dot(matVec);
    double alpha = extDOToldRes / oldDirDOTmatVec;
    if (verbose) { os<<"alpha="<<alpha<<"  "; }
    if (verbose) { os<<"oldDirDOTmatVec="<<oldDirDOTmatVec<<"  "; }
    sol->axpy(alpha, oldDir, oldSol);
    res->axpy(-alpha, matVec, oldRes);
    double resNorm = res->L2Norm();
    if (verbose) { os<<"resNorm="<<resNorm<<"  "; }
    if (resNorm < tol) { os<<"\n"; break; }
    pre->solve(res, ext);
    double extDOTres = ext->dot(res);
    double beta = extDOTres / extDOToldRes;
    if (verbose) { os<<"beta="<<beta<<"\n"; }
    dir->axpy(beta, oldDir, ext);
    oldSol->copyVector(sol);
    oldRes->copyVector(res);
    oldDir->copyVector(dir);
  } // end for
}
