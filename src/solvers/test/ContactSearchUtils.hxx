
#include "fe_interface.h"
#include "cell_hex8.h"

bool myContainsPoint(::Elem* e, const ::Point& p, double tol) {
  ::FEType fe_type(e->default_order());

  const ::Point mapped_point = ::FEInterface::inverse_map(e->dim(),
      fe_type, e, p, tol, false);

  return ::FEInterface::on_reference_element(mapped_point, e->type(), tol);
}

void computeSlave2MasterNodes(const double precision, const unsigned int slaveId, const unsigned int masterId, 
    AMP::Mesh::MeshManager::Adapter::shared_ptr slaveMeshAdapter, 
    AMP::Mesh::MeshManager::Adapter::shared_ptr masterMeshAdapter,
    std::vector<size_t> const & slaveNodes, std::vector<size_t> const & slave2MasterElem,
    std::vector<std::vector<size_t> > & slave2MasterNodes,
    std::vector<std::vector<double> > & slave2MasterFactors) {

  slave2MasterNodes.resize(slaveNodes.size());
  slave2MasterFactors.resize(slaveNodes.size());

  for(size_t i = 0; i < slaveNodes.size(); i++) {
    size_t elemId = slave2MasterElem[i];
    AMP::Mesh::LibMeshElement el = masterMeshAdapter->getElement(elemId);
    size_t bndSideId = 0;
    for(; bndSideId < el.numSides(); bndSideId++) {
      AMP::Mesh::LibMeshSide side = el.getSide(bndSideId, true);
      bool allNodesOnMasterBnd = true;
      for(size_t k = 0; k < side.numNodes(); k++) {
        size_t ndId = side.getNodeID(k);
        AMP::Mesh::LibMeshNode nd = masterMeshAdapter->getNode(ndId);
        if(!(masterMeshAdapter->isOnBoundary(nd, masterId))) {
          allNodesOnMasterBnd = false;
          break;
        }
      }//end for k
      if(allNodesOnMasterBnd) {
        for(size_t k = 0; k < side.numNodes(); k++) {
          size_t ndId = side.getNodeID(k);
          slave2MasterNodes[i].push_back(ndId);
        }//end for k
        AMP::Mesh::LibMeshNode nd = slaveMeshAdapter->getNode(slaveNodes[i]);
        ::Point pt(nd.x(), nd.y(), nd.z()); 
        ::Elem* elPtr = &(el.getElem());
        ::FEType fe_type(elPtr->default_order());
        ::Point ref_point = ::FEInterface::inverse_map(elPtr->dim(), fe_type, elPtr, pt, precision, true);
        assert(::FEInterface::on_reference_element(ref_point, elPtr->type(), precision));
        switch(bndSideId) {
          case 0 : {
                     ref_point(2) = -1;
                     break;
                   }
          case 1 : {
                     ref_point(1) = -1;
                     break;
                   }
          case 2 : {
                     ref_point(0) = 1;
                     break;
                   }
          case 3 : {
                     ref_point(1) = 1;
                     break;
                   }
          case 4 : {
                     ref_point(0) = -1;
                     break;
                   }
          case 5 : {
                     ref_point(2) = 1;
                     break;
                   }                   
          default :
                   AMP_ERROR("This should not happen.");
        }//end switch
        for(int j = 0; j < 3; j++) {
          if(ref_point(j) < -1) {
            ref_point(j) = -1;
          }
          if(ref_point(j) > 1) {
            ref_point(j) = 1;
          }
        }//end for j
        int xCoords[] = {-1, 1, 1, -1, -1, 1, 1, -1};
        int yCoords[] = {-1, -1, 1, 1, -1, -1, 1, 1};
        int zCoords[] = {-1, -1, -1, -1, 1, 1, 1, 1};
        for(size_t k = 0; k < side.numNodes(); k++) {
          unsigned int locNodeId = ::Hex8::side_nodes_map[bndSideId][k];
          double xFac = (2.0 - fabs(ref_point(0) - xCoords[locNodeId]))/2.0;
          double yFac = (2.0 - fabs(ref_point(1) - yCoords[locNodeId]))/2.0;
          double zFac = (2.0 - fabs(ref_point(2) - zCoords[locNodeId]))/2.0;
          slave2MasterFactors[i].push_back(xFac*yFac*zFac);
        }//end for k
        break;
      }
    }//end for bndSideId
    assert(bndSideId < el.numSides());
    //    std::cout<<"Slave node: "<<(slaveNodes[i])<<
    //      " is mapped to Side: "<<bndSideId<<" on Master element: "<<(slave2MasterElem[i])<<std::endl;
  }//end for i

}

void computeSlave2MasterElem(const unsigned int slaveId, const unsigned int masterId,
    const double precision, double* rgH, const unsigned int rgDim,
    AMP::Mesh::MeshManager::Adapter::shared_ptr slaveMeshAdapter,
    AMP::Mesh::MeshManager::Adapter::shared_ptr masterMeshAdapter, 
    double* minXYZ, double* maxXYZ, std::vector<std::vector<size_t> > const & rg2ElemMap, 
    std::vector<size_t> & slaveNodes, std::vector<size_t> & slave2MasterElem) {

  AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd =
    slaveMeshAdapter->beginOwnedBoundary( slaveId );
  AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd =
    slaveMeshAdapter->endOwnedBoundary( slaveId );

  unsigned int numNodesOnSlaveSurface = 0;
  for( ; bnd != end_bnd; ++bnd) {
    numNodesOnSlaveSurface++;
    double currXYZ[3];
    currXYZ[0] = bnd->x();
    currXYZ[1] = bnd->y();
    currXYZ[2] = bnd->z();
    ::Point currPt(currXYZ[0], currXYZ[1], currXYZ[2]);
    for(int j = 0; j < 3; j++) {
      if(currXYZ[j] <= minXYZ[j]) {
        std::cout<<"j = "<<j<<" currXYZ = "<<currXYZ[j]<<" minXYZ = "<<minXYZ[j]<<std::endl;
      }
      if(currXYZ[j] >= maxXYZ[j]) {
        std::cout<<"j = "<<j<<" currXYZ = "<<currXYZ[j]<<" maxXYZ = "<<maxXYZ[j]<<std::endl;
      }
      assert(currXYZ[j] > minXYZ[j]);
      assert(currXYZ[j] < maxXYZ[j]);
    }//end for j
    unsigned int xi = static_cast<unsigned int>(floor((currXYZ[0] - minXYZ[0])/rgH[0]));
    unsigned int yi = static_cast<unsigned int>(floor((currXYZ[1] - minXYZ[1])/rgH[1]));
    unsigned int zi = static_cast<unsigned int>(floor((currXYZ[2] - minXYZ[2])/rgH[2]));
    unsigned int rgId = (zi*rgDim*rgDim) + (yi*rgDim) + xi;
    std::vector<size_t> selectedElems;
    for(size_t k = 0; k < rg2ElemMap[rgId].size(); k++) {
      size_t elemId = rg2ElemMap[rgId][k];
      ::Elem* el = &((masterMeshAdapter->getElement(elemId)).getElem());
      if(myContainsPoint(el, currPt, precision)) {
        // std::cout<<"Slave node "<<(bnd->globalID())<<" is in Master element "<<(elemId)<<std::endl; 
        selectedElems.push_back(elemId);
      }
    }//end for k
    if(!(selectedElems.empty())) {
      slaveNodes.push_back(bnd->globalID());
      if(selectedElems.size() > 1) {
        std::cout<<"Slave node can be mapped to multiple master elements"<<std::endl;
        bool hasMasterNodes = false;
        for(size_t i = 0; i < selectedElems.size(); i++) {
          AMP::Mesh::LibMeshElement el = masterMeshAdapter->getElement(selectedElems[i]);
          for(size_t j = 0; j < el.numNodes(); j++) {
            size_t ndId = el.getNodeID(j);
            AMP::Mesh::LibMeshNode nd = masterMeshAdapter->getNode(ndId);
            if(masterMeshAdapter->isOnBoundary(nd, masterId)) {
              hasMasterNodes = true;
              break;
            }
          }//end for j
          if(hasMasterNodes) {
            slave2MasterElem.push_back(selectedElems[i]);
            break;
          }
        }//end for i
        assert(hasMasterNodes);
      } else {
        bool hasMasterNodes = false;
        AMP::Mesh::LibMeshElement el = masterMeshAdapter->getElement(selectedElems[0]);
        for(size_t j = 0; j < el.numNodes(); j++) {
          size_t ndId = el.getNodeID(j);
          AMP::Mesh::LibMeshNode nd = masterMeshAdapter->getNode(ndId);
          if(masterMeshAdapter->isOnBoundary(nd, masterId)) {
            hasMasterNodes = true;
            break;
          }
        }//end for j
        assert(hasMasterNodes);
        slave2MasterElem.push_back(selectedElems[0]);
      }
    } 
  }//end for bnd

  std::cout<<"Num Nodes on Slave Surface = "<<numNodesOnSlaveSurface<<std::endl;

  assert(slaveNodes.size() == slave2MasterElem.size());

  /*
     for(size_t i = 0; i < slaveNodes.size(); i++) {
     std::cout<<"Slave node: "<<(slaveNodes[i])<<
     " is mapped to Master element: "<<(slave2MasterElem[i])<<std::endl;
     }//end for i
     */

}

void computeRG2ElemMap(const double precision, const unsigned int rgDim, 
    AMP::Mesh::MeshManager::Adapter::shared_ptr masterMeshAdapter, double* minXYZ, double* maxXYZ,
    std::vector<std::vector<size_t> > & rg2ElemMap, double* rgH) {

  rg2ElemMap.resize(rgDim*rgDim*rgDim);

  for(int j = 0; j < 3; j++) {
    rgH[j] = (maxXYZ[j] - minXYZ[j])/(static_cast<double>(rgDim));
  }//end for j

  AMP::Mesh::MeshManager::Adapter::ElementIterator  el = masterMeshAdapter->beginElement();
  AMP::Mesh::MeshManager::Adapter::ElementIterator  end_el = masterMeshAdapter->endElement();

  for( ; el != end_el; ++el) {
    double currMinXYZ[3];
    double currMaxXYZ[3]; 
    for(size_t i = 0; i < el->numNodes(); i++) {
      size_t ndId = el->getNodeID(i);
      AMP::Mesh::MeshManager::Adapter::Node nd = masterMeshAdapter->getNode(ndId);
      if(i == 0) {
        currMinXYZ[0] = currMaxXYZ[0] = nd.x(); 
        currMinXYZ[1] = currMaxXYZ[1] = nd.y(); 
        currMinXYZ[2] = currMaxXYZ[2] = nd.z(); 
      } else {
        if(nd.x() < currMinXYZ[0]) {
          currMinXYZ[0] = nd.x();
        }
        if(nd.x() > currMaxXYZ[0]) {
          currMaxXYZ[0] = nd.x();
        }
        if(nd.y() < currMinXYZ[1]) {
          currMinXYZ[1] = nd.y();
        }
        if(nd.y() > currMaxXYZ[1]) {
          currMaxXYZ[1] = nd.y();
        }
        if(nd.z() < currMinXYZ[2]) {
          currMinXYZ[2] = nd.z();
        }
        if(nd.z() > currMaxXYZ[2]) {
          currMaxXYZ[2] = nd.z();
        }
      }
    }//end for i
    for(int j = 0; j < 3; j++) {
      assert(currMinXYZ[j] > minXYZ[j]);
      assert(currMaxXYZ[j] < maxXYZ[j]);
    }//end for j
    unsigned int xyzs[3];
    unsigned int xyze[3];
    for(int j = 0; j < 3; j++) {
      xyzs[j] = static_cast<unsigned int>(floor((currMinXYZ[j] - minXYZ[j])/rgH[j]));
      xyze[j] = static_cast<unsigned int>(ceil((currMaxXYZ[j] - minXYZ[j])/rgH[j]));
      assert(xyze[j] <= rgDim);
    }//end for j
    for(unsigned int zi = xyzs[2]; zi < xyze[2]; zi++) {
      for(unsigned int yi = xyzs[1]; yi < xyze[1]; yi++) {
        for(unsigned int xi = xyzs[0]; xi < xyze[0]; xi++) {
          unsigned int rgId = (zi*rgDim*rgDim) + (yi*rgDim) + xi;
          rg2ElemMap[rgId].push_back(el->globalID());
        }//end for xi
      }//end for yi
    }//end for zi
  }//end for el

  size_t maxElemsPerCell = 0;
  for(size_t i = 0; i < rg2ElemMap.size(); i++) {
    if((rg2ElemMap[i].size()) > maxElemsPerCell) {
      maxElemsPerCell = (rg2ElemMap[i].size());
    }
  }//end for i

  std::cout<<"Max Elements per RG Cell: "<<maxElemsPerCell<<std::endl;
}

void computeRGboundingBox(const double precision, 
    AMP::Mesh::MeshManager::Adapter::shared_ptr masterMeshAdapter, double* minXYZ, double* maxXYZ) {

  AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator  nd = masterMeshAdapter->beginOwnedNode();
  AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator  end_nd = masterMeshAdapter->endOwnedNode();

  minXYZ[0] = maxXYZ[0] = nd->x();
  minXYZ[1] = maxXYZ[1] = nd->y();
  minXYZ[2] = maxXYZ[2] = nd->z();
  for(; nd != end_nd; ++nd ) {
    if(nd->x() < minXYZ[0]) {
      minXYZ[0] = nd->x();
    }
    if(nd->x() > maxXYZ[0]) {
      maxXYZ[0] = nd->x();
    }
    if(nd->y() < minXYZ[1]) {
      minXYZ[1] = nd->y();
    }
    if(nd->y() > maxXYZ[1]) {
      maxXYZ[1] = nd->y();
    }
    if(nd->z() < minXYZ[2]) {
      minXYZ[2] = nd->z();
    }
    if(nd->z() > maxXYZ[2]) {
      maxXYZ[2] = nd->z();
    }
  }//end for nd

  for(int j = 0; j < 3; j++) {
    minXYZ[j] = minXYZ[j] - precision;
    maxXYZ[j] = maxXYZ[j] + precision;
  }//end for j

}



