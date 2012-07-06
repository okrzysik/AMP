#include "operators/contact/NodeToSegmentConstraintsOperator.h"
#include "ampmesh/dendro/DendroSearch.h"
#include "ampmesh/MeshID.h"

#include <limits>
#include <algorithm>
#include <vector>

namespace AMP {
  namespace Operator {

void NodeToSegmentConstraintsOperator::reset(const boost::shared_ptr<OperatorParameters> & params) {

  AMP_INSIST( (params != NULL), "NULL parameter" );
  AMP_INSIST( ((params->d_db) != NULL), "NULL database" );

  AMP::Mesh::Mesh::shared_ptr mesh = params->d_Mesh;
  AMP::AMP_MPI comm = mesh->getComm(); 

  AMP::Mesh::Mesh::shared_ptr masterMesh = mesh->Subset(AMP::Mesh::MeshID(d_MasterMeshID));
  AMP::Mesh::Mesh::shared_ptr slaveMesh = mesh->Subset(AMP::Mesh::MeshID(d_SlaveMeshID));

  /** get the boundary slave vertices coordinates */
  AMP::Mesh::MeshIterator slaveMeshIterator = slaveMesh->getBoundaryIDIterator(AMP::Mesh::Vertex, d_SlaveBoundaryID);
  AMP::Mesh::MeshIterator slaveMeshIterator_begin = slaveMeshIterator.begin(), 
      slaveMeshIterator_end = slaveMeshIterator.end();

  size_t nSlaveVertices = slaveMeshIterator.size();

  const unsigned int dofsPerNode = 3;
  std::vector<size_t> tmpSlaveNodeIDs(dofsPerNode);
  std::vector<size_t> slaveNodeIDs(dofsPerNode*nSlaveVertices);
  std::vector<size_t>::iterator slaveNodeIDsIterator = slaveNodeIDs.begin();

  std::vector<double> tmpSlaveVertexCoord(3);
  std::vector<double> slaveVerticesCoord(3*nSlaveVertices);
  std::vector<double>::iterator slaveVerticesCoordIterator = slaveVerticesCoord.begin();

  for (slaveMeshIterator = slaveMeshIterator_begin; slaveMeshIterator != slaveMeshIterator_end; ++slaveMeshIterator) {
    d_DOFmanager->getDOFs(slaveMeshIterator->globalID(), tmpSlaveNodeIDs);
    AMP_ASSERT( tmpSlaveNodeIDs.size() == dofsPerNode );
    std::copy(tmpSlaveNodeIDs.begin(), tmpSlaveNodeIDs.end(), slaveNodeIDsIterator);
    for (size_t i = 0; i < dofsPerNode; ++i) { ++slaveNodeIDsIterator; }
    tmpSlaveVertexCoord = slaveMeshIterator->coord();
    AMP_ASSERT( tmpSlaveVertexCoord.size() == 3 );
    std::copy(tmpSlaveVertexCoord.begin(), tmpSlaveVertexCoord.end(), slaveVerticesCoordIterator);
    for (size_t i = 0; i < 3; ++i) { ++slaveVerticesCoordIterator; }
  } // end loop over the slave vertices on boundary

  AMP_ASSERT( slaveNodeIDsIterator == slaveNodeIDs.end() );
  AMP_ASSERT( slaveVerticesCoordIterator == slaveVerticesCoord.end() );

  /** do a dendro search for the boundary slave vertices on the master mesh */
  DendroSearch dendroSearchOnMaster(comm, masterMesh);
  dendroSearchOnMaster.search(slaveVerticesCoord);

  std::vector<double> slaveVerticesProjLocalCoord;
  std::vector<size_t> masterNodeIDs, masterNodeOwnerRanks;
  std::vector<int> flags;

//  dendroSearchOnMaster.projectOnBoundaryID(d_MasterBoundaryID, dofsPerNode, d_DOFmanager,
//      masterNodeOwnerRanks, masterNodeIDs, slaveVerticesProjLocalCoord, flags);

  AMP_ASSERT( nSlaveVertices == flags.size() );
  AMP_ASSERT( nSlaveVertices == slaveNodeIDs.size() / dofsPerNode );
  AMP_ASSERT( nSlaveVertices == masterNodeIDs.size() / (4 * dofsPerNode) );
  AMP_ASSERT( nSlaveVertices == masterNodeOwnerRanks.size() );
  AMP_ASSERT( nSlaveVertices == slaveVerticesProjLocalCoord.size() / 2 );

  /** build the constraints */
  const unsigned int nConstraints = std::count(flags.begin(), flags.end(), DendroSearch::FoundOnBoundary);
  AMP_ASSERT( std::count(flags.begin(), flags.end(), DendroSearch::FoundNotOnBoundary) == 0 ); // DendroSearch::FoundNotOnBoundary is not acceptable

  size_t npes = comm.getSize();
  std::vector<size_t> row, col, ranks;
  std::vector<double> val;
  ranks.resize(nConstraints);
  std::fill(ranks.begin(), ranks.end(), npes);
  row.resize(4*nConstraints);
  std::fill(row.begin(), row.end(), std::numeric_limits<unsigned int>::max());
  col.resize(4*nConstraints);
  std::fill(col.begin(), col.end(), std::numeric_limits<unsigned int>::max());
  val.resize(4*nConstraints);
  std::fill(val.begin(), val.end(), 0.0);

  slaveNodeIDsIterator = slaveNodeIDs.begin();
  std::vector<size_t>::const_iterator masterNodeIDsIterator = masterNodeIDs.begin();
  std::vector<size_t>::iterator rowIterator = row.begin(), colIterator = col.begin();
  std::vector<double>::iterator valIterator = val.begin();
  double const * slaveVerticesProjLocalCoord_ptr = &(slaveVerticesProjLocalCoord[0]);
  double basis_functions_values_on_face[4];
  std::vector<int>::const_iterator flagsIterator = flags.begin(),
      flagsIterator_end = flags.end();
  for ( ; flagsIterator != flagsIterator_end; ++flagsIterator) {
    AMP_ASSERT( (*flagsIterator == DendroSearch::NotFound) || (*flagsIterator == DendroSearch::FoundOnBoundary) );
    if (*flagsIterator == DendroSearch::FoundOnBoundary) {
      hex8_element_t::get_basis_functions_values_on_face(slaveVerticesProjLocalCoord_ptr, basis_functions_values_on_face);
      for (size_t i = 0; i < 2; ++i) { ++slaveVerticesProjLocalCoord_ptr; }
      for (unsigned int n = 0; n < dofsPerNode; ++n) {
        for (size_t i = 0; i < 4; ++i) {
          *rowIterator = *slaveNodeIDsIterator;
          *colIterator = *masterNodeIDsIterator;
          *valIterator = basis_functions_values_on_face[i];
          ++rowIterator;
          ++colIterator;
          ++valIterator;
          ++masterNodeIDsIterator;
        } // end for i
        ++slaveNodeIDsIterator;
      } // end for n
    } else {
      for (size_t i = 0; i < 2; ++i) { ++slaveVerticesProjLocalCoord_ptr; }
      for (unsigned int n = 0; n < dofsPerNode; ++n) {
        for (size_t i = 0; i < 4; ++i) {
          ++rowIterator;
          ++colIterator;
          ++valIterator;
          ++masterNodeIDsIterator;
        } // end for i
        ++slaveNodeIDsIterator;
      } // end for n
    } // end if
  } // end for
  AMP_ASSERT( slaveVerticesProjLocalCoord_ptr == &(slaveVerticesProjLocalCoord[0])+2*nSlaveVertices );
  AMP_ASSERT( masterNodeIDsIterator == masterNodeIDs.end() );
  AMP_ASSERT( slaveNodeIDsIterator == slaveNodeIDs.end() );
  AMP_ASSERT( rowIterator == row.end() );
  AMP_ASSERT( colIterator == col.end() );
  AMP_ASSERT( valIterator == val.end() );

  /** setup for apply */
  d_SendCnts.resize(npes);
  std::fill(d_SendCnts.begin(), d_SendCnts.end(), 0);
  for (size_t i = 0; i < nConstraints; ++i) {
    ++d_SendCnts[ranks[i]]; 
  } // end for i
  for (size_t i = 0; i < npes; ++i) {
    d_SendCnts[i] *= 5; 
  } // end for i

  d_SendDisps.resize(npes);
  d_SendDisps[0] = 0;
  for (size_t i = 1; i < npes; ++i) {
    d_SendDisps[i] = d_SendDisps[i-1] + d_SendCnts[i-1]; 
  } // end for i
  AMP_ASSERT( d_SendDisps[npes-1] == 5*nConstraints );

  std::vector<int> tmpSendCnts(npes, 0);
  std::vector<size_t> sendBuff(d_SendDisps[npes-1]+d_SendCnts[npes-1], 0.0);
  for (size_t i = 0; i < nConstraints; ++i) {
    size_t sendToRank = ranks[i];
    sendBuff[d_SendDisps[sendToRank] + tmpSendCnts[sendToRank]] = i;
    ++tmpSendCnts[sendToRank];
    for (size_t j = 0; j < 4; ++j) {
      sendBuff[d_SendDisps[sendToRank] + tmpSendCnts[sendToRank]] = col[4*i+j];
      ++tmpSendCnts[sendToRank];
    } // end for j
  } // end for i 
  AMP_ASSERT( std::equal(tmpSendCnts.begin(), tmpSendCnts.end(), d_SendCnts.end()) ); 

  d_RecvCnts.resize(npes);
  comm.allToAll(1, &(d_SendCnts[0]), &(d_RecvCnts[0]));
  d_RecvDisps.resize(npes);
  d_RecvDisps[0] = 0;
  for (size_t i = 1; i < npes; ++i) {
    d_RecvDisps[i] = d_RecvDisps[i-1] + d_RecvCnts[i-1];
  } // end for i
  d_RecvMasterGlobalIDs.resize(d_RecvDisps[npes-1]+d_RecvCnts[npes-1]);
  comm.allToAll((!(sendBuff.empty()) ? &(sendBuff[0]) : NULL), &(d_SendCnts[0]), &(d_SendDisps[0]),
      (!(d_RecvMasterGlobalIDs.empty()) ? &(d_RecvMasterGlobalIDs[0]) : NULL), &(d_RecvCnts[0]), &(d_RecvDisps[0]), true);

  for (size_t i = 0; i < npes; ++i) {
    d_TransposeSendCnts[i] = (d_SendCnts[i] / 5) * 6;
    d_TransposeSendDisps[i] = (d_SendDisps[i] / 5) * 6;
    d_TransposeRecvCnts[i] = (d_RecvCnts[i] / 5) * 6;
    d_TransposeRecvDisps[i] = (d_RecvDisps[i] / 5) * 6;
  } // end for i
  std::swap_ranges(d_RecvCnts.begin(), d_RecvCnts.end(), d_SendCnts.begin());
  std::swap_ranges(d_RecvDisps.begin(), d_RecvDisps.end(), d_SendDisps.begin());

}

void NodeToSegmentConstraintsOperator::apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
    AMP::LinearAlgebra::Vector::shared_ptr &r, const double a, const double b) {

  /** send and receive the master values */
  AMP::AMP_MPI comm = u->getComm();
  size_t npes = comm.getSize();
//  size_t rank = comm.getRank();

  std::vector<double> sendMasterValues(d_SendDisps[npes-1]+d_SendCnts[npes-1]);
  for (size_t i = 0; i < npes; ++i) {
    for (size_t j = 0; j < d_SendCnts[i]; j += 5) {
      size_t k = d_SendDisps[i] + j;
      sendMasterValues[k] = static_cast<double>(d_RecvMasterGlobalIDs[k]);
      u->getValuesByGlobalID(4, &(d_RecvMasterGlobalIDs[k+1]), &(sendMasterValues[k+1]));
    } // end for j
  } // end for i

  std::vector<double> recvMasterValues(d_RecvDisps[npes-1]+d_RecvCnts[npes-1]);
  comm.allToAll((!(sendMasterValues.empty()) ? &(sendMasterValues[0]) : NULL), &(d_SendCnts[0]), &(d_SendDisps[0]),
      (!(recvMasterValues.empty()) ? &(recvMasterValues[0]) : NULL), &(d_RecvCnts[0]), &(d_RecvDisps[0]), true);

  /** compute slave values */
  for (size_t i = 0; i < npes; ++i) {
    for (size_t j = 0; j < d_RecvCnts[i]; j += 5) {
      size_t k = d_RecvDisps[i] + j;
      size_t constraintsLocalIndex = static_cast<size_t>(recvMasterValues[k]);
      size_t slaveGlobalID = d_SlaveGlobalIDs[constraintsLocalIndex];
      double value = 0.0;
      for (size_t l = 0; l < 4; ++l) {
        value += d_MasterShapeFunctionsValues[4*constraintsLocalIndex+l] * recvMasterValues[k+1+l]; 
      } // end for j
      u->setValueByGlobalID(slaveGlobalID, value);
    } // end for j
  } // end for i

}

void NodeToSegmentConstraintsOperator::applyTranspose(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
    AMP::LinearAlgebra::Vector::shared_ptr  &r, const double a, const double b) {

  /** send and receive the slave values */
  AMP::AMP_MPI comm = u->getComm();
  size_t npes = comm.getSize();

/*

  std::vector<double> sendSlaveValueAndShapeFunctionsValues(d_TransposeSendDisps[npes-1]+d_TransposeSendCnts[npes-1]);
  for (size_t i = 0; i < npes; ++i) {
    for (size_t j = 0; j < d_TransposeSendCnts[i]; j += 6) {
      size_t k = d_TransposeSendDisps[i] + j;
      sendSlaveValueAndShapeFunctionsValues[k] = static_cast<double>(d_RecvMasterGlobalIDs[k]);
      u->getValueByGlobalID(&(d_SlaveGlobalIDs[k+1]), &(sendMasterValues[k+1]));
    } // end for j
  } // end for i
  */
}

  } // end namespace Operator
} // end namespace AMP
