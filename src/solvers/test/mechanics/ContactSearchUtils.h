
#ifndef included_AMP_ContactSearchUtils
#define included_AMP_ContactSearchUtils

#include <vector>

static void computeRGboundingBox( const double precision,
                                  AMP::Mesh::MeshManager::Adapter::shared_ptr masterMeshAdapter,
                                  double *minXYZ,
                                  double *maxXYZ );

static void computeRG2ElemMap( const double precision,
                               const unsigned int rgDim,
                               AMP::Mesh::MeshManager::Adapter::shared_ptr masterMeshAdapter,
                               double *minXYZ,
                               double *maxXYZ,
                               std::vector<std::vector<size_t>> &rg2ElemMap,
                               double *rgH );

static void computeSlave2MasterElem( const unsigned int slaveId,
                                     const unsigned int masterId,
                                     const double precision,
                                     double *rgH,
                                     const unsigned int rgDim,
                                     AMP::Mesh::MeshManager::Adapter::shared_ptr slaveMeshAdapter,
                                     AMP::Mesh::MeshManager::Adapter::shared_ptr masterMeshAdapter,
                                     double *minXYZ,
                                     double *maxXYZ,
                                     std::vector<std::vector<size_t>> const &rg2ElemMap,
                                     std::vector<size_t> &slaveNodes,
                                     std::vector<size_t> &slave2MasterElem );

static void computeSlave2MasterNodes( const double precision,
                                      const unsigned int slaveId,
                                      const unsigned int masterId,
                                      AMP::Mesh::MeshManager::Adapter::shared_ptr slaveMeshAdapter,
                                      AMP::Mesh::MeshManager::Adapter::shared_ptr masterMeshAdapter,
                                      std::vector<size_t> const &slaveNodes,
                                      std::vector<size_t> const &slave2MasterElem,
                                      std::vector<std::vector<size_t>> &slave2MasterNodes,
                                      std::vector<std::vector<double>> &slave2MasterFactors );

bool myContainsPoint(libMesh::Elem *e, const ::Point &p, double tol );

#include "ContactSearchUtils.hpp"

#endif
