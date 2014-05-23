
#include "operators/boundary/libmesh/PressureBoundaryOperator.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"
#include "vectors/VectorBuilder.h"
#include "discretization/simpleDOF_Manager.h"
#include "ampmesh/Mesh.h"

#include "enum_order.h"
#include "enum_fe_family.h"
#include "enum_quadrature_type.h"
#include "string_to_enum.h"
#include "auto_ptr.h"
#include "fe_type.h"
#include "fe_base.h"
#include "quadrature.h"
#include "cell_hex8.h"
#include "node.h"

namespace AMP {
  namespace Operator {

    PressureBoundaryOperator :: PressureBoundaryOperator(const boost::shared_ptr<OperatorParameters> & params)
      : BoundaryOperator(params) {
        AMP_ASSERT((params->d_db)->keyExists("BoundaryID"));
        short int bndId = (params->d_db)->getInteger("BoundaryID");

        //ASSUMPTION: Each boundary face element is associated with an unique volume element. This will be true
        //if the face element is truly on the boundary of the mesh.
        //NOTE: The processor that owns a face element may not own the
        //corresponding volume element.  

        AMP::AMP_MPI meshComm = d_Mesh->getComm();
        const int npes = meshComm.getSize();

        std::vector<std::vector<double> > volElemMap(npes);
        std::vector<std::vector<unsigned int> > sideMap(npes);
        std::vector<std::vector<AMP::Mesh::MeshElementID> > idMap(npes);

        AMP::Mesh::MeshIterator el = d_Mesh->getIterator(AMP::Mesh::Volume, 0);
        AMP::Mesh::MeshIterator end_el = el.end();
        for(; el != end_el; ++el) {
          std::vector<AMP::Mesh::MeshElement> sides = el->getElements(AMP::Mesh::Face);
          for(size_t s = 0; s < sides.size(); ++s) {
            if(sides[s].isOnBoundary(bndId)) {
              AMP::Mesh::MeshElementID sideId = sides[s].globalID();
              unsigned int owner = sideId.owner_rank();
              std::vector<AMP::Mesh::MeshElement> vertices = el->getElements(AMP::Mesh::Vertex);
              for(size_t j = 0; j < vertices.size(); ++j) {
                std::vector<double> pt = vertices[j].coord();
                for(size_t k = 0; k < pt.size(); ++k) {
                  volElemMap[owner].push_back(pt[k]);
                }//end k
                AMP::Mesh::MeshElementID nodeId = vertices[j].globalID();
                idMap[owner].push_back(nodeId);
              }//end for j
              sideMap[owner].push_back(s);
            }
          }//end s
        }//end el

        std::vector<int> sendCnts(npes);
        for(int i = 0; i < npes; ++i) {
          sendCnts[i] = sideMap[i].size();
        }//end i

        std::vector<int> recvCnts(npes);
        meshComm.allToAll<int>(1, &(sendCnts[0]), &(recvCnts[0]));

        std::vector<int> sendDisps(npes);
        std::vector<int> recvDisps(npes);
        sendDisps[0] = 0;
        recvDisps[0] = 0;
        for(int i = 1; i < npes; ++i) {
          sendDisps[i] = sendDisps[i - 1] + sendCnts[i - 1];
          recvDisps[i] = recvDisps[i - 1] + recvCnts[i - 1];
        }//end i

        std::vector<double> sendVolElemList;
        std::vector<unsigned int> sendSideList;
        std::vector<AMP::Mesh::MeshElementID> sendIDlist;
        for(int i = 0; i < npes; ++i) {
          sendSideList.insert(sendSideList.end(), sideMap[i].begin(), sideMap[i].end());
          sendIDlist.insert(sendIDlist.end(), idMap[i].begin(), idMap[i].end());
          sendVolElemList.insert(sendVolElemList.end(), volElemMap[i].begin(), volElemMap[i].end());
          volElemMap[i].clear();
          sideMap[i].clear();
          idMap[i].clear();
        }//end i
        volElemMap.clear();
        sideMap.clear();
        idMap.clear();

        std::vector<unsigned int> recvSideList(recvDisps[npes - 1] + recvCnts[npes - 1]);

        unsigned int* sendSideListPtr = NULL;
        unsigned int* recvSideListPtr = NULL;
        if(!(sendSideList.empty())) {
          sendSideListPtr = &(sendSideList[0]);
        }
        if(!(recvSideList.empty())) {
          recvSideListPtr = &(recvSideList[0]);
        }
        meshComm.allToAll<unsigned int>(sendSideListPtr, &(sendCnts[0]), &(sendDisps[0]), 
            recvSideListPtr, &(recvCnts[0]), &(recvDisps[0]), true);
        sendSideList.clear();

        for(int i = 0; i < npes; ++i) {
          sendCnts[i] *= 8;
          sendDisps[i] *= 8;
          recvCnts[i] *= 8;
          recvDisps[i] *= 8;
        }//end i

        std::vector<AMP::Mesh::MeshElementID> recvIDlist(8*(recvSideList.size()));

        AMP::Mesh::MeshElementID* sendIDlistPtr = NULL;
        AMP::Mesh::MeshElementID* recvIDlistPtr = NULL;
        if(!(sendIDlist.empty())) {
          sendIDlistPtr = &(sendIDlist[0]);
        }
        if(!(recvIDlist.empty())) {
          recvIDlistPtr = &(recvIDlist[0]);
        }
        meshComm.allToAll<AMP::Mesh::MeshElementID>(sendIDlistPtr, &(sendCnts[0]), &(sendDisps[0]), 
            recvIDlistPtr, &(recvCnts[0]), &(recvDisps[0]), true);
        sendIDlist.clear();

        for(int i = 0; i < npes; ++i) {
          sendCnts[i] *= 3;
          sendDisps[i] *= 3;
          recvCnts[i] *= 3;
          recvDisps[i] *= 3;
        }//end i

        std::vector<double> recvVolElemList(24*(recvSideList.size()));

        double* sendVolElemListPtr = NULL;
        double* recvVolElemListPtr = NULL;
        if(!(sendVolElemList.empty())) {
          sendVolElemListPtr = &(sendVolElemList[0]);
        }
        if(!(recvVolElemList.empty())) {
          recvVolElemListPtr = &(recvVolElemList[0]);
        }
        meshComm.allToAll<double>(sendVolElemListPtr, &(sendCnts[0]), &(sendDisps[0]), 
            recvVolElemListPtr, &(recvCnts[0]), &(recvDisps[0]), true);
        sendVolElemList.clear();
        sendCnts.clear();
        recvCnts.clear();
        sendDisps.clear();
        recvDisps.clear();

        libMeshEnums::Order feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>("FIRST");
        libMeshEnums::FEFamily feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>("LAGRANGE");
        libMeshEnums::QuadratureType qruleType = Utility::string_to_enum<libMeshEnums::QuadratureType>("QGAUSS");
        boost::shared_ptr < ::FEType > feType ( new ::FEType(feTypeOrder, feFamily) );
        libMeshEnums::Order qruleOrder = feType->default_quadrature_order();
        boost::shared_ptr < ::QBase > qrule( (::QBase::build(qruleType, 2, qruleOrder)).release() );

        AMP_ASSERT((params->d_db)->keyExists("Value"));
        const double val = (params->d_db)->getDouble("Value");

        std::vector<double> pressure(12*(recvSideList.size()));

        for(size_t i = 0; i < recvSideList.size(); ++i) {
          ::Elem* elem = new ::Hex8;
          for(int j = 0; j < 8; ++j) {
            elem->set_node(j) = new ::Node(recvVolElemList[(24*i) + (3*j) + 0],
                recvVolElemList[(24*i) + (3*j) + 1], recvVolElemList[(24*i) + (3*j) + 2], j);
          }//end j

          boost::shared_ptr < ::FEBase > fe( (::FEBase::build(3, (*feType))).release() );
          fe->attach_quadrature_rule( qrule.get() );
          fe->reinit(elem, recvSideList[i]);

          AMP_ASSERT(qrule->n_points() == 4);

          const std::vector<Point>& normals = fe->get_normals();
          AMP_ASSERT(normals.size() == 4);

          for(size_t qp = 0; qp < qrule->n_points(); ++qp) {
            for(int d = 0; d < 3; ++d) {
              double tractionVal = val*(normals[qp](d));
              pressure[(12*i) + (3*qp) + d] = tractionVal;
            }//end d
          }//end qp

          for(size_t j = 0; j < elem->n_nodes(); ++j) {
            delete (elem->get_node(j));
            elem->set_node(j) = NULL;
          }//end for j
          delete elem;
          elem = NULL;
        }//end i

        boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));
        AMP_ASSERT((params->d_db)->keyExists("Variable"));
        std::string varName = params->d_db->getString("Variable");
        tmp_db->putString("Variable", varName);
        AMP_ASSERT((params->d_db)->keyExists("ResidualMode"));
        tmp_db->putBool("ResidualMode", ((params->d_db)->getBool("ResidualMode")));

        boost::shared_ptr<TractionBoundaryOperatorParameters> tracOpParams(new 
            TractionBoundaryOperatorParameters(tmp_db));
        tracOpParams->d_Mesh = d_Mesh;
        tracOpParams->d_traction = pressure;
        tracOpParams->d_volumeElements = recvVolElemList;
        tracOpParams->d_sideNumbers = recvSideList;
        tracOpParams->d_nodeID = recvIDlist;

        d_tractionOp.reset(new TractionBoundaryOperator(tracOpParams));        
      }

  }
}


