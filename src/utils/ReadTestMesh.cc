#ifdef USE_AMP_MESH
#include "ReadTestMesh.h"

#include "elem.h"
#include "cell_hex8.h"
#include "boundary_info.h"

#include <cstdio>

namespace AMP {

  void readBinaryTestMesh(std::string mesh_file, boost::shared_ptr< ::Mesh> mesh) {
    FILE* fp = fopen(mesh_file.c_str(), "rb");

    int num_nodes, n;
    n = fread(&num_nodes, sizeof(int), 1, fp);

    mesh->reserve_nodes(num_nodes);

    std::vector<double> points(3*num_nodes);

    n = fread(&(points[0]), sizeof(double), (3*num_nodes), fp);

    for(int i = 0; i < num_nodes; i++) {
      mesh->add_point(::Point(points[(3*i) + 0], points[(3*i) + 1], points[(3*i) + 2]), i);
    }

    points.clear();

    int num_elem;
    n = fread(&num_elem, sizeof(int), 1, fp);

    mesh->reserve_elem(num_elem);

    std::vector<int> elemNodeMap(8*num_elem);

    n = fread(&(elemNodeMap[0]), sizeof(int), (8*num_elem), fp);

    for(int i = 0; i < num_elem; i++) {
      ::Elem* newElem = new ::Hex8;
      newElem->set_id(i);
      ::Elem* elem = mesh->add_elem(newElem);
      for(int j = 0; j < 8; j++) {
        elem->set_node(j) = mesh->node_ptr(elemNodeMap[(8*i) + j]);
      }
    }

    elemNodeMap.clear();

    int num_boundaryNodeIds;
    n = fread(&num_boundaryNodeIds, sizeof(int), 1, fp);

    for(int bid = 1; bid <= num_boundaryNodeIds; bid++) {
      int bndDofSize;
      n = fread(&bndDofSize, sizeof(int), 1, fp);
      if(bndDofSize > 0) {
        int idx;
        for(int i = 0; i < bndDofSize; i++) {
          n = fread(&idx, sizeof(int), 1, fp);
          mesh->boundary_info->add_node(mesh->node_ptr(idx), bid);
        }
      }
    } 

    int num_boundarySideIds;
    n = fread(&num_boundarySideIds, sizeof(int), 1, fp);

    for(int bid = 1; bid <= num_boundarySideIds; bid++) {
      int bndDofSize;
      n = fread(&bndDofSize, sizeof(int), 1, fp);
      if(bndDofSize > 0) {
        int idxE;
        int idxS;
        for(int i = 0; i < bndDofSize; i++) {
          n = fread(&idxE, sizeof(int), 1, fp);
          n = fread(&idxS, sizeof(int), 1, fp);
          mesh->boundary_info->add_side(mesh->elem(idxE), idxS, bid);
        }
      }
    } 

    AMP_INSIST(n==0, "fread returned non-zero error code");
    fclose(fp);
  }

  void readTestMesh(std::string mesh_file, boost::shared_ptr< ::Mesh > mesh) {
    FILE* fp = fopen(mesh_file.c_str(), "r");
    char str[256];
    int n;

    //Mesh {
    n = fscanf(fp, "%s {", str);

    int num_nodes;
    //NumberOfNodes = <val>
    n = fscanf(fp, "%s = %d", str, &num_nodes);

    mesh->reserve_nodes(num_nodes);

    for(int i = 0; i < num_nodes; i++) {
      double point[3];
      //PointK = x, y, z
      n = fscanf(fp, "%s = %lf, %lf, %lf", str, &(point[0]), &(point[1]), &(point[2]));
      mesh->add_point(::Point(point[0], point[1], point[2]), i);
    }//end for i

    int num_elem;
    //NumberOfElements = <val>
    n = fscanf(fp, "%s = %d", str, &num_elem);

    mesh->reserve_elem(num_elem);

    std::vector<std::vector<int> > elemNodeMap;

    for(int i = 0; i < num_elem; i++) {
      std::vector<int> nodesForElem(8);
      //ElemK = i0, i1, i2, i3, i4, i5, i6, i7
      n = fscanf(fp, "%s = %d,", str, &(nodesForElem[0]));
      for(int j = 1; j < 7; j++) {
        n = fscanf(fp, "%d,", &(nodesForElem[j]));
      }//end for j
      n = fscanf(fp, "%d", &(nodesForElem[7]));
      elemNodeMap.push_back(nodesForElem);
    }//end for i

    for(int i = 0; i < num_elem; i++) {
      ::Elem* newElem = new ::Hex8;
      newElem->set_id(i);
      ::Elem* elem = mesh->add_elem(newElem);
      for(int j = 0; j < 8; j++) {
        elem->set_node(j) = mesh->node_ptr(elemNodeMap[i][j]);
      }//end for j
    }//end for i

    int num_boundaryNodeIds;
    //NumberOfBoundaryNodeIds = <val>
    n = fscanf(fp, "%s = %d", str, &num_boundaryNodeIds);

    for(int bid = 1; bid <= num_boundaryNodeIds; bid++) {
      int bndDofSize;
      //NumberOfBoundaryNodesK = <val> 
      n = fscanf(fp, "%s = %d", str, &bndDofSize);
      if(bndDofSize > 0) {
        //BoundaryNodeIdK =
        n = fscanf(fp, "%s =", str);
        int idx;
        for(int i = 0; i < (bndDofSize - 1); i++) {
          n = fscanf(fp, "%d,", &idx);
          mesh->boundary_info->add_node(mesh->node_ptr(idx), bid);
        }//end for i
        n = fscanf(fp, "%d", &idx);
        mesh->boundary_info->add_node(mesh->node_ptr(idx), bid);
      }
    }//end for bid 

    int num_boundarySideIds;
    //NumberOfBoundarySideIds = <val>
    n = fscanf(fp, "%s = %d", str, &num_boundarySideIds);

    for(int bid = 1; bid <= num_boundarySideIds; bid++) {
      int bndDofSize;
      //NumberOfBoundarySidesK = <val> 
      n = fscanf(fp, "%s = %d", str, &bndDofSize);
      if(bndDofSize > 0) {
        //BoundarySideIdK =
        n = fscanf(fp, "%s =", str);
        int idxE;
        int idxS;
        for(int i = 0; i < (bndDofSize - 1); i++) {
          n = fscanf(fp, "%d,", &idxE);
          n = fscanf(fp, "%d,", &idxS);
          mesh->boundary_info->add_side(mesh->elem(idxE), idxS, bid);
        }//end for i
        n = fscanf(fp, "%d,", &idxE);
        n = fscanf(fp, "%d", &idxS);
        mesh->boundary_info->add_side(mesh->elem(idxE), idxS, bid);
      }
    }//end for bid 

    AMP_INSIST(n==0, "fread returned non-zero error code");
    fclose(fp);
  }

  void readTestMesh(boost::shared_ptr<AMP::InputDatabase> mesh_file_db, 
      boost::shared_ptr< ::Mesh > mesh) {
    boost::shared_ptr<AMP::Database> mesh_db = mesh_file_db->getDatabase("Mesh");
    int num_elem = mesh_db->getInteger("NumberOfElements");
    int num_nodes = mesh_db->getInteger("NumberOfNodes");
    int num_boundaryNodeIds = mesh_db->getInteger("NumberOfBoundaryNodeIds");
    int num_boundarySideIds = mesh_db->getInteger("NumberOfBoundarySideIds");

    mesh->reserve_elem(num_elem);
    mesh->reserve_nodes(num_nodes);

    for(int i = 0; i < num_nodes; i++) {
      char key[100];
      sprintf(key, "Point%d", i);
      double point[3];
      mesh_db->getDoubleArray(key, point, 3);
      mesh->add_point(::Point(point[0], point[1], point[2]), i);
    }//end for i

    std::vector<std::vector<int> > elemNodeMap;

    for(int i = 0; i < num_elem; i++) {
      char key[100];
      sprintf(key, "Elem%d", i);
      int nodesForElem[8];
      mesh_db->getIntegerArray(key, nodesForElem, 8);
      std::vector<int> tmpArr(8);
      for(int j = 0; j < 8; j++) {
        tmpArr[j] = nodesForElem[j];
      }//end for j
      elemNodeMap.push_back(tmpArr);
    }//end for i

    for(int i = 0; i < num_elem; i++) {
      ::Elem* elem = mesh->add_elem(new ::Hex8);
      for(int j = 0; j < 8; j++) {
        elem->set_node(j) = mesh->node_ptr(elemNodeMap[i][j]);
      }//end for j
    }//end for i

    for(int bid = 1; bid <= num_boundaryNodeIds; bid++) {
      char key[100];
      sprintf(key, "BoundaryNodeId%d", bid);
      char newKey[100];
      sprintf(newKey, "NumberOfBoundaryNodes%d", bid);
      int bndDofSize = mesh_db->getInteger(newKey);
      if(bndDofSize > 0) {
        int* bndDofIndices = new int[bndDofSize];
        mesh_db->getIntegerArray(key, bndDofIndices, bndDofSize);
        for(int i = 0; i < bndDofSize; i++) {
          mesh->boundary_info->add_node(mesh->node_ptr(bndDofIndices[i]), bid);
        }//end for i
        delete [] bndDofIndices;
      }
    }//end for bid 

    for(int bid = 1; bid <= num_boundarySideIds; bid++) {
      char key[100];
      sprintf(key, "BoundarySideId%d", bid);
      char newKey[100];
      sprintf(newKey, "NumberOfBoundarySides%d", bid);
      int bndDofSize = mesh_db->getInteger(newKey);
      if(bndDofSize > 0) {
        int* bndDofIndices = new int[2*bndDofSize];
        mesh_db->getIntegerArray(key, bndDofIndices, (2*bndDofSize));
        for(int i = 0; i < bndDofSize; i++) {
          mesh->boundary_info->add_side(mesh->elem(bndDofIndices[2*i]), bndDofIndices[(2*i) + 1], bid);
        }//end for i
        delete [] bndDofIndices;
      }
    }//end for bid 

  }

  }

#endif

