#ifndef included_AMP_Unit_test_Mesh_Generators_h
#define included_AMP_Unit_test_Mesh_Generators_h

#include "libmesh_config.h"
#include "libmesh.h"
#include "mesh.h"
#include "cell_hex8.h"
#include "boundary_info.h"

#include "ampmesh/MeshAdapter.h"
#include "ampmesh/LibMeshAdapter.h"
#include "ampmesh/MeshManager.h"
#include <exodusII_io.h>

namespace AMP {
namespace unit_test {

// Base class for Mesh Generators
class MeshGenerator
{
public:
    // Routine to build the mesh
    virtual void build_mesh() { AMP_ERROR("ERROR"); }
    // Routine to get the pointer to the mesh
    virtual AMP::Mesh::MeshAdapter::shared_ptr getMesh() {
        if ( mesh.get()==NULL ) 
            this->build_mesh();
        return mesh;
    }
protected:
    AMP::Mesh::MeshAdapter::shared_ptr  mesh;
};


// Class to create a cube in Libmesh
template <int SIZE>
class  LibMeshCubeGenerator : public MeshGenerator
{
public:
    enum { NUM_NODES = (SIZE+1)*(SIZE+1)*(SIZE+1) ,
        NUM_ELEMS = SIZE*SIZE*SIZE };

    virtual void build_mesh()
    {
        mesh = AMP::Mesh::MeshAdapter::shared_ptr ( new AMP::Mesh::LibMeshAdapter () );
        mesh->generateCube ( SIZE );
    }
};


// Class to read in a default exodus file
class  ExodusReaderGenerator : public MeshGenerator
{
public:
    void build_mesh() {
        mesh = AMP::Mesh::MeshAdapter::shared_ptr ( new AMP::Mesh::LibMeshAdapter () );
        mesh->readExodusIIFile ( "clad_1x_1pellet.e" );
    }
};


// Mesh Adapter generator
template <typename ADAPTER>
class   AdapterFromManagerTmpl : public MeshGenerator
{
public:
    static typename AMP::Mesh::MeshManagerTmpl<ADAPTER>::shared_ptr   d_Manager;
    static int                                    d_which;

    virtual void build_mesh() {
        AMP::Mesh::MeshManager::MeshIterator cur = d_Manager->beginMeshes();
        int i = 1;
        while ( cur != d_Manager->endMeshes() ) {
            if ( i == d_which ) {
                mesh = *cur;
                return;
            }
            i++;
            cur++;
        }
    }
};


template <typename ADAPTER>
typename AMP::Mesh::MeshManagerTmpl<ADAPTER>::shared_ptr  AdapterFromManagerTmpl<ADAPTER>::d_Manager;

template <typename ADAPTER>
int  AdapterFromManagerTmpl<ADAPTER>::d_which;

typedef AdapterFromManagerTmpl<AMP::Mesh::MeshAdapter>   AdapterFromManager;


// ThreeElement generator
class   ThreeElementLGenerator : public MeshGenerator
{
public:

    static std::vector<unsigned int> getBndDofIndices() {
        std::vector<unsigned int> bndDofIndices(4);
        bndDofIndices[0] = 0;
        bndDofIndices[1] = 3;
        bndDofIndices[2] = 4;
        bndDofIndices[3] = 7;
        return bndDofIndices;
    }

    static std::vector<std::vector<unsigned int> > getElemNodeMap() {
        std::vector<std::vector<unsigned int> > elemNodeMap(3);

        elemNodeMap[0].resize(8);
        elemNodeMap[1].resize(8);
        elemNodeMap[2].resize(8);

        elemNodeMap[0][0] = 0;
        elemNodeMap[0][1] = 1;
        elemNodeMap[0][2] = 2;
        elemNodeMap[0][3] = 3;
        elemNodeMap[0][4] = 4;
        elemNodeMap[0][5] = 5;
        elemNodeMap[0][6] = 6;
        elemNodeMap[0][7] = 7;

        elemNodeMap[1][0] = 1;
        elemNodeMap[1][1] = 8;
        elemNodeMap[1][2] = 9;
        elemNodeMap[1][3] = 2;
        elemNodeMap[1][4] = 5;
        elemNodeMap[1][5] = 10;
        elemNodeMap[1][6] = 11;
        elemNodeMap[1][7] = 6;

        elemNodeMap[2][0] = 2;
        elemNodeMap[2][1] = 9;
        elemNodeMap[2][2] = 12;
        elemNodeMap[2][3] = 13;
        elemNodeMap[2][4] = 6;
        elemNodeMap[2][5] = 11;
        elemNodeMap[2][6] = 14;
        elemNodeMap[2][7] = 15;

        return elemNodeMap;
    }

    virtual void build_mesh() {
        const unsigned int mesh_dim = 3;
        const unsigned int num_elem = 3;
        const unsigned int num_nodes = 16;

        boost::shared_ptr< ::Mesh > local_mesh(new ::Mesh(mesh_dim));
        local_mesh->reserve_elem(num_elem);
        local_mesh->reserve_nodes(num_nodes);

        local_mesh->add_point(::Point(0.0, 0.0, 0.0), 0);
        local_mesh->add_point(::Point(0.5, 0.0, 0.0), 1);
        local_mesh->add_point(::Point(0.5, 0.5, 0.0), 2);
        local_mesh->add_point(::Point(0.0, 0.5, 0.0), 3);
        local_mesh->add_point(::Point(0.0, 0.0, 0.5), 4);
        local_mesh->add_point(::Point(0.5, 0.0, 0.5), 5);
        local_mesh->add_point(::Point(0.5, 0.5, 0.5), 6);
        local_mesh->add_point(::Point(0.0, 0.5, 0.5), 7);
        local_mesh->add_point(::Point(1.0, 0.0, 0.0), 8);
        local_mesh->add_point(::Point(1.0, 0.5, 0.0), 9);
        local_mesh->add_point(::Point(1.0, 0.0, 0.5), 10);
        local_mesh->add_point(::Point(1.0, 0.5, 0.5), 11);
        local_mesh->add_point(::Point(1.0, 1.0, 0.0), 12);
        local_mesh->add_point(::Point(0.5, 1.0, 0.0), 13);
        local_mesh->add_point(::Point(1.0, 1.0, 0.5), 14);
        local_mesh->add_point(::Point(0.5, 1.0, 0.5), 15);

        std::vector<std::vector<unsigned int> > elemNodeMap = getElemNodeMap();

        for(size_t i = 0; i < elemNodeMap.size(); i++) {
            ::Elem* elem = local_mesh->add_elem(new ::Hex8);
            for(int j = 0; j < 8; j++) {
                elem->set_node(j) = local_mesh->node_ptr(elemNodeMap[i][j]);
            }
        }

        const short int boundaryId = 1;
        std::vector<unsigned int> bndDofIndices = getBndDofIndices(); 

        for(size_t i = 0; i < bndDofIndices.size(); i++) {
            local_mesh->boundary_info->add_node(local_mesh->node_ptr(bndDofIndices[i]), boundaryId);
        }

        local_mesh->prepare_for_use(true);

        mesh = AMP::Mesh::MeshAdapter::shared_ptr ( new AMP::Mesh::MeshAdapter (local_mesh) );
    }


};


//AMP::Mesh::MeshAdapter::shared_ptr  unit_test_common_mesh;

 
}
}

#endif
