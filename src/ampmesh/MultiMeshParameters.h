#ifndef included_AMP_MultiMeshParameters
#define included_AMP_MultiMeshParameters

#include "ampmesh/MeshParameters.h"


namespace AMP {
namespace Mesh {


/**
 * \class Mesh
 * \brief A class used to pass the parameters for creating a mesh
 */
class MultiMeshParameters : public MeshParameters
{
public:
    // Constructors
    MultiMeshParameters() : MeshParameters() {}
    MultiMeshParameters( const AMP::shared_ptr<AMP::Database> db ) : MeshParameters( db ) {}

protected:
    //! A vector containing the mesh parameters for the sum meshes
    std::vector<MeshParameters::shared_ptr> params;

    //! A vector containing the number of elements in each submesh
    std::vector<size_t> N_elements;

    //! See AMP::Mesh::Mesh for Mesh class
    friend class Mesh;
    friend class MultiMesh;
};


} // Mesh namespace
} // AMP namespace

#endif
