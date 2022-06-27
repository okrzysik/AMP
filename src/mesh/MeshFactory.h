#ifndef included_AMP_Mesh_MeshFactory
#define included_AMP_Mesh_MeshFactory

#include "AMP/utils/FactoryStrategy.hpp"


namespace AMP::Mesh {

class Mesh;
class MeshParameters;


//! Mesh factory class
class MeshFactory : public FactoryStrategy<Mesh, std::shared_ptr<MeshParameters>>
{
public:
    static std::shared_ptr<Mesh> create( std::shared_ptr<MeshParameters> parameters );
};


} // namespace AMP::Mesh

#endif
