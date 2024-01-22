#ifndef included_AMP_Mesh_MeshFactory
#define included_AMP_Mesh_MeshFactory

#include "AMP/utils/FactoryStrategy.hpp"


namespace AMP::IO {
class RestartManager;
}


namespace AMP::Mesh {

class Mesh;
class MeshParameters;


template<class MESH>
std::unique_ptr<MESH> createMeshFromParams( std::shared_ptr<MeshParameters> params )
{
    return std::make_unique<MESH>( params );
}
template<class MESH>
std::unique_ptr<MESH> createMeshFromRestart( int64_t fid, AMP::IO::RestartManager *manager )
{
    return std::make_unique<MESH>( fid, manager );
}


//! Mesh factory class
class MeshFactory
{
public:
    //! get a singleton instance of the factory
    static MeshFactory &getFactory()
    {
        static MeshFactory singletonInstance;
        return singletonInstance;
    }

    //! Create the mesh from the parameters
    static std::shared_ptr<Mesh> create( std::shared_ptr<MeshParameters> parameters );


    //! Create the mesh from the restart file
    static std::shared_ptr<Mesh> create( int64_t fid, AMP::IO::RestartManager *manager );

    //! Register a mesh with the factory
    template<class MESH>
    static void registerMesh( const std::string &name )
    {
        FactoryStrategy<Mesh, std::shared_ptr<MeshParameters>>::registerFactory(
            name, createMeshFromParams<MESH> );
        FactoryStrategy<Mesh, int64_t, AMP::IO::RestartManager *>::registerFactory(
            name, createMeshFromRestart<MESH> );
    }
};


} // namespace AMP::Mesh

#endif
