#if defined( USE_AMP_MESH ) && defined( USE_EXT_LIBMESH )
    #ifndef included_AMP_ReadTestMesh
        #define included_AMP_ReadTestMesh

        #include "AMP/utils/Database.h"
        #include "AMP/utils/UtilityMacros.h"
        #include <memory>

DISABLE_WARNINGS
        #include "libmesh/mesh.h"
ENABLE_WARNINGS

        #include <cstring>

namespace AMP {

void readTestMesh( std::shared_ptr<AMP::Database> mesh_file_db,
                   std::shared_ptr<libMesh::Mesh> mesh );

void readTestMesh( std::string mesh_file, std::shared_ptr<libMesh::Mesh> mesh );

void readBinaryTestMesh( std::string mesh_file, std::shared_ptr<libMesh::Mesh> mesh );
} // namespace AMP

    #endif
#endif
