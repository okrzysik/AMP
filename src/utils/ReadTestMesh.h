#if defined(USE_AMP_MESH) && defined(USE_EXT_LIBMESH)
#ifndef included_AMP_ReadTestMesh
#define included_AMP_ReadTestMesh

#include "utils/shared_ptr.h"

#include "libmesh/mesh.h"

#include "utils/InputDatabase.h"

#include <cstring>

namespace AMP {

  void readTestMesh(AMP::shared_ptr<AMP::InputDatabase> mesh_file_db, AMP::shared_ptr< ::Mesh > mesh);

  void readTestMesh(std::string mesh_file, AMP::shared_ptr< ::Mesh > mesh);

  void readBinaryTestMesh(std::string mesh_file, AMP::shared_ptr< ::Mesh> mesh);

}

#endif
#endif
