#if defined(USE_AMP_MESH) && defined(USE_LIBMESH)
#ifndef included_AMP_ReadTestMesh
#define included_AMP_ReadTestMesh

#include "boost/shared_ptr.hpp"

#include "mesh.h"

#include "utils/InputDatabase.h"

#include <cstring>

namespace AMP {

  void readTestMesh(boost::shared_ptr<AMP::InputDatabase> mesh_file_db, boost::shared_ptr< ::Mesh > mesh);

  void readTestMesh(std::string mesh_file, boost::shared_ptr< ::Mesh > mesh);

  void readBinaryTestMesh(std::string mesh_file, boost::shared_ptr< ::Mesh> mesh);

}

#endif
#endif
