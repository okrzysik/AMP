#if defined( USE_AMP_MESH ) && defined( USE_EXT_LIBMESH )
#ifndef included_AMP_ReadTestMesh
#define included_AMP_ReadTestMesh

#include "utils/InputDatabase.h"
#include "utils/UtilityMacros.h"
#include "utils/shared_ptr.h"

DISABLE_WARNINGS
#include "libmesh/mesh.h"
ENABLE_WARNINGS

#include <cstring>

namespace AMP {

void readTestMesh( AMP::shared_ptr<AMP::InputDatabase> mesh_file_db, AMP::shared_ptr<::Mesh> mesh );

void readTestMesh( std::string mesh_file, AMP::shared_ptr<::Mesh> mesh );

void readBinaryTestMesh( std::string mesh_file, AMP::shared_ptr<::Mesh> mesh );
} // namespace AMP

#endif
#endif
