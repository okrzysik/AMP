#include "ampmesh/MeshParameters.h"
#include "utils/AMP_MPI.h"

namespace AMP {
namespace Mesh {


/********************************************************
 * Constructors                                          *
 ********************************************************/
MeshParameters::MeshParameters() : comm( AMP_COMM_NULL ), MAX_GCW_WIDTH( 1 ) {}
MeshParameters::MeshParameters( const AMP::shared_ptr<AMP::Database> db )
    : d_db( db ), comm( AMP_COMM_NULL ), MAX_GCW_WIDTH( 1 )
{
}


/********************************************************
 * De-constructor                                        *
 ********************************************************/
MeshParameters::~MeshParameters() {}


/********************************************************
 * Set the desired communicator                          *
 ********************************************************/
void MeshParameters::setComm( AMP::AMP_MPI comm_in ) { comm = comm_in; }


/********************************************************
 * Return the database                                   *
 ********************************************************/
AMP::shared_ptr<AMP::Database> MeshParameters::getDatabase() { return d_db; }


} // namespace Mesh
} // namespace AMP
