#include "ampmesh/MeshParameters.h"
#include "utils/AMP_MPI.h"

namespace AMP { 
namespace Mesh {


/********************************************************
* Constructors                                          *
********************************************************/
MeshParameters::MeshParameters( )
{
    d_db = boost::shared_ptr<AMP::Database>();
    MAX_GCW_WIDTH = 1;
    comm = AMP::AMP_MPI(AMP_COMM_NULL);
}
MeshParameters::MeshParameters( const boost::shared_ptr<AMP::Database>  db )
{
    d_db = db;
    MAX_GCW_WIDTH = 1;
    comm = AMP::AMP_MPI(AMP_COMM_NULL);
}


/********************************************************
* De-constructor                                        *
********************************************************/
MeshParameters::~MeshParameters( )
{
}


/********************************************************
* Set the desired communicator                          *
********************************************************/
void MeshParameters::setComm ( AMP::AMP_MPI comm_in )
{
    comm = comm_in;
}


/********************************************************
* Return the database                                   *
********************************************************/
boost::shared_ptr<AMP::Database> MeshParameters::getDatabase ( )
{
    return d_db;
}


} // Mesh namespace
} // AMP namespace

