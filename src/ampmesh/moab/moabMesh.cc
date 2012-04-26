#include "ampmesh/moab/moabMesh.h"
#include "ampmesh/MeshElementVectorIterator.h"
#include "ampmesh/MultiIterator.h"
#include "utils/MemoryDatabase.h"
#include "utils/AMPManager.h"
#include "utils/Utilities.h"

#ifdef USE_AMP_VECTORS
    #include "vectors/Vector.h"
    #include "vectors/Variable.h"
    #include "vectors/VectorBuilder.h"
#endif
#ifdef USE_AMP_DISCRETIZATION
    #include "discretization/DOF_Manager.h"
    #include "discretization/simpleDOF_Manager.h"
#endif


// moab includes
#include "moab/Core.hpp"


namespace AMP {
namespace Mesh {


/********************************************************
* Constructors                                          *
********************************************************/
moabMesh::moabMesh( const MeshParameters::shared_ptr &params_in ):
    Mesh(params_in)
{
    d_core = boost::shared_ptr<moab::Core>( new moab::Core() );
    AMP_ERROR("Not finished");
}


/********************************************************
* De-constructor                                        *
********************************************************/
moabMesh::~moabMesh()
{
}


/********************************************************
* Function to copy the mesh                             *
********************************************************/
Mesh moabMesh::copy() const
{
    return moabMesh(*this);
}


/********************************************************
* Function to initialize the moabMesh object             *
********************************************************/
void moabMesh::initialize()
{
    AMP_ERROR("Not finished");
}


/********************************************************
* Function to estimate the mesh size                    *
********************************************************/
size_t moabMesh::estimateMeshSize( const MeshParameters::shared_ptr &params )
{
    AMP_ERROR("Not finished");
}


/********************************************************
* Return the number of elements                         *
********************************************************/
size_t moabMesh::numLocalElements( const GeomType type ) const
{
    AMP_ERROR("Not finished");
}
size_t moabMesh::numGlobalElements( const GeomType type ) const
{
    AMP_ERROR("Not finished");
}
size_t moabMesh::numGhostElements( const GeomType type, int gcw ) const
{
    AMP_ERROR("Not finished");
}


/********************************************************
* Return an iterator over the given geometric type      *
********************************************************/
MeshIterator moabMesh::getIterator( const GeomType type, const int gcw ) const
{
    AMP_ERROR("Not finished");
}


/********************************************************
* Return an iterator over the given boundary ids        *
* Note: we have not programmed this for ghosts yet      *
********************************************************/
MeshIterator moabMesh::getSurfaceIterator ( const GeomType type, const int gcw ) const
{
    AMP_ERROR("Not finished");
    return MeshIterator();
}


/********************************************************
* Return an iterator over the given boundary ids        *
* Note: we have not programmed this for ghosts yet      *
********************************************************/
std::vector<int> moabMesh::getBoundaryIDs ( ) const
{
    AMP_ERROR("Not finished");
}
MeshIterator moabMesh::getBoundaryIDIterator ( const GeomType type, const int id, const int gcw ) const
{
    AMP_ERROR("Not finished");
}


/********************************************************
* Displace a mesh                                       *
********************************************************/
void moabMesh::displaceMesh( std::vector<double> x_in )
{
    AMP_ERROR("Not finished");
}
#ifdef USE_AMP_VECTORS
void moabMesh::displaceMesh( const AMP::LinearAlgebra::Vector::const_shared_ptr x )
{
    #ifdef USE_AMP_DISCRETIZATION
        AMP_ERROR("Not finished");
    #else
        AMP_ERROR("displaceMesh requires DISCRETIZATION");
    #endif
}
#endif


} // Mesh namespace
} // AMP namespace

