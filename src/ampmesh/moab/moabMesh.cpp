#include "AMP/ampmesh/moab/moabMesh.h"
#include "AMP/ampmesh/MeshElementVectorIterator.h"
#include "AMP/ampmesh/MultiIterator.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"

#ifdef USE_AMP_VECTORS
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#endif
#ifdef USE_AMP_DISCRETIZATION
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#endif


// moab includes
#include "moab/Core.hpp"


namespace AMP {
namespace Mesh {


/********************************************************
 * Constructors                                          *
 ********************************************************/
moabMesh::moabMesh( std::shared_ptr<const MeshParameters> params_in ) : Mesh( params_in )
{
    d_core = std::shared_ptr<moab::Core>( new moab::Core() );
    AMP_ERROR( "Not finished" );
}


/********************************************************
 * Destructor                                            *
 ********************************************************/
moabMesh::~moabMesh() {}


/********************************************************
 * Return the class name                                 *
 ********************************************************/
std::string moabMesh::meshClass() const { return "moabMesh"; }


/********************************************************
 * Function to copy the mesh                             *
 ********************************************************/
Mesh moabMesh::copy() const { return moabMesh( *this ); }


/********************************************************
 * Function to initialize the moabMesh object             *
 ********************************************************/
void moabMesh::initialize() { AMP_ERROR( "Not finished" ); }


/********************************************************
 * Function to estimate the mesh size                    *
 ********************************************************/
size_t moabMesh::estimateMeshSize( std::shared_ptr<const MeshParameters> params )
{
    AMP_ERROR( "Not finished" );
    return 0;
}


/********************************************************
 * Return the number of elements                         *
 ********************************************************/
size_t moabMesh::numLocalElements( const GeomType type ) const
{
    AMP_ERROR( "Not finished" );
    return 0;
}
size_t moabMesh::numGlobalElements( const GeomType type ) const
{
    AMP_ERROR( "Not finished" );
    return 0;
}
size_t moabMesh::numGhostElements( const GeomType type, int gcw ) const
{
    AMP_ERROR( "Not finished" );
    return 0;
}


/********************************************************
 * Return an iterator over the given geometric type      *
 ********************************************************/
MeshIterator moabMesh::getIterator( const GeomType type, const int gcw ) const
{
    AMP_ERROR( "Not finished" );
    return MeshIterator();
}


/********************************************************
 * Return an iterator over the given boundary ids        *
 * Note: we have not programmed this for ghosts yet      *
 ********************************************************/
MeshIterator moabMesh::getSurfaceIterator( const GeomType type, const int gcw ) const
{
    AMP_ERROR( "Not finished" );
    return MeshIterator();
}


/********************************************************
 * Return an iterator over the given boundary ids        *
 * Note: we have not programmed this for ghosts yet      *
 ********************************************************/
std::vector<int> moabMesh::getBoundaryIDs() const
{
    AMP_ERROR( "Not finished" );
    return std::vector<int>();
}
MeshIterator
moabMesh::getBoundaryIDIterator( const GeomType type, const int id, const int gcw ) const
{
    AMP_ERROR( "Not finished" );
    return MeshIterator();
}


/********************************************************
 * Displace a mesh                                       *
 ********************************************************/
void moabMesh::displaceMesh( std::vector<double> x_in ) { AMP_ERROR( "Not finished" ); }
#ifdef USE_AMP_VECTORS
void moabMesh::displaceMesh( const AMP::LinearAlgebra::Vector::const_shared_ptr x )
{
#ifdef USE_AMP_DISCRETIZATION
    AMP_ERROR( "Not finished" );
#else
    AMP_ERROR( "displaceMesh requires DISCRETIZATION" );
#endif
}
#endif


} // namespace Mesh
} // namespace AMP
