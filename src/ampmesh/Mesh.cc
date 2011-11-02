#include "ampmesh/Mesh.h"
#include "utils/Utilities.h"

namespace AMP {
namespace Mesh {


/********************************************************
* Constructors                                          *
********************************************************/
Mesh::Mesh( const MeshParameters::shared_ptr &params_in )
{
    params = params_in;
    GeomDim = null;
    comm = params->comm;
    d_db = params->d_db;
}
Mesh::Mesh( const Mesh::shared_ptr &old_mesh )
{
    AMP_ERROR("Not Implimented Yet");
}


/********************************************************
* De-constructor                                        *
********************************************************/
Mesh::~Mesh()
{
}


/********************************************************
* Assignment operator                                   *
********************************************************/
Mesh Mesh::operator=(const Mesh& rhs)
{
    return rhs.copy();
}
Mesh Mesh::copy() const
{
    return Mesh(*this);
}


/********************************************************
* Functions that aren't implimented for teh base class  *
********************************************************/
boost::shared_ptr<Mesh> Mesh::Subset( MeshIterator::shared_ptr & )
{
    AMP_ERROR("Not Implimented Yet");
    return boost::shared_ptr<Mesh>();
}
boost::shared_ptr<Mesh> Mesh::Subset( Mesh & )
{
    AMP_ERROR("Not Implimented Yet");
    return boost::shared_ptr<Mesh>();
}
MeshIterator Mesh::getIterator( GeomType &, int )
{
    AMP_ERROR("Not Implimented Yet");
    return MeshIterator();
}
MeshIterator Mesh::getSurfaceIterator( GeomType & )
{
    AMP_ERROR("Not Implimented Yet");
    return MeshIterator();
}
MeshIterator Mesh::getIterator( SetOP &, MeshIterator::shared_ptr &, MeshIterator::shared_ptr &)
{
    AMP_ERROR("Not Implimented Yet");
    return MeshIterator();
}
size_t Mesh::numLocalElements( GeomType &type )
{
    AMP_ERROR("Not Implimented Yet");
    return 0;
}
size_t Mesh::numTotalElements( GeomType &type )
{
    AMP_ERROR("Not Implimented Yet");
    return 0;
}


} // Mesh namespace
} // AMP namespace

