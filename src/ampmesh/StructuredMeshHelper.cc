#include "StructuredMeshHelper.h"

namespace AMP {
namespace Mesh {


/************************************************************
* Function to return the coordinates of a cube mesh         *
************************************************************/
void StructuredMeshHelper::getXYZCoordinates(AMP::Mesh::Mesh::shared_ptr mesh, 
        std::vector<double>& x_out, std::vector<double>& y_out, std::vector<double>& z_out )
{
    AMP_ASSERT(mesh!=NULL);
    std::set<double> x, y, z;
    if ( mesh.get() != NULL ) {
        AMP::Mesh::MeshIterator it = mesh->getIterator( AMP::Mesh::Vertex, 0 );
        for (size_t i=0; i<it.size(); i++) {
            std::vector<double> coord = it->coord();
            AMP_ASSERT(coord.size()==3);
            x.insert( coord[0] );
            y.insert( coord[1] );
            z.insert( coord[2] );
            ++it;
        }
    }
    mesh->getComm().setGather(x);
    mesh->getComm().setGather(y);
    mesh->getComm().setGather(z);
    x_out.resize(0);
    y_out.resize(0);
    z_out.resize(0);
    x_out.reserve(x.size());
    y_out.reserve(y.size());
    z_out.reserve(z.size());
    x_out.push_back(*(x.begin()));
    y_out.push_back(*(y.begin()));
    z_out.push_back(*(z.begin()));
    for (std::set<double>::iterator it=x.begin(); it!=x.end(); ++it) {
        if ( !Utilities::approx_equal(x_out.back(),*it,1e-12) )
            x_out.push_back(*it);
    }
    for (std::set<double>::iterator it=y.begin(); it!=y.end(); ++it) {
        if ( !Utilities::approx_equal(y_out.back(),*it,1e-12) )
            y_out.push_back(*it);
    }
    for (std::set<double>::iterator it=z.begin(); it!=z.end(); ++it) {
        if ( !Utilities::approx_equal(z_out.back(),*it,1e-12) )
            z_out.push_back(*it);
    }
    size_t Nx = x.size()-1;
    size_t Ny = y.size()-1;
    size_t Nz = z.size()-1;
    AMP_ASSERT(Nx*Ny*Nz==mesh->numGlobalElements(AMP::Mesh::Volume));
}



/************************************************************
* Functions to iterators over particular sets of faces      *
************************************************************/
AMP::Mesh::MeshIterator  StructuredMeshHelper::getXYFaceIterator(
    AMP::Mesh::Mesh::shared_ptr mesh, int gcw )
{
    return getFaceIterator( mesh, gcw, 2 );
}
AMP::Mesh::MeshIterator  StructuredMeshHelper::getXZFaceIterator(
    AMP::Mesh::Mesh::shared_ptr mesh, int gcw )
{
    return getFaceIterator( mesh, gcw, 1 );
}
AMP::Mesh::MeshIterator  StructuredMeshHelper::getYZFaceIterator(
    AMP::Mesh::Mesh::shared_ptr mesh, int gcw )
{
    return getFaceIterator( mesh, gcw, 0 );
}
AMP::Mesh::MeshIterator  StructuredMeshHelper::getFaceIterator(
    AMP::Mesh::Mesh::shared_ptr mesh, int gcw, int direction)
{
    AMP::Mesh::MeshIterator iterator = mesh->getIterator( AMP::Mesh::Face, gcw );
    std::vector<AMP::Mesh::MeshElement> face_list;
    face_list.reserve(iterator.size());
    for(size_t i=0; i<iterator.size(); ++i ) {
        std::vector<AMP::Mesh::MeshElement> nodes = iterator->getElements(AMP::Mesh::Vertex);
        std::vector<double> center = iterator->centroid();
        bool is_valid = true;
        for (size_t j=0; j<nodes.size(); ++j) {
            std::vector<double> coord = nodes[j].coord();
            if ( !AMP::Utilities::approx_equal(coord[direction],center[direction],1e-12) )
                is_valid = false;
        }
        if ( is_valid )
            face_list.push_back(*iterator);
        ++iterator;
    }
    boost::shared_ptr<std::vector<AMP::Mesh::MeshElement> > elements( 
        new std::vector<AMP::Mesh::MeshElement>() );
    *elements = face_list;
    return AMP::Mesh::MultiVectorIterator( elements );
}


AMP::Mesh::MeshIterator  StructuredMeshHelper::getGapFaceIterator(AMP::Mesh::Mesh::shared_ptr subChannel, int ghostWidth)
{
    AMP_ERROR("Not finished");
    return AMP::Mesh::MeshIterator();
}


}
}
