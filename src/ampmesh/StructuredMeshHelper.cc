#include "StructuredMeshHelper.h"

namespace AMP {
namespace Mesh {



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
