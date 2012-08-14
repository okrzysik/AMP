#include "StructuredMeshHelper.h"

namespace AMP {
namespace Mesh {


AMP::Mesh::MeshIterator  StructuredMeshHelper::getXYFaceIterator(AMP::Mesh::Mesh::shared_ptr subChannel, int ghostWidth)
{

    std::multimap<double,AMP::Mesh::MeshElement> xyFace;

    AMP::Mesh::MeshIterator iterator = subChannel->getIterator( AMP::Mesh::Face, ghostWidth );

    for(size_t i=0; i<iterator.size(); ++i ) {
        std::vector<AMP::Mesh::MeshElement> nodes = iterator->getElements(AMP::Mesh::Vertex);
        std::vector<double> center = iterator->centroid();
        bool is_valid = true;
        for (size_t j=0; j<nodes.size(); ++j) {
            std::vector<double> coord = nodes[j].coord();
            if ( !AMP::Utilities::approx_equal(coord[2],center[2], 1e-6) )
                is_valid = false;
        }
        if ( is_valid ) {
            xyFace.insert(std::pair<double,AMP::Mesh::MeshElement>(center[2],*iterator));
        }
        ++iterator;
    }

    boost::shared_ptr<std::vector<AMP::Mesh::MeshElement> > elements( 
        new std::vector<AMP::Mesh::MeshElement>() );
    elements->reserve(xyFace.size());
    for (std::multimap<double,AMP::Mesh::MeshElement>::iterator it=xyFace.begin(); it!=xyFace.end(); ++it)
        elements->push_back( it->second );

    return AMP::Mesh::MultiVectorIterator( elements );
}


AMP::Mesh::MeshIterator  StructuredMeshHelper::getGapFaceIterator(AMP::Mesh::Mesh::shared_ptr subChannel, int ghostWidth)
{
    AMP_ERROR("Not finished");
    return AMP::Mesh::MeshIterator();
}


}
}
