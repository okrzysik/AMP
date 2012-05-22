#include "createLibmeshElements.h"
#include "utils/Utilities.h"

#include "face_quad4.h"
#include "cell_hex8.h"
#include "node.h"


namespace AMP {
namespace Discretization {


// Default constuctor
createLibmeshElements::createLibmeshElements() 
{
}


// De-constuctor
createLibmeshElements::~createLibmeshElements() 
{
    reinit( AMP::Mesh::MeshIterator() );
}


// Re-initialize the class
void createLibmeshElements::reinit( const AMP::Mesh::MeshIterator &iterator_in )
{
    // Destroy the existing libmesh elements
    for (size_t i=0; i<d_elements.size(); i++) {
        for(size_t j=0; j<d_elements[i]->n_nodes(); j++) {
            delete d_elements[i]->get_node(j);
            d_elements[i]->set_node(j) = NULL;
        }
        delete d_elements[i];
        d_elements[i] = NULL;
    }
    d_elements.resize(0);
    d_ids.resize(0);
    // Create the new libmesh elements
    AMP::Mesh::MeshIterator iterator = iterator_in.begin();
    d_ids.resize(iterator.size());
    d_elements.resize(iterator.size(),NULL);
    for (size_t i=0; i<iterator.size(); i++) {
        d_ids[i] = iterator->globalID();
        int dim = (int) iterator->elementType();
        std::vector<AMP::Mesh::MeshElement> nodes = iterator->getElements(AMP::Mesh::Vertex);
        if ( dim==3 && nodes.size()==8 ) {
            // We are dealing with a hex8 element
            d_elements[i] = new ::Hex8;
        } else if ( dim==2 && nodes.size()==4 ) {
            // We are dealing with a hex8 element
            d_elements[i] = new ::Quad4;
        } else {
            AMP_ERROR("Unknown element type");
        }
        for(size_t j=0; j<nodes.size(); j++) {
            std::vector<double> pt = nodes[j].coord();
            if ( pt.size()==3 )
                d_elements[i]->set_node(j) = new ::Node(pt[0], pt[1], pt[2], j);
            else if ( pt.size()==2 )
                d_elements[i]->set_node(j) = new ::Node(pt[0], pt[1], 0, j);
            else
                AMP_ERROR("Unsupported physical dimension");
        }
        ++iterator;
    }
    // Sort the ids and elements for fast search
    AMP::Utilities::quicksort( d_ids, d_elements );
}


// Search and return the desired libmesh element
::Elem* createLibmeshElements::getElement( const AMP::Mesh::MeshElementID &id )
{
    size_t index = AMP::Utilities::findfirst( d_ids, id );
    if ( index==d_ids.size() ) { index--; }
    AMP_INSIST(d_ids[index]==id,"Desired element was not found");
    return d_elements[index];
}


}
}


