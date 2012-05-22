#ifndef included_AMP_createLibmeshElements
#define included_AMP_createLibmeshElements

/* AMP files */
#include "ampmesh/Mesh.h"
#include "elem.h"

namespace AMP {
namespace Discretization {


/**
  This is a helper class to create libmesh elements given a MeshIterator.
  It will cache the elements and allow for fast O(log(n)) access to the 
  libmesh element given the MeshElementID. 
*/
class createLibmeshElements
{
public:
    //! Empty constructor
    createLibmeshElements();

    //! De-constructor
    ~createLibmeshElements();

    /**
     *  This function initializes / re-intializes the class to create the libmesh elements
     *  for all of the MeshElements in the given iterator.
     *  \param[in]  iterator  MeshElementIterator containing the elements of interest.
     */
    void reinit( const AMP::Mesh::MeshIterator &iterator );


    /**
     *  This function returns the libmesh element given a MeshElementID
     *  \param[in]  id  MeshElementID for the element of interest
     */
    ::Elem* getElement( const AMP::Mesh::MeshElementID &id );

private:
    std::vector<AMP::Mesh::MeshElementID> d_ids;
    std::vector< ::Elem* > d_elements;
};


}
}

#endif

